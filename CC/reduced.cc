#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#define DEBUG 0
using namespace std;

extern "C" void read_mmaped_memory(char * memory, long size);
static const unsigned long int PAGESIZE=getpagesize();
const int APP_DEMON=100;    // Tag for communication app to daemon
const int DEMON_APP=101;    // Reverse of the above
const int MAX_FILENAME_LENGTH=500; 
const int MAX_CR_ENV_BUF=2*sizeof(int)+sizeof(long int)+MAX_FILENAME_LENGTH;
const int MAX_BUF_SIZE=134217728; //    Buffer cannot be larger than 128 MB
const int NORM_BUF_SIZE=8388608;
const int DUMMY_POSITION=4096;   //  We write from this position in order to
// avoid position==0. Zero indicates "no entry" in the data arrays. The number
// 4096 fits most filesystems blocks sizes.

int g_curr_buf_size=NORM_BUF_SIZE; // initial buffer size 8 MB
static char* g_buffer=new char[NORM_BUF_SIZE];
string g_str_rank;

MPI_Comm C_COMM=0,g_Intercomm;
int g_io_size=0, g_my_rank=-1;
//
// FILES data structures:
//
static const int max_files=50; 
char * g_filepaths[max_files+1]; // Array containing pointers to files' paths
long g_filetable[max_files+1][3];// 0: file local descriptor; 1: max record size
                                 // 2: position of EOF (i.e. file size)
static int tj_file_numbers[max_files]={0}; //for process managing opening files
long * rec_info0[max_files+1];        // entries with addresses to rec info
#ifdef DEBUG
long * rec_info1[max_files+1];        // array of addresses of arrays of 
                                     // records lengths
#endif
int    rec_lens[max_files+1];        // lengths of data structures in rec_info's

/////////////////////////////////////////////////////////////////////

bool check_buffer(int size) {
if (size>MAX_BUF_SIZE) {cerr << "Too big chunk of data" << endl; abort();}
if (size>g_curr_buf_size) {
  delete [] g_buffer;
  g_buffer=new char[size];
  g_curr_buf_size=size;
  return true;
  }
if (size<=NORM_BUF_SIZE && g_curr_buf_size>NORM_BUF_SIZE) {
  delete [] g_buffer;
  g_buffer=new char[NORM_BUF_SIZE];
  g_curr_buf_size=NORM_BUF_SIZE;
  return true;
  }
return false;
}
// This subroutine inreases, if needed, the length of array containing
// record info
bool check_array_size(int fileno,long add) {
  const int m=sizeof(long);
  if (add<rec_lens[fileno]) return false;
  int new_len=rec_lens[fileno]*1.2;  // increase by 1.2:
  if (add>=new_len) new_len=add+1;
  if (new_len>512*1024)
       cout << "Warning, over 0.5 million of record entries" << endl;
  if (new_len>1024*1024*16) 
      {cerr << "Error, over 16 millions of record entries" << endl; abort(); }
  new_len=(new_len+PAGESIZE/m-1) & ~(PAGESIZE/m-1); // round to pagesize
  if (add>=new_len) {cerr << "Error 4345hd" << endl; abort();}
  int old_len=rec_lens[fileno];
  rec_lens[fileno]=new_len;
  void* new_ad=mremap(rec_info0[fileno],old_len*m,new_len*m,MREMAP_MAYMOVE);
  if (new_ad==MAP_FAILED) 
       {perror("mremap failed in database"); abort();}
  rec_info0[fileno]=(long*) new_ad;
  memset(rec_info0[fileno]+old_len,0,m*(new_len-old_len));
#ifdef DEBUG
        new_ad=mremap(rec_info1[fileno],old_len*m,new_len*m,MREMAP_MAYMOVE);
  if (new_ad==MAP_FAILED) 
       {perror("mremap failed in database"); abort();}
  rec_info1[fileno]=(long*) new_ad;
  memset(rec_info1[fileno]+old_len,0,m*(new_len-old_len));
#endif
  return true;
  }

#ifdef __cplusplus
        extern "C"{
#endif

/////////////////////////////////////////////////////////////////////////////

int file_no_create() {
for (int i=0;i<max_files;i++) {
  if (tj_file_numbers[i]==0) {
    tj_file_numbers[i]=1;
    return i;
    }
  }
printf("Number of files exceeded\n");
abort();
}

/////////////////////////////////////////////////////////////////////////////

void file_no_delete(int fileno) {
if (tj_file_numbers[fileno]!=1) {
  printf("Error in file_no_delete\n");
  abort();
  }
tj_file_numbers[fileno]=0;
}

/////////////////////////////////////////////////////////////////////////////

void fafcreatem_(char *filename,long int *fileNumber, long int * max_rec_size) {
  if (g_my_rank) {
    cerr << "This subroutine is for master only" << endl; abort();
    }
  int command=0; // open
  char *buf=new char[MAX_CR_ENV_BUF];
  *((int*)(buf))=command;
  char *file=buf+2*sizeof(int)+sizeof(long int);
  int buflen;
  if (*max_rec_size!=-1) { 
    *fileNumber=file_no_create();
    buflen=2*sizeof(int)+sizeof(long int)+strlen(filename)+1;
    if (buflen>MAX_CR_ENV_BUF) {cerr << "Too long filename" << endl; abort();}
    *((int*)(buf+sizeof(int)))=*fileNumber;
    *((long int*)(buf+2*sizeof(int)))=*max_rec_size;
    strcpy(file,filename); 
  } else {
    buflen=2*sizeof(int)+sizeof(long int)+1;
    *((int*)(buf+sizeof(int)))=-1;
    *((long int*)(buf+2*sizeof(int)))=-1;
    *file=0;
  }
  if (*max_rec_size!=-1) { // Otherwise we have dummy open
    for (int id=0;id<g_io_size;++id) {
      MPI_Ssend(buf,buflen,MPI_BYTE,id,APP_DEMON,g_Intercomm);
      }
    }
  MPI_Bcast(buf,MAX_CR_ENV_BUF,MPI_BYTE,0,C_COMM);
  delete [] buf;
}
void fafcreates_(long int *fileNumber) {
  if (!g_my_rank) {
    cerr << "This subroutine is for slave only" << endl; abort();
    }
  char *buf=new char[MAX_CR_ENV_BUF];
  MPI_Bcast(buf,MAX_CR_ENV_BUF,MPI_BYTE,0,C_COMM);
  int command=*((int*)(buf));
  *fileNumber=*((int*)(buf+sizeof(int)));
  long int max_rec_size=*((long int*)(buf+2*sizeof(int)));
  char *filename=buf+2*sizeof(int)+sizeof(long int);
  if (command!=0) {cerr << "Error in slave opener" << endl; abort();}
  delete [] buf;
}
/////////////////////////////////////////////////////////////////////////////

void fafclosem_(const long int *fileID, const long int *clos, long int *fileStatus) {
  if (g_my_rank) {
    cerr << "This subroutine is for master only" << endl; abort();
    }
   int command=1; // close
   int buflen=2*sizeof(int);
   char *buf=new char[buflen];
   *((int*)(buf))=command;
   *((int*)(buf+sizeof(int)))=*fileID;
   if (*clos!=-1) { // if *clos==-1 do not close - dummy subroutine
     if (*fileID<0) {cerr << "Wrong file number supplied" << endl; abort();}
     file_no_delete(*fileID);
     for (int id=0;id<g_io_size;++id) {
       MPI_Ssend(buf,buflen,MPI_BYTE,id,APP_DEMON,g_Intercomm);
       }
     }
   else *((int*)(buf+sizeof(int)))=-1; // fileID==-1
   MPI_Bcast(buf,2*sizeof(int),MPI_BYTE,0,C_COMM);
   delete [] buf;
   *fileStatus=0;
}
void fafcloses_() {
  if (!g_my_rank) {
    cerr << "This subroutine is for slave only" << endl; abort();
    }
   char *buf=new char[2*sizeof(int)];
   MPI_Bcast(buf,2*sizeof(int),MPI_BYTE,0,C_COMM);
   int command=*((int*)(buf));
   int fileID=*((int*)(buf+sizeof(int)));
   if (command!=1) {
     cerr << "This call does not match master's close" << endl;
     sleep(100);
     abort();}
   delete [] buf;
}

/////////////////////////////////////////////////////////////////////////////

void xxafinit_(int * my_rank, MPI_Comm * G_COMM, MPI_Comm * I_COMM, int * io_size) {
// Just transferring the communicators and related values
  g_my_rank = *my_rank;
  C_COMM = * G_COMM;
  g_Intercomm = * I_COMM;
  g_io_size = *io_size;
  }

/////////////////////////////////////////////////////////////////////////////

void afterminate_(long int *exitCode){
  const long int zero=0;
  const long int dumm=-1;
  long int stat;
  int command=4; // finish
   if (!g_my_rank) {        // g_my_rank is a rank in MPI_COMM_WORLD on slaves 
                            // and master, but intercomm rank on I/O daemons
     for (long int i=0;i<max_files;i++) {
       if (tj_file_numbers[i]!=0) fafclosem_(&i,&zero,&stat);
       else                       fafclosem_(&i,&dumm,&stat);
     }
     for (int id=0;id<g_io_size;++id) {
       MPI_Ssend(&command,sizeof(int),MPI_BYTE,id,APP_DEMON,g_Intercomm);
     }
   } else
     for (long int i=0;i<max_files;i++) fafcloses_();
   MPI_Comm_free(&g_Intercomm); 
   MPI_Comm_free(&C_COMM); 
}

/////////////////////////////////////////////////////////////////////////////

void fafwrite_(long int * fd, char * data, long int * unitsize, long int * size,
               long int * dummy, long int * record, long int * info){
  int command=2; // write
  int id=*record%g_io_size;
  int buflen=2*sizeof(int)+2*sizeof(long int);
  char *buf=new char[buflen];
  *((int*)(buf))=command;
  *((int*)(buf+sizeof(int)))=*fd;
  *((long int*)(buf+2*sizeof(int)))=*record;
  *((long int*)(buf+2*sizeof(int)+sizeof(long)))=*unitsize**size;
  MPI_Ssend(buf,buflen,MPI_BYTE,id,APP_DEMON,g_Intercomm);
  MPI_Ssend(data,*unitsize*(*size),MPI_BYTE,id,APP_DEMON,g_Intercomm);
  delete [] buf;
}

/////////////////////////////////////////////////////////////////////////////

void fafread_(long int * fd, char * data, long int * unitsize, long int * size,
               long int * dummy, long int * record, long int * info){
  int command=3; // read
  int id=*record%g_io_size;
  MPI_Status status;
  int buflen=3*sizeof(int)+sizeof(long int);
  int loc_size;
  char *buf=new char[buflen];
  *((int*)(buf))=command;
  *((int*)(buf+sizeof(int)))=*fd;
  *((long int*)(buf+2*sizeof(int)))=*record;
  *((int*)(buf+sizeof(long int)+2*sizeof(int)))=*unitsize**size;
  MPI_Ssend(buf,buflen,MPI_BYTE,id,APP_DEMON,g_Intercomm);
  MPI_Recv(data,*unitsize**size,MPI_BYTE,id,DEMON_APP,g_Intercomm,&status);
  MPI_Get_count(&status,MPI_BYTE,&loc_size);
  *info=loc_size;
  delete [] buf;
}

/////////////////////////////////////////////////////////////////////////////
void afdx_() {
  int command;
  nice(1);
  for (int i=0;i<max_files;++i) {
    g_filetable[i][0]=-1;
    g_filetable[i][1]=-1;
    g_filetable[i][2]=-1;
    }
  const int ENV_SIZE=100+MAX_FILENAME_LENGTH;
  MPI_Status status;
// Create string-line suffix based on numerical value of my rank
  char * tmp=new char[200];
  memset(tmp,0,200);
  sprintf(tmp,"%04d",g_my_rank);
  g_str_rank=tmp;
  delete [] tmp;
  tmp=0;
//
  char *filename=0;
  for (;;) {
    int command=-1,flag=0;
    void*l_b;
    off_t error;
    MPI_Recv(g_buffer,MAX_BUF_SIZE,MPI_BYTE,MPI_ANY_SOURCE,APP_DEMON,
             g_Intercomm,&status);
    int peer=status.MPI_SOURCE;
    int env_size=0;
    MPI_Get_count(&status,MPI_BYTE,&env_size);
    if (env_size<4) {cerr << "Envelope too small" << endl; abort();}
    string myfilename("");
    int fd;
    size_t size,size1;
    int lrec;
    off_t rec;
    ssize_t read1;
    command=*((int*)g_buffer);
    int rounded_size;
    int offset;
    int filenumber;
    long int max_rec_size;
    long trsize=0;
    switch (command) {
    case 0:            // open
      filenumber=*((int*)(g_buffer+sizeof(int)));
      if (filenumber>=max_files) {cerr << "I/O error hfs333" << endl; abort();}
      max_rec_size=*((long int*)(g_buffer+2*sizeof(int)));
      filename=(g_buffer+2*sizeof(int)+sizeof(long int));
      myfilename=string(filename)+string("_")+g_str_rank;
      if (g_filetable[filenumber][1]!=-1) 
                                {cerr << "I/O error 79taab" << endl; abort();}
      if (max_rec_size<256) {
      cout << "Array Files I/O - near-null record size: " << filename
           << "System too small for efficient parallel execution"  << endl;
        abort();
        }
//      max_rec_size=(max_rec_size+PAGESIZE-1) & ~(PAGESIZE-1);
      g_filetable[filenumber][1]=max_rec_size;
                                       // we start writing from this position
      g_filetable[filenumber][2]=DUMMY_POSITION; 
                                       // - avoid record with pos. no. 0
      rec_lens[filenumber]=1024*8; // initial array has 64 KB
      rec_info0[filenumber]=(long*)mmap(0,rec_lens[filenumber]*sizeof(long),
                           PROT_READ|PROT_WRITE|PROT_NONE,MAP_PRIVATE|MAP_ANONYMOUS,0,0);
      if (rec_info0[filenumber]==MAP_FAILED)
                         {perror("database mmap error"); abort();}
      memset(rec_info0[filenumber],0,rec_lens[filenumber]*sizeof(long));
#ifdef DEBUG
      rec_info1[filenumber]=(long*)mmap(0,rec_lens[filenumber]*sizeof(long),
                           PROT_READ|PROT_WRITE|PROT_NONE,MAP_PRIVATE|MAP_ANONYMOUS,0,0);
      if (rec_info1[filenumber]==MAP_FAILED)
                         {perror("database mmap error"); abort();}
      memset(rec_info1[filenumber],0,rec_lens[filenumber]*sizeof(long));
#endif
      filename=0;
      tmp=new char[strlen(myfilename.c_str())+1];
      strcpy(tmp,myfilename.c_str());
      fd=open(tmp,O_RDWR|O_CREAT|O_TRUNC,00666);
      if (fd<0) {perror(""); cerr << "I/O error dfgdfg" << endl; abort();}
      if (g_filetable[filenumber][0]!=-1) 
                                {cerr << "I/O error 39taab" << endl; abort();}
      g_filetable[filenumber][0]=fd;
      g_filepaths[filenumber]=tmp;
      tmp=0;
      read1=write(fd,g_buffer,DUMMY_POSITION); // write some thrash, in order to move curr pos.
      if (read1!=DUMMY_POSITION) {perror("Write dummy error");abort();}
      break;
    case 1:            // close
      filenumber=*((int*)(g_buffer+sizeof(int)));
      if (g_filetable[filenumber][0]==-1)
          {cerr << "I/O error bhh7hb " << filenumber << endl; abort(); abort();}
      if (close(g_filetable[filenumber][0])) 
                         {perror("Close error 3454"); abort();}
      if (unlink(g_filepaths[filenumber])) 
                          {perror("unlink error sdhgd6"); abort();};
      delete [] g_filepaths[filenumber];
      g_filepaths[filenumber]=0;
      if (munmap(rec_info0[filenumber],rec_lens[filenumber]*8))
                       {perror("database munmap error"); abort();}
      rec_info0[filenumber]=0;
#ifdef DEBUG
      if (munmap(rec_info1[filenumber],rec_lens[filenumber]*8))
                       {perror("database munmap1 error"); abort();}
      rec_info1[filenumber]=0;
#endif
      rec_lens[filenumber]=0;
      g_filetable[filenumber][0]=-1;
      g_filetable[filenumber][1]=-1;
      g_filetable[filenumber][2]=-1;
    break;
    case 2:            // write
      filenumber=*((int*)(g_buffer+sizeof(int)));
      lrec=g_filetable[filenumber][1];
      if (lrec<=0) {cout << g_filepaths[filenumber] << " " << filenumber << endl; abort();}
      rec=*((long int*)(g_buffer+2*sizeof(int)));
      trsize=*((long int*)(g_buffer+2*sizeof(int)+sizeof(long)));
      MPI_Probe(peer,APP_DEMON,g_Intercomm,&status);
      int intsize;
      MPI_Get_count(&status,MPI_BYTE,&intsize);
      if (intsize!=trsize) {cerr << "Unexpected condition w" << endl; abort();}
      size=intsize;
      check_buffer(size);
      if (size>lrec) {cerr << "I/O error xxsdhgd6: " << g_filepaths[filenumber] << " " << size << " " << lrec << endl; abort();}
      if (g_filetable[filenumber][0]<0) 
         {cerr << "I/O error xxsdvvd6" << endl; abort();}
      MPI_Recv(g_buffer,size,MPI_BYTE,peer,APP_DEMON,g_Intercomm,&status);
      check_array_size(filenumber,rec/g_io_size);
      if (rec_info0[filenumber][rec/g_io_size]==0) {
        rec_info0[filenumber][rec/g_io_size]=g_filetable[filenumber][2];
        g_filetable[filenumber][2]+=size;
#ifdef DEBUG
        rec_info1[filenumber][rec/g_io_size]=size;
#endif
        read1=write(g_filetable[filenumber][0],g_buffer,size);
        if (read1<0) {perror("Write error"); abort();}
        else if (read1!=size) {
          int prev=read1;
          read1=write(g_filetable[filenumber][0],g_buffer+prev,size-prev);
          if (read1!=(size-prev)) {perror("Write error1");abort();}
          }
        }
      else {
#ifdef DEBUG
        if (rec_info1[filenumber][rec/g_io_size]!=size) {
          cerr << "Attempting to overwrite the record with different size" 
               << endl;
          abort();
          }
#endif
        long myrec=rec_info0[filenumber][rec/g_io_size];
        read1=pwrite(g_filetable[filenumber][0],g_buffer,size,myrec);
        if (read1<0) {perror("Write error");abort();}
        else if (read1!=size) {
          int prev=read1;
          read1=pwrite(g_filetable[filenumber][0],g_buffer+prev,size-prev,
                       myrec+prev);
          if (read1!=(size-prev)) {perror("Write error1");abort();}
          }
        }
    break;
    case 3:            // read
      filenumber=*((int*)(g_buffer+sizeof(int)));
      lrec=g_filetable[filenumber][1];
      rec=*((long int*)(g_buffer+2*sizeof(int)));
      size=*((int*)(g_buffer+2*sizeof(int)+sizeof(long int)));
//    check_buffer(1); // mmap will be used
      if (size>lrec) {cerr << "I/O error zxsdhgd6" << endl; abort();}
      if (g_filetable[filenumber][0]<0) 
         {cerr << "I/O error qxsdvvd6" << endl; abort();}
      if (rec/g_io_size>rec_lens[filenumber]) 
      {cerr << "Attempt to read a record beyond end of file" << endl; abort();}
      if (rec_info0[filenumber][rec/g_io_size]==0) {
        cerr << "Non-existing record" << endl;
        abort();
        }
      else {
#ifdef DEBUG
        if (rec_info1[filenumber][rec/g_io_size]!=size) {
          cerr << "Attempting to read a record with different size "
               << rec_info1[filenumber][rec/g_io_size] << " "
               << size << g_filepaths[filenumber] << " " 
               << filenumber << " "
               << rec << endl;
          abort();
          }
#endif
        off_t myrec=rec_info0[filenumber][rec/g_io_size];
        off_t myrec1=myrec & ~(PAGESIZE-1);
        size1=size+myrec-myrec1;
        l_b=mmap(0,size1+256,PROT_READ|PROT_WRITE,MAP_PRIVATE,
                 g_filetable[filenumber][0],myrec1);
        if (l_b==MAP_FAILED) {perror("Map failed"); abort();}
        MPI_Ssend((char*)l_b+(myrec-myrec1),size,MPI_BYTE,peer,DEMON_APP,g_Intercomm);
        if (munmap(l_b,size1+256)) {perror("Unmap failed"); abort();}
        }
    break;
    case 4:
    goto exit_loop;
    default:
      cerr << "I/O error 33sdvvd6" << endl; abort();
    }
    }
  exit_loop:
  return;
  }
void create_c_communicators_(long int * color) {
// color == 1 Master + slave
// // color == 2 I/O
  int my_rank,size,info;
  long int all_colors[5000]={0};
  int idd=-1;
  info=MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  info=MPI_Comm_size(MPI_COMM_WORLD,&size);
  info=MPI_Allgather(color,1,MPI_LONG,all_colors,1,MPI_LONG,MPI_COMM_WORLD);
  info=MPI_Comm_split(MPI_COMM_WORLD,*color,my_rank,&C_COMM);
  for (int i=0;i<size;++i) if (all_colors[i]==2) { idd=i; break; }
  if (*color==1)
    info=MPI_Intercomm_create(C_COMM,0,MPI_COMM_WORLD,idd,46534,&g_Intercomm);
  else if (*color==2)
    info=MPI_Intercomm_create(C_COMM,0,MPI_COMM_WORLD,0  ,46534,&g_Intercomm);
  else abort();
  info=MPI_Comm_rank(g_Intercomm,&g_my_rank);
  if (*color==1)
    info=MPI_Comm_remote_size(g_Intercomm,&g_io_size);
  else if (*color==2)
    info=MPI_Comm_size(g_Intercomm,&g_io_size);
  }
void remove_c_communicators_(long int * color) {
  MPI_Comm_free(&g_Intercomm);
  MPI_Comm_free(&C_COMM);
  }
#ifdef __cplusplus
}
#endif
