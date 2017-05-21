#include <pvm3.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#define DEBUG 0
using namespace std;

extern "C" void read_mmaped_memory(char * memory, long size);
static const unsigned long int PAGESIZE=getpagesize();
const int APP_DEMON=10100;    // Tag for communication app to daemon
const int DEMON_APP=10101;    // Reverse of the above
const int WR_DATA=10200;    // Tag for communication app to daemon
const int RD_DATA=10201;    // Reverse of the above
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

int g_io_size=0, g_my_rank=-1,master_id,nslv=0;
int* dmnid = NULL;
int* slvid = NULL;
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

void fafcreatem_(char *filename, long int *fileNumber, long int * max_rec_size) {
  if (g_my_rank) {
    cerr << "This subroutine is for master only" << endl; abort();
    }
  int command=0; // open
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&command,1,1);
  pvm_pklong(max_rec_size,1,1);
  if (*max_rec_size!=-1) { // Otherwise we have dummy open, no message to IO d.
    *fileNumber=file_no_create();
    pvm_pklong(fileNumber,1,1);
    pvm_pkstr(filename);
    for (int id=0;id<g_io_size;++id) {
        pvm_send(dmnid[id],APP_DEMON);
      }
    }
  for (int id=0;id<nslv;++id) {
      pvm_send(slvid[id],APP_DEMON);
    }
}
void fafcreates_(long int *fileNumber) {
  if (!g_my_rank) {
    cerr << "This subroutine is for slave only" << endl; abort();
    }
  int command;
  long int max_rec_size;
  pvm_recv(master_id,APP_DEMON);
  pvm_upkint(&command,1,1);
  pvm_upklong(&max_rec_size,1,1);
  if (command!=0) {cerr << "Error in slave opener" << endl; abort();}
}
/////////////////////////////////////////////////////////////////////////////

void fafclosem_(long int *fileID, long int *clos, long int *fileStatus) {
  if (g_my_rank) {
    cerr << "This subroutine is for master only" << endl; abort();
    }
   int command=1; // close
   pvm_initsend(PvmDataRaw);
   pvm_pkint(&command,1,1);
   pvm_pklong(fileID,1,1);
   pvm_pklong(clos,1,1);
   if (*clos!=-1) { // if *clos==-1 do not close - dummy subroutine
     if (*fileID<0) {cerr << "Wrong file number supplied" << endl; abort();}
     file_no_delete(*fileID);
     for (int id=0;id<g_io_size;++id) {
        pvm_send(dmnid[id],APP_DEMON);
       }
     }
   for (int id=0;id<nslv;++id) {
       pvm_send(slvid[id],APP_DEMON);
    }
   *fileStatus=0;
}
void fafcloses_() {
  int command;
  long fileID,clos;
  if (!g_my_rank) {
    cerr << "This subroutine is for slave only" << endl; abort();
    }
   pvm_recv(master_id,APP_DEMON);
   pvm_upkint(&command,1,1);
   pvm_upklong(&fileID,1,1);
   pvm_upklong(&clos,1,1);
   if (command!=1) {
     cerr << "This call does not match master's close" << endl;
     sleep(100);
     abort();}
}

/////////////////////////////////////////////////////////////////////////////

void afinit_(int * my_rank, int * io_size, int * f_dmnid, int * mstrid, int * numsl, int * f_slvid) {
// Just transferring rank,no and tid of IO daemons and slaves
  g_my_rank = *my_rank;
  g_io_size = *io_size;
  dmnid = new int[g_io_size];
  for (int i=0; i<g_io_size;i++){
     dmnid[i] = f_dmnid[i];
  }
  master_id = *mstrid;
  nslv = *numsl;
  slvid = new int[nslv];
  for (int i=0; i<nslv;i++){
     slvid[i] = f_slvid[i+1];
  }
  }

/////////////////////////////////////////////////////////////////////////////

void afterminate_(long int *exitCode){
  long int zero=0;
  long int dumm=-1;
  long int stat;
  int command=4; // finish
   if (!g_my_rank) {        // g_my_rank is a rank in MPI_COMM_WORLD on slaves 
     for (long int i=0;i<max_files;i++) {
       if (tj_file_numbers[i]!=0) fafclosem_(&i,&zero,&stat);
       else                       fafclosem_(&i,&dumm,&stat);
     }
     pvm_initsend(PvmDataRaw);
     pvm_pkint(&command,1,1);
     for (int id=0;id<g_io_size;++id) {
        pvm_send(dmnid[id],APP_DEMON);
     }
   } else
     for (long int i=0;i<max_files;i++) fafcloses_();
   delete [] dmnid; 
   delete [] slvid; 
}

/////////////////////////////////////////////////////////////////////////////

void fafwrite_(long int * fd, char * data, long int * unitsize, long int * size,
               long int * dummy, long int * record, long int * info){
  int command=2; // write
  int id=*record%g_io_size;
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&command,1,1);
  pvm_pklong(fd,1,1);
  pvm_pklong(record,1,1);
  long int bytes = *unitsize * *size;
  pvm_pklong(&bytes,1,1);
  pvm_send(dmnid[id],APP_DEMON);
  pvm_psend(dmnid[id],WR_DATA,data,bytes,PVM_BYTE);
}

/////////////////////////////////////////////////////////////////////////////

void fafread_(long int * fd, char * data, long int * unitsize, long int * size,
               long int * dummy, long int * record, long int * info){
  int command=3; // read
  int id=*record%g_io_size;
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&command,1,1);
  pvm_pklong(fd,1,1);
  pvm_pklong(record,1,1);
  long int bytes = *unitsize * *size;
  pvm_pklong(&bytes,1,1);
  pvm_send(dmnid[id],APP_DEMON);
  int atid,atag,alen;
  pvm_precv(dmnid[id],RD_DATA,data,bytes,PVM_BYTE,&atid,&atag, &alen);
  *info=alen;
}

/////////////////////////////////////////////////////////////////////////////


void afdx_(int * my_rank, int * io_size) {
  int command,cc;
  g_my_rank = *my_rank;
  g_io_size = *io_size;
  nice(1);
  for (int i=0;i<max_files;++i) {
    g_filetable[i][0]=-1;
    g_filetable[i][1]=-1;
    g_filetable[i][2]=-1;
    }
// Create string-line suffix based on numerical value of my rank
  char * tmp=new char[200];
  memset(tmp,0,200);
  sprintf(tmp,"%04d",g_my_rank);
  g_str_rank=tmp;
  delete [] tmp;
  tmp=0;
//
  char *filename=new char[300];
  for (;;) {
    int command=-1,flag=0;
    void*l_b;
    off_t error;
    int bufid = pvm_recv(-1,APP_DEMON);
    string myfilename("");
    int fd;
    size_t size,size1;
    int lrec;
    off_t rec;
    ssize_t read1;
    pvm_upkint(&command,1,1);
    int rounded_size;
    int offset;
    long filenumber;
    long max_rec_size;
    long trsize=0,trec;
    int bytes,msgtag,source;
    switch (command) {
    case 0:            // open
      pvm_upklong(&max_rec_size,1,1);
      pvm_upklong(&filenumber,1,1);
      if (filenumber>=max_files) {cerr << "I/O error hfs333" << endl; abort();}
      pvm_upkstr(filename);

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
//      filename=0;
      tmp=new char[strlen(myfilename.c_str())+1];
      strcpy(tmp,myfilename.c_str());
      fd=open(tmp,O_RDWR|O_CREAT|O_TRUNC,00666);
      if (fd<0) {perror(""); cerr << "I/O error dfgdfg" << endl; abort();}
      if (g_filetable[filenumber][0]!=-1) 
                                {cerr << "I/O error 39taab" << endl; abort();}
      g_filetable[filenumber][0]=fd;
      g_filepaths[filenumber]=tmp;
      tmp=0;
      read1=write(fd,g_buffer,DUMMY_POSITION); // write some trash, in order to move curr pos.
      if (read1!=DUMMY_POSITION) {perror("Write dummy error");abort();}
      break;
    case 1:            // close
      pvm_upklong(&filenumber,1,1);
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
      pvm_upklong(&filenumber,1,1);
      lrec=g_filetable[filenumber][1];
      if (lrec<=0) {cout << g_filepaths[filenumber] << " i  " << filenumber << endl; abort();}
      pvm_upklong(&trec,1,1);
      pvm_upklong(&trsize,1,1);
      rec=trec;
      size=trsize;
      pvm_bufinfo(bufid,&bytes, &msgtag, &source);
      check_buffer(size);
      if (size>lrec) {cerr << "I/O error xxsdhgd6: " << g_filepaths[filenumber] << " " << size << " " << lrec << endl; abort();}
      if (g_filetable[filenumber][0]<0) 
         {cerr << "I/O error xxsdvvd6" << endl; abort();}
      int atid,atag,alen;
      pvm_precv(source,WR_DATA,g_buffer,size,PVM_BYTE,&atid,&atag,&alen);
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
      pvm_upklong(&filenumber,1,1);
      lrec=g_filetable[filenumber][1];
      pvm_upklong(&trec,1,1);
      pvm_upklong(&trsize,1,1);
      rec=trec;
      size=trsize;
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
        pvm_bufinfo(bufid,&bytes, &msgtag, &source);
        pvm_psend(source,RD_DATA,(char*)l_b+(myrec-myrec1),size,PVM_BYTE);
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
#ifdef __cplusplus
}
#endif
