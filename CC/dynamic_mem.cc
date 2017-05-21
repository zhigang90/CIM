#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>
#include <signal.h>
#include <string.h>
#include <iostream>
#include <sstream>
static const unsigned long int DEBUG=1;
static const unsigned long int HARD_DEBUG=1;
static const unsigned long int PAGESIZE=getpagesize();
static inline void assure(const bool condition, const char* const message) {
if (!condition) return;
printf(message);
abort();
} /*
namespace memory {
void * operator new (size_t amount) {
  char * real_start=::new char[amount+PAGESIZE];
  long unsigned int int_real_start=
      reinterpret_cast<long unsigned int>(real_start);
  unsigned long int int_data=(int_real_start + PAGESIZE-1) & ~(PAGESIZE-1);
  assure((int_data-int_real_start)%8!=0,"Wrong memory boundary\n");
  void * data=reinterpret_cast<void*>(int_data);
  return data;
}
void tomek() {}
} */
//******************************************************************************
//
//
//
//  Chunk definitions:
//
//
//
//******************************************************************************
class chunk {
public:
  inline double* allocate(long int amount);// return c++ address
  long int freespace();                 // return amount of free space
  long int free(double * address, long int amount);// return number of allocations after freeing
  chunk(long int size) : top(0), data_size(size), no_alloc(0)  {
          if (DEBUG) {
            real_start=::new double[size+PAGESIZE/sizeof(double)];
            unsigned long int int_real_start=
                                    reinterpret_cast<unsigned long int>(real_start);
            unsigned long int int_data=(int_real_start + PAGESIZE-1) & ~(PAGESIZE-1);
	    assure((int_data-int_real_start)%8!=0,"Wrong memory boundary\n");
            data=reinterpret_cast<double*>(int_data);
            }
          else {
            real_start=::new double[size];
            data=real_start;
            }
  }
  ~chunk() {delete []real_start;}
  long int show_size() {return data_size;}
private:
  double * real_start;
  double * data;
  long int data_size;
  long int top;
  long int no_alloc;
};
//******************************************************************************
inline double* chunk::allocate(long int amount) {
double* address;
assure(amount>(data_size-top),"Dynamic memory internal chunk error\n");
no_alloc++;
address=&data[top];
top+=amount;
return address;
}
//******************************************************************************
inline long int chunk::free(double * address,long int amount) {
long int index=address-data;
assure(index>=top,"Dynamic memory internal error 1\n");
assure((top-index)!=amount,"Dynamic memory internal error 2\n");
no_alloc--;
top=index;
return no_alloc;
}
//******************************************************************************
inline long int chunk::freespace() {                 // return amount of free space
return data_size-top;
}
//******************************************************************************
//
//
//    The main class handling memory.
//
//
//
//******************************************************************************
class allocator {
public:
  enum {max_alloc=20000, max_chunk=2000, chunk_size=131072, max_marks=100}; //1048576=2**20
inline void getmem(const long int &iamount, long int& iaddress);
inline long int retmem(long int iamount);
inline void mmark();
inline void retmark();
inline allocator();
inline void initializer(long int max, double* addr);
inline void dumpall();
inline long int showfree() const;
inline long int allocated(double *address) const;
inline long int lock(double *address) ;
inline long int unlock(double *address) ;
inline void get_ref(double ** ref) {
  assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
  *ref=refaddress;
  }
// exp:
void * operator new (size_t amount) {
  char * real_start=::new char[amount+PAGESIZE];
  unsigned long int int_real_start=
      reinterpret_cast<unsigned long int>(real_start);
  unsigned long int int_data=(int_real_start + PAGESIZE-1) & ~(PAGESIZE-1);
  assure((int_data-int_real_start)%8!=0,"Wrong memory boundary\n");
  void * data=reinterpret_cast<void*>(int_data);
  return data;
}
//
private:
  long int           allocation_no;
  long int           chunk_no;
  long int           alloc2chunkno   [max_alloc];
  long int           alloc2size      [max_alloc];
  long int           alloc2_ord_size [max_alloc];
  double*       alloc2address   [max_alloc];
  chunk*        chunkno2addr    [max_chunk];
  char          locked          [max_alloc];
  double*       refaddress;
  long int           marks[max_marks];
  char          marker;
  long int           max_amount;
  long int           real_amount;
  long int           amount;
};
//******************************************************************************
inline long int allocator::showfree() const {
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
return max_amount-real_amount;
}
//******************************************************************************
inline void allocator::dumpall() {
printf("allocation_no: %ld\n",allocation_no);
printf("chunk_no:      %ld\n",chunk_no);
printf("marker:        %d\n",marker);
printf("max_amount:    %ld\n",max_amount);
printf("amount:        %ld\n",amount);
printf("real_amount:   %ld\n",real_amount);
long int total=0;
for (long int i=0;i<chunk_no;i++) total+=chunkno2addr[i]->freespace();
printf("Space wasted:  %ld\n",total);
printf(
"***************************   Allocations:   ***********************************\n");
for (long int i=0;i<allocation_no;i++)
   printf("Allocation: %ld, amount allocated: %ld\n",i+1,alloc2_ord_size[i]);
printf(
"********************************************************************************\n");
fflush(stdout);
}
inline void allocator::initializer(long int max, double* addr) {
static long int a=0;
a++;
if (a>1 && max<real_amount) {
        printf("Error! Decreasing memory below present allocation!\n");
        abort();
        }
max_amount=max;
if (a>1 && refaddress!=addr) {
        printf("Error! Redefining starting address!\n");
        abort(); 
        }
refaddress=addr;
}
//******************************************************************************
inline void allocator::mmark() {
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
assure(marks[marker]!=-1,"Dynamic memory: the marks error 3\n");
marks[marker]=allocation_no;
marker++;
assure(marks[marker]!=-1,"Dynamic memory: the marks error 4\n");
assure(marker>=max_marks,"Dynamic memory: the marks amount exceeded\n");
}
//******************************************************************************
inline void allocator::retmark() {
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
assure(marks[marker]!=-1,"Dynamic memory: the marks error 1\n");
assure(marker==0,"Dynamic memory: the marks error 2\n");
marker--;
long int start=allocation_no;
long int stop=marks[marker];
marks[marker]=-1;
long int iamount=start-stop;
if (iamount>0) retmem(iamount);
}
//******************************************************************************
inline long int allocator::retmem(long int iamount) {
long int chunk;
long int start,stop,oldamount;
long int chunk_alloc;
long int i;
long int last_free=0;
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
start=allocation_no-1;
stop=allocation_no-1-iamount;
assure(stop<-1,"Internal memory error #7\n");
for (i=start;i>stop;i--) {
  if (locked[i]) {
    unlock(alloc2address[i]);
    }
    if (HARD_DEBUG) 
      if (mprotect(alloc2address[i]+alloc2size[i]-PAGESIZE/8,PAGESIZE,PROT_READ|PROT_WRITE)<0) {perror("mprotect error"); exit(33);}
  amount-=alloc2size[i];
  oldamount=alloc2size[i];
  alloc2size[i]=0;
  alloc2_ord_size[i]=0;
  chunk=alloc2chunkno[i];
  alloc2chunkno[i]=0;
  chunk_alloc=chunkno2addr[chunk]->free(alloc2address[i],oldamount);
  last_free=alloc2address[i]-refaddress+1;
  alloc2address[i]=0;
  if (chunk_alloc==0) {
    assure((chunk_no-1)!=chunk,"Internal memory error #6\n");
    real_amount-=chunkno2addr[chunk]->show_size();
    assure(real_amount<0,"Internal memory error #8\n");
    delete chunkno2addr[chunk];
    chunkno2addr[chunk]==NULL;
    chunk_no--;
    }
  if (allocation_no==0)  break;
  allocation_no--;
  }
return last_free;
}
//******************************************************************************
inline void allocator::getmem(const long int & iamount, long int & iaddress) {
  assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
  if (amount+iamount>max_amount) {
  for (long int i=0;i<allocation_no;++i) printf("%6ld, %9ld\n",i,alloc2size[i]);
  }
  assure(amount+iamount>max_amount,"Dynamic memory memory overflow\n");
  assure(max_alloc==allocation_no,"Dynamic memory memory requests overflow\n");
  assure(iamount<=0,"Negative or zero amount of memory requested\n");
  long int page_amount=((iamount + 128 + PAGESIZE/sizeof(double)-1) & 
                       ~(PAGESIZE/sizeof(double)-1));
  if (HARD_DEBUG) page_amount+=PAGESIZE/sizeof(double);
  long int reqest_amount;
  if (DEBUG) reqest_amount=page_amount;
  else       reqest_amount=iamount;
for (long int i=0;i<chunk_no;i++) {
  if (chunkno2addr[i]->freespace()>=reqest_amount) {
    alloc2chunkno[allocation_no]=i;
    alloc2size[allocation_no]      = reqest_amount;
    alloc2_ord_size[allocation_no] = iamount;
    alloc2address[allocation_no]=chunkno2addr[i]->allocate(reqest_amount);
    assure((reinterpret_cast<long int>(alloc2address[allocation_no])-
            reinterpret_cast<long int>(refaddress))%8!=0, 
            "Severe memory allocation error!\n");
    iaddress=alloc2address[allocation_no]-refaddress+1;
    if (HARD_DEBUG) 
      if (mprotect(alloc2address[allocation_no]+reqest_amount-PAGESIZE/8,PAGESIZE,PROT_NONE)<0) {perror("mprotect error"); exit(33);}
  fflush(stdout);
    allocation_no++;
    amount+=reqest_amount;
    return;
    }
  }
assure(chunk_no>=max_chunk,"Too many memory chunks\n");
//if (reqest_amount<=6291456) chunkno2addr[chunk_no]=new chunk(6291456);
long int isize;
if (reqest_amount<=chunk_size) isize=chunk_size;
else                     isize=reqest_amount;
if (real_amount+isize>max_amount) 
printf("real_amount %ld, amount %ld, requested %ld\n",max_amount-real_amount,max_amount-amount,reqest_amount);
assure(real_amount+isize>max_amount,
		"Free space is not available, because is dispersed...\n");
chunkno2addr[chunk_no]=::new chunk(isize);
real_amount+=isize;
alloc2chunkno[allocation_no]=chunk_no;
alloc2size[allocation_no]      = reqest_amount;
alloc2_ord_size[allocation_no] = iamount;
alloc2address[allocation_no]=chunkno2addr[chunk_no]->allocate(reqest_amount);
iaddress=alloc2address[allocation_no]-refaddress+1;
if (HARD_DEBUG) 
  if (mprotect(alloc2address[allocation_no]+reqest_amount-PAGESIZE/8,PAGESIZE,PROT_NONE)<0) {perror("mprotect error"); exit(33);}
  fflush(stdout);
allocation_no++;
chunk_no++;
amount+=reqest_amount;
return;
}
inline long int allocator::allocated(double *address) const {
for (long int i=0;i<allocation_no;i++) {
  if (alloc2address[i]==address) return alloc2_ord_size[i];
  }
return -1;
}
inline long int allocator::lock(double *address) {
long int res;
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
if (!DEBUG) return 0;
for (long int i=0;i<allocation_no;i++) {
  if (alloc2address[i]==address) {
    if (res=mprotect(address,alloc2_ord_size[i]*8,PROT_READ)) {
      perror("Couldn't mprotect");
      printf("'mprotect' error, possible system-dependent problem.\n");
      abort();
      }
    locked[i]=1;
    return alloc2_ord_size[i];
    }
  }
for (long int i=0;i<allocation_no;i++) {
  printf("%p  %p\n",alloc2address[i],address);
  }
fflush(stdout);
abort();
}
inline long int allocator::unlock(double *address) {
long int res;
assure(refaddress==NULL,"Dynamic memory allocator is not initialized\n");
if (!DEBUG) return 0;
for (long int i=0;i<allocation_no;i++) {
  if (alloc2address[i]==address) {
    if (res=mprotect(address,alloc2_ord_size[i]*8,PROT_READ|PROT_WRITE)) {
      perror("Couldn't mprotect");
      printf("Memory error no 10\n");
      abort();
      }
    locked[i]=0;
    return alloc2_ord_size[i];
    }
  }
abort();
}
//******************************************************************************
inline allocator::allocator() : allocation_no(0), chunk_no(0), refaddress(0),
                         marker(0), max_amount(0), real_amount(0)
{
for (long int i=0;i<max_alloc;i++) alloc2chunkno[i]=0;
for (long int i=0;i<max_alloc;i++) alloc2size[i]=0;
for (long int i=0;i<max_alloc;i++) alloc2_ord_size[i]=0;
for (long int i=0;i<max_alloc;i++) locked[i]=0;
for (long int i=0;i<max_alloc;i++) alloc2address[i]=0;
for (long int i=0;i<max_chunk;i++) chunkno2addr[i]=0;
for (long int i=0;i<max_marks;i++) marks[i]=-1;
}
// static allocator memo;
static allocator * memo= new allocator;
//
inline void unlocker() {
    if (HARD_DEBUG) {
      long int res;
      if (res=mprotect(memo,sizeof(allocator),PROT_READ|PROT_WRITE)) {
        perror("Couldn't mprotect");
        printf("Memory error no 110\n");
        abort();
        }
      }
}
inline void locker() {
    if (HARD_DEBUG) {
      long int res;
      if (res=mprotect(memo,sizeof(allocator),PROT_READ)) {
        perror("Couldn't mprotect");
        printf("Memory error no 110\n");
        abort();
        }
      }
}
//char*xx=operator memory::new char;
extern "C" {
#ifdef UNDERSCORE2
  void dynamic_init__(long int *max,double* pointer) {
#else
  void dynamic_init_(long int *max,double* pointer) {
#endif
      unlocker();
    memo->initializer(*max,pointer);
      locker();
  }
#ifdef UNDERSCORE2
  void dynamic_getmem__(long int *amount,long int *bladdress) {
#else
  void dynamic_getmem_(long int *amount,long int *bladdress) {
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
    memo->getmem(*amount,*bladdress);
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_retmem__(long int *pointer){
#else
  void dynamic_retmem_(long int *pointer){
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
    memo->retmem(*pointer);
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_retmem_matrix__(long int *pointer,long int *returned){
#else
  void dynamic_retmem_matrix_(long int *pointer,long int *returned){
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
    *returned=memo->retmem(*pointer);
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_dump__() {
#else
  void dynamic_dump_() {
#endif
    unlocker();
  memo->dumpall();
    locker();
  }
#ifdef UNDERSCORE2
  void dynamic_show_free__(long int *freeamount) {
#else
  void dynamic_show_free_(long int *freeamount) {
#endif
  *freeamount=memo->showfree();
  }
#ifdef UNDERSCORE2
  void dynamic_mmark__() {
#else
  void dynamic_mmark_() {
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
  memo->mmark();
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_retmark__() {
#else
  void dynamic_retmark_() {
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
  memo->retmark();
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_show_alloc__(double * address, long int * iamount) {
#else
  void dynamic_show_alloc_(double * address, long int * iamount) {
#endif
    unlocker();
  *iamount=memo->allocated(address);
    locker();
  }
#ifdef UNDERSCORE2
  void dynamic_lock__(double * address, long int * iamount) {
#else
  void dynamic_lock_(double * address, long int * iamount) {
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
  *iamount=memo->lock(address);
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void dynamic_unlock__(double * address, long int * iamount) {
#else
  void dynamic_unlock_(double * address, long int * iamount) {
#endif
    void (*previous_handler)(int);
    previous_handler=signal(SIGALRM,SIG_IGN);
    unlocker();
  *iamount=memo->unlock(address);
    locker();
    signal(SIGALRM,previous_handler);
  }
#ifdef UNDERSCORE2
  void show_pointer__(char * message, double * address, int len) {
#else
  void show_pointer_(char * message, double * address, int len) {
#endif
  char * loc_message = new char[len+1];
  for (int i=0;i<len;++i) loc_message[i]=message[i];
  loc_message[len]=0;
  printf("%s: %p, my_pid: %d\n",loc_message,(void*)address,getpid());
  fflush(stdout);
  delete[] loc_message;
  }
// for given matrix address the bl index is returned
  void pointer2bl_(double * address, long * index) { 
  double * ref;
  memo->get_ref(&ref);
  *index=address-ref+1;
  }
// Fortran call: call partition_size(path,isize)
// where "path" is a character variable and isize is a default integer.
// It returns (in the "isize" variable) the size of the partition (in bytes) 
// which contains the file or directory "path".
#include <sys/vfs.h>
#ifdef UNDERSCORE2
void partition_size__(char * path, long int * size, int len) {
#else
void partition_size_(char * path, long int * size, int len) {
#endif
struct statfs result;
char path1[len+1];
int i;
for (i=0;i<len;++i) path1[i]=path[i];
path1[len]='\000';
for (i=len;i>=0;--i) if (path1[i]!='/') 
                         path1[i]='\000' ; else break;
if (statfs(path1,&result)==0) 
   *size=result.f_bavail*result.f_bsize;
else
   *size=0;
}
void czeroit_(void*space,long * size) {
  memset(space,0,*size*sizeof(double));
  }
void cizeroit_(void*space,long * size) {
  memset(space,0,*size*sizeof(long));
  }
#ifdef UNDERSCORE2
void general_protect__(char * start, char * stop,
                       char * protected_address, char* action, long * flag) {
#else
void general_protect_(char * start, char * stop,
                       char * protected_address, char* action, long * flag) {
#endif
    using std::cout;
    using std::endl;
    typedef unsigned long ul;
    int res;
    char * p_be = (char *)(((ul) start + PAGESIZE-1) & ~(PAGESIZE-1));
    char * p_en = (char *)(((ul) stop )              & ~(PAGESIZE-1));
    if (p_en==p_be || (stop-start)<=0) { // the area for prot == 0
      *flag=0;
      cout << "Empty area supplied " 
           << (void*) start
	   << "       "
           << (void*) stop
	   << endl;
      return;
      }
    if ((p_be-protected_address)>0 || (protected_address-p_en)>=0) {
      *flag=0;
      cout << "Prot address not inside adopted area" << endl;
      return;
      }
  if (*action=='r') res=mprotect(p_be,p_en-p_be,PROT_READ);
  else if (*action=='w') res=mprotect(p_be,p_en-p_be,PROT_READ|PROT_WRITE);
  else abort();
  if (res!=0) *flag=0; else *flag=1;
  if (res!=0) { 
    perror("Couldn't mprotect");
    }
  }
void sleeper(int sig) {
fprintf(stderr,"Got signal%d\n",sig);
sleep(1000);
abort();
}
#ifdef UNDERSCORE2
void signal_default__(const long * signal_no) {
#else
void signal_default_(const long * signal_no) {
#endif
  signal(*signal_no,SIG_DFL);
//  signal(*signal_no,sleeper);
  }
#include <stdarg.h>
char * g_buffer=NULL;
int g_size=0;
//
// packer(number of parameters, bl index of resulting dest, bl(1) )
//
void packer_(long int * no_param, long int * dest, char * ref, ...) {
 using std::cerr;
 using std::endl;
 va_list ap;
 if (*no_param>5) {cerr << "Too many 'va' parameters" << endl; exit(2);}
 char* data[*no_param];
 long int sizes[*no_param];
 va_start(ap, ref);
 for (int i=0;i<*no_param;++i) {
   data[i]=va_arg(ap, char*);
   sizes[i]=*(va_arg(ap,long int*));
   }
 va_end(ap);
//
 int tot_size=0;
 for (int i=0;i<*no_param;++i) tot_size+=sizes[i];
 if (tot_size>g_size || g_buffer==NULL) {
   if (g_buffer!=NULL) delete[] g_buffer;
   g_buffer=new char[tot_size];
   g_size=tot_size;
   }
//
 for (int i=0,l=0;i<*no_param;++i) {
   for (int k=0;k<sizes[i];++k,++l) {
     g_buffer[l]=data[i][k];
     }
   if (l>g_size) {cerr << "Error in packer" << endl; exit(2);}
   }
 long int difference=g_buffer-ref;
 if (difference%8!=0) {cerr << "Error div in packer" << endl; exit(2);}
 *dest=difference/8+1;
}
//
//
void bufpointer_(long int * dest, char * ref, long int * size) {
 using std::cerr;
 using std::endl;
  if (*size>g_size || g_buffer==NULL) {
    if (g_buffer!=NULL) delete[] g_buffer;
    g_buffer=new char[*size];
    g_size=*size;
    }
  long int difference=g_buffer-ref;
  if (difference%8!=0) {cerr << "Error div in bufpointer" << endl; exit(2);}
  *dest=difference/8+1;
  }
//
//
void unpacker_(long int * no_param, ...) {
 using std::cerr;
 using std::endl;
 va_list ap;
 if (*no_param>5) {cerr << "Too many 'va' parameters" << endl; exit(2);}
 char* data[*no_param];
 long int sizes[*no_param];
 va_start(ap, no_param);
 for (int i=0;i<*no_param;++i) {
   data[i]=va_arg(ap, char*);
   sizes[i]=*(va_arg(ap,long int*));
   }
 va_end(ap);
//
 int tot_size=0;
 for (int i=0;i<*no_param;++i) tot_size+=sizes[i];
 if (tot_size>g_size || g_buffer==NULL) {
   cerr << "The buffer is too small for reading in unpacker" << endl;
   }
//
 for (int i=0,l=0;i<*no_param;++i) {
   for (int k=0;k<sizes[i];++k,++l) {
     data[i][k]=g_buffer[l];
     }
   if (l>g_size) {cerr << "Error in unpacker" << endl; exit(2);}
   }
 if (g_buffer!=NULL) {
   delete[] g_buffer;
   g_buffer=NULL;
   }
}
/*
void getcpumask_(char * strin, long *len) {
using namespace std;
cpu_set_t cpu_no;
CPU_ZERO(&cpu_no);
pid_t pid=getpid();
if (sched_getaffinity(pid,sizeof(cpu_no),&cpu_no)) {
  perror("sched_setaffinity failed: ");
  }
ostringstream xxxx;
for (int i=*len-1;i>=0;--i) {
  if (CPU_ISSET(i,&cpu_no)) xxxx << 1;
  else                      xxxx << 0;
  }
strcpy(strin,xxxx.str().c_str());
}
#ifdef UNDERSCORE2
void set_affinity__(int *aff) {
#else
void set_affinity_(int *aff) {
#endif
  using namespace std;
  cpu_set_t cpu_no;
  CPU_ZERO(&cpu_no);
  pid_t pid=getpid();
  long no_cpus=2;
  if ((no_cpus=sysconf(_SC_NPROCESSORS_CONF))<0) perror("sysconf");
  CPU_SET((*aff%no_cpus),&cpu_no);
  if (sched_setaffinity(pid,sizeof(cpu_no),&cpu_no)) {
    perror("sched_setaffinity failed: ");
    }
  }
#ifdef UNDERSCORE2
void set_allaffinity__() {
#else
void set_allaffinity_() {
#endif
  using namespace std;
  cpu_set_t cpu_no;
  CPU_ZERO(&cpu_no);
  pid_t pid=getpid();
  long no_cpus=4;
  if ((no_cpus=sysconf(_SC_NPROCESSORS_CONF))<0) perror("sysconf");
  for (int i=0;i<no_cpus;++i) CPU_SET(i,&cpu_no);
  if (sched_setaffinity(pid,sizeof(cpu_no),&cpu_no)) {
    perror("sched_setaffinity failed: ");
    }
  }
*/
void read_mmaped_memory(char * memory, long size) {
  volatile int i;
  for (long k=0;k<size;++k) i=memory[k];
  }
#ifdef UNDERSCORE2
void pid_string__(char * string, int len) {
#else
void pid_string_ (char * string, int len) {
#endif
  pid_t pid=getpid();
  snprintf(string,len,"%d",pid);
  }
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#ifdef UNDERSCORE2
void suicide__(const long * signal_no) {
#else
void suicide_(const long * signal_no) {
#endif
kill(getpid(),*signal_no);
}
}
