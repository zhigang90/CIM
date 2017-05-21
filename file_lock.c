#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>


int fd_waitandlock(int fd) {
struct flock fl;

   fl.l_type = F_WRLCK;
   fl.l_whence = SEEK_SET;
   fl.l_start = 0;
   fl.l_len = 0;

   return fcntl(fd, F_SETLKW, &fl);
}



int fd_unlock(int fd) {
struct flock fl;

   fl.l_type = F_UNLCK;
   fl.l_whence = SEEK_SET;
   fl.l_start = 0;
   fl.l_len = 0;

   return fcntl(fd, F_SETLKW, &fl);
}



int open_by_name(char *fname) {
int fd; 
    fd=open(fname, O_WRONLY | O_CREAT);
    if (fd<0) {fprintf(stderr, "Cannot open file %s\n", fname); return -1;}
    return fd;
}



int waitandlock_(char *fname, int *len) {
   fname[*len]='\0';
   return fd_waitandlock(open_by_name(fname));
}



int unlock_(char *fname, int *len) {
   fname[*len]='\0';
   return fd_unlock(open_by_name(fname));
}
