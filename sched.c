/*
 * Here are collected some utility fuctions for interacting with the operating
 * system. These routine works mostly on Linux, some works also on MAC, very
 * few are more than empty stubs on Windows
 *
 * there are functions for retrieving the PID, operating system name and
 * version, number of available CPUs. Plus there are some functions for
 * interacting with the linux scheduler in order to retrieve and set the
 * CPU affinity mask and the scheduling policy. The scheduler interaction
 * routines are used only in the parallel version
 */

#ifdef WINDOWS

#pragma once
#define WIN32_LEAN_AND_MEAN   // Exclude rarely-used stuff from Windows headers

#else

#ifndef MACINTEL

#define _GNU_SOURCE
#include <sched.h>
#ifndef SCHED_BATCH    // SCHED_BATCH should be defined in the <sched.h> header
#define SCHED_BATCH 3  // file, but if it is not, we set it to 3, that should be
#endif                 // the value for the batch scheduling policy on Linux

#endif

#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>



#ifdef UNDERSCORE2
int get_pid__(void);
#else
int get_pid_(void);
#endif
#ifdef WININTEL
int GET_NCPUS(void);
void GET_UNAME( char* );
#else
#ifdef UNDERSCORE2
int get_ncpus__(void);
void get_uname__( char* );
#else
int get_ncpus_(void);
void get_uname_( char* );
#endif
#endif
int yield_(void);
int getaffinity_( int*, int* );
int setaffinity_( int*, int* );
int setbatch_( void );
int setother_( void );

#ifdef WININTEL
void GET_UNAME( char *name )
#else
#ifdef UNDERSCORE2
void get_uname__( char *name )
#else
void get_uname_( char *name )
#endif
#endif
{
/*
 *   query for kernel name and release
 */

#ifndef WINDOWS
  int iret=0;

  struct utsname sname;
  char *blank = " \0";

  iret = uname( &sname );

  if ( iret == 0 ) {
    name[0]='\0';
    strncat( name, sname.sysname, (size_t) 255 - strlen( name ) );
    strncat( name, blank, (size_t) 255 - strlen( name ) );
  //strncat( name, sname.nodename, (size_t) 255 - strlen( name ) );
  //strncat( name, blank, (size_t) 255 - strlen( name ) );
    strncat( name, sname.release, (size_t) 255 - strlen( name ) );
    strncat( name, blank, (size_t) 255 - strlen( name ) );
  //strncat( name, sname.version, (size_t) 255 - strlen( name ) );
  //strncat( name, blank, (size_t) 255 - strlen( name ) );
    strncat( name, sname.machine, (size_t) 255 - strlen( name ) );
    name[strlen(name)]=' ';
  }
#endif

}

#ifdef UNDERSCORE2
int get_pid__(void)
#else
int get_pid_(void)
#endif
{
  int pid=0;
#ifndef WINDOWS
  pid = (int) getpid();
#endif
  return pid;
}

int setbatch_( void )
{
#ifndef MACINTEL

#ifndef WINDOWS
  struct sched_param sparam;
  sparam.sched_priority = 0;
  return sched_setscheduler( 0, SCHED_BATCH, &sparam );
#else
  // empty stub on Windows
  return 0;
#endif

#else

  // empty stub on MAC
  return 0;

#endif
}

int setother_( void )
{
#ifndef MACINTEL
#ifndef WINDOWS
  struct sched_param sparam;
  sparam.sched_priority = 0;
  return sched_setscheduler( 0, SCHED_OTHER, &sparam );
#else
  // empty stub on Windows
  return 0;
#endif
#else
  // empty stub on MAC
  return 0;
#endif
}

int getaffinity_( int *ncpu, int *icpu )
{
#ifndef MACINTEL
#ifndef WINDOWS
  cpu_set_t mask;
  size_t cpusetsize = (size_t) sizeof( mask );
  int iret, i;

  CPU_ZERO( &mask );

  iret = sched_getaffinity( 0, cpusetsize, &mask );

  if( iret == 0 ) {
    for ( i = 0; i <= *ncpu-1; i++ ){
      icpu[i] = CPU_ISSET( i, &mask );
    }
  }

  return iret;
#else
  // empty stub on Windows
  return 0;
#endif
#else
  // empty stub on MAC
  return 0;
#endif
}

int setaffinity_( int *ncpu, int *icpu )
{
#ifndef MACINTEL
#ifndef WINDOWS
  cpu_set_t mask;
  size_t cpusetsize = (size_t) sizeof( mask );
  int iret, i;

  CPU_ZERO( &mask );
  for ( i = 0; i <= *ncpu-1; i++ ){
    if( icpu[i] != 0) {
      CPU_SET( i, &mask );
    }
  }

  iret = sched_setaffinity( 0, cpusetsize, &mask );

  return iret;
#else
  // empty stub on Windows
  return 0;
#endif
#else
  // empty stub on MAC
  return 0;
#endif
}

int yield_(void)
{
#ifndef MACINTEL
#ifndef WINDOWS
   return sched_yield();
#else
  // empty stub on Windows
  return 0;
#endif
#else
  // empty stub on MAC
  return 0;
#endif
}

#ifdef WININTEL
int GET_NCPUS(void)
#else
#ifdef UNDERSCORE2
int get_ncpus__(void)
#else
int get_ncpus_(void)
#endif
#endif
{
/*
 *  get the number of processors of the current host
 */

  int ncpus=0, ncpur;
#ifdef WINDOWS
  char *var = NULL;
#endif

#ifndef MACINTEL

#ifdef WINDOWS  // on windows we test the environment
  var = getenv("NUMBER_OF_PROCESSORS");
  if( !var ) {
    sscanf( var, "%d", &ncpur ); 
   }
  else {
    ncpur = 0;
  }
#else    // on Linux we call sysconf
  ncpur = sysconf( _SC_NPROCESSORS_ONLN );
#endif

#else  // on MAC we call MPProcessors (part of Carbon framework)
  ncpur = MPProcessors();
#endif

  if( ncpur >= 0 ) {
    ncpus = ncpur;
  }

  return ncpus;
}

#ifdef MAIN
int main(int argc, char **argv)
{
  cpu_set_t mask;
  size_t cpusetsize = (size_t) sizeof( mask );
  int iret, imax, err;
  struct sched_param sparam;
  struct utsname sname;

  sparam.sched_priority = 0;

  iret = sched_getaffinity( 0, cpusetsize, &mask );

  printf( "iret %d\n", iret );
  printf( "cpu 0 %d\n", CPU_ISSET( 0, &mask ) );
  printf( "cpu 1 %d\n", CPU_ISSET( 1, &mask ) );
  printf( "cpu 2 %d\n", CPU_ISSET( 2, &mask ) );
  printf( "cpu 3 %d\n", CPU_ISSET( 3, &mask ) );
  printf( "cpu 4 %d\n", CPU_ISSET( 4, &mask ) );
  printf( "cpu # %d\n", get_ncpus_() );
  printf( "max %d\n", CPU_SETSIZE );

  iret = sched_getscheduler( 0 );
  printf( "iret %d\n", iret );
  printf( " SCHED_FIFO %d\n", SCHED_FIFO );
  printf( " SCHED_RR %d\n", SCHED_RR );
  printf( " SCHED_OTHER %d\n", SCHED_OTHER );
  printf( " SCHED_BATCH %d\n", SCHED_BATCH );

  iret = sched_setscheduler( 0, SCHED_BATCH, &sparam );
  err=errno;
  printf( "iret %d\n", iret );
  printf( "errno %d\n", err );
  printf( " EINVAL %d\n", EINVAL );
  printf( " EPERM %d\n", EPERM );
  printf( " ESRCH %d\n", ESRCH );

  iret = sched_getscheduler( 0 );
  printf( "iret %d\n", iret );

  iret = uname( &sname );
  printf( "iret for uname %d\n", iret );
  if ( iret == 0 ) {
    printf("sysname %s\n", sname.sysname );
    printf("nodename %s\n", sname.nodename );
    printf("release %s\n", sname.release );
    printf("version %s\n", sname.version );
    printf("machine %s\n", sname.machine );
#ifdef _GNU_SOURCE
    printf("domainname %s\n", sname.domainname );
#endif

  }
}

#endif
