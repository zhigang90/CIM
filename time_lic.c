/*
 * This file contains the PQS license-handling routines.
 *
 * This is a new version written in C, mainly for compatibility with
 * PQSMOL, so that we do not need a Fortran compiler for compiling it
 * (we would need a C compiler for PQS anyway, for the mac address bit).
 *
 */


#ifdef WINDOWS

#pragma once
#define WIN32_LEAN_AND_MEAN   // Exclude rarely-used stuff from Windows headers
#include <Windows.h>
#include <Iphlpapi.h>
#include <Assert.h>
#include <tchar.h>

#else

  #ifdef MACINTEL

  #include <CoreFoundation/CoreFoundation.h>
  
  #include <IOKit/IOKitLib.h>
  #include <IOKit/network/IOEthernetInterface.h>
  #include <IOKit/network/IONetworkInterface.h>
  #include <IOKit/network/IOEthernetController.h>

  
  #else
  
  #include <sys/ioctl.h>
  #include <net/if.h>
  #include <unistd.h>
  
  #endif

#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifndef VERSION
#define VERSION "Development"
#endif

#ifndef REVISION
#define REVISION "Beta"
#endif

#ifndef IDENT
#define IDENT "None"
#endif

#ifndef ARCHIT
#define ARCHIT "Generic"
#endif

char *MD5="some_crap\0";

#ifndef WINDOWS
int pqs_setenv_( const char *name_inp, int *l_name, const char *value_inp, int *l_value );
#endif

#ifdef WININTEL
int  CHKLICENSE( const char *product );
int LOCKCODE_FILE_( const char *product );
#else
int  chklicense_( const char *product );
int lockcode_file_( const char *product );
#endif
int lockcode_( const char *product );
int intxor( int xorini, void *field, int nbytes );

#ifdef WINDOWS
long macaddr_(unsigned char *MAC); // Windows version
#else
  #ifdef MACINTEL
  long macaddr_( UInt8* MAC);  // Mac Intel version
  static kern_return_t FindEthernetInterfaces(io_iterator_t *matchingServices);
  static kern_return_t GetMACAddress(io_iterator_t intfIterator, UInt8 *MACAddress, UInt8 bufferSize);
  #else
  long macaddr_(char *MAC);   // Linux version
  #endif
#endif
void lockcode_line(const char *product, char *lock_line);
void lockcode_file_text(const char *product, char *lock_line);

unsigned int checksum( void *buffer, int length, unsigned int sum32);
void char_encode( unsigned int value, char *ascii);

#ifndef WINDOWS
int pqs_setenv_( const char *name_inp, int *l_name, const char *value_inp, int *l_value )
{
/*
 * interface routine for manipulating environment variables
 * this might be necessary in the parallel version if the
 * environment is not properly set, and if the values are
 * different from the default
 */
  int iret = -1;
  char name[257], value[257];

       // sanity check

  if ( *l_name > 256 ) return iret;
  if ( *l_value > 256 ) return iret;

  strncpy( name, name_inp, (size_t) *l_name);
  name[*l_name]='\0';
  strncpy( value, value_inp, (size_t) *l_value);
  value[*l_value]='\0';

  iret = setenv( name, value, 0 ); // set the overwrite parameter to 0,
                                   // i.e., do not ovevrwrite, but set
                                   // the variable only if the if it is not
                                   // already set
  return iret;
}
#endif

#ifdef WININTEL
int  CHKLICENSE( const char *product )
#else
int  chklicense_( const char *product )
#endif
{
/*
 * PQS license checker. If a valid license is found, the routine
 * returns the maximum number of processes allowed. If a license is
 * not found, it returns 0
 *
 * product    product string to validate license for
 * 
 * A license line looks like this:
 *
 * 259271766 243082838 243082582 PQS parallel 0308 0999 0008 ca4VfU1UcZ1UcZ1U
 *    (1)       (2)        (3)   (4)  (5)      (6)  (7)  (8)        (9)
 *
 *   1)  encoding of expiration date
 *   2)  encoding of issue date
 *   3)  encoding of maximum number of processors (optional field)
 *   4)  comment string
 *   5)  ditto
 *   6)  issue date
 *   7)  expiration date
 *   8)  maximum number of processors
 *   9)  checksum of fields 1--8, including spaces
 *
 *   If all the above fields are presents the function checks if the checksum (9)
 *   matches the remaining of the line, to make sure the license has not been 
 *   manipulated, then it tests fields (1), (2) and, optionally, (3).
 *
 *   Field 3 (maximum number of processors) is optional. If it is not found, the
 *   corresponding value is set to 9999
 *
 *   If only fields (1) and (2) are presents (old style license), then the 
 *   checksum test is not performed.
 *
 *   All in all, this license system is pretty weak. In fact it has already
 *   happened that in large cluster, the same license line is validated by more
 *   than one node. This might be due to the fact that MAC addresses in clusters
 *   are usually very similar, or it is just a plain weakeness of the algorithm.
 */

#ifdef WINDOWS
  char pqsroot[256] = "c:\\Program Files\\PQS\\PQS 3.3";  // default PQS_ROOT
  char DIR_SEP[2] = "\\";
#else
  char pqsroot[256] = "/usr/local/share/PQS";  // default PQS_ROOT
  char DIR_SEP[2] = "/";
#endif

  char *env = NULL, filename[256], line[256], linet[256], *field;
  char *argvi[256];  // array for parsing license line
  int argc, argi, lent, lenl, i;
  unsigned int sum32;
  char encode[17], encoder[17];

  FILE* licfile;
  int licf = 0, licp;
  int iexpd, issue, lock, month, year, monthl, yearl, datel, mpl, maxproc;
  char dstrl[5], mpstrl[5];
  struct tm *tptr;
  time_t sec;

                    // overwrite the default if PQS_ROOT environment is defined

  env = getenv( "PQS_ROOT" );
  if( env != NULL ) { 
    strcpy( pqsroot, env );
  }
  else {
    fprintf( stderr, "\n ***Warning: PQS_ROOT not defined\n\n" );
  }

                  // append a DIR_SEP to PQS_ROOT, if needed

  if( pqsroot[strlen(pqsroot)-1] != DIR_SEP[0] ) { strcat( pqsroot, DIR_SEP ); }

                   // open license file

  strcpy( filename, pqsroot );
  strcat( filename, "pqs_lic" );
  if ( ( licfile = fopen(  filename, "r" ) ) != NULL ) {

                   // get lockcode

    lock = lockcode_( product );

                   // get current date

    time( &sec );
    tptr = localtime( &sec );
    month = (int) tptr->tm_mon + 1;
    year  = (int) tptr->tm_year + 1900;
    year = year - 2000;  //  PQS counts years starting from 2000

                  // read license file

    while( !feof( licfile ) ) {
      fgets( line, 255, licfile );

      if( !feof(licfile) ) {  // new license line

        licp = 0;

        lenl = strlen(line);
        if( lenl > 256 ) { 
           lenl=256;
           line[255] = '\0';
        }
        if( lenl < 10 ) { continue; } // line is too short, do not bother

        i = strspn( line, " " );  // trim spaces at the beginning of the line
        if( i > 0 ) {
          field = &line[i];
          strcpy( linet, field);
          strcpy( line, linet );
        }

        if( !isdigit( line[0] ) && line[0] != '-' ) { continue; } // line must start
                                                                  // with a number
        i = strlen( line ) - 1;        // trim spaces at the end of the line
        while( isspace( line[i] ) && i > 0 ) {
          line[i] = '\0';
          i--;
        }

        if( strlen( line ) < 10 ) { continue; } // line is too short, do not bother

        strncpy( linet, line, 255 );    // splits line into fields
        line[255] = '\0';
        argc = 0;
        argvi[0] = NULL;
        field = strtok( linet, " " );
        while( field != NULL ) {
          argc++;
          argvi[argc-1] = field;
          field = strtok( NULL, " " );
        }
        if( argc < 2 ) { continue; }  // we need at least two fields 

        if( argc > 2 ) {   // checksum test
          argi = argc - 1;
          while( strlen( argvi[argi] ) < 16 && argi > 2 ) { argi--; }
          strncpy( encoder, argvi[argi], 16 );
          encoder[16] = '\0';
          lent = argvi[argi] - argvi[0] - 1;
          lent = lent /4;
          sum32 = checksum( line, lent*4, 0 );
          char_encode( ~sum32, encode );
          encode[16] = '\0';
          if( strcmp( encode, encoder) ) { continue; }
        }

        sscanf( line, "%d%d", &iexpd, &issue );

                // recover issue date

        monthl = 0;
        yearl = 999;
        datel = intxor( lock, (void *) &issue, sizeof(int) );
        strncpy( dstrl, (char *) &datel, 4 );
        dstrl[4] = '\0';

               // a valid date is composed only of digits and spaces

        if( !isdigit( dstrl[0] ) && dstrl[0] != ' ' ) { continue; }
        if( !isdigit( dstrl[1] ) && dstrl[1] != ' ' ) { continue; }
        if( !isdigit( dstrl[2] ) && dstrl[2] != ' ' ) { continue; }
        if( !isdigit( dstrl[3] ) && dstrl[3] != ' ' ) { continue; }

        sscanf( dstrl, "%2d%2d", &monthl, &yearl );

        if( monthl <= 0 || monthl > 12 ) { continue; }  // invalid month

        if( year < yearl ||                                   // curent date is
          ( year == yearl && month < monthl ) ) { continue; } // earlier than issue

                // recover expiration date

        monthl = 0;
        yearl = -999;
        datel = intxor( lock, (void *) &iexpd, sizeof(int) );
        strncpy( dstrl, (char *) &datel, 4 );
        dstrl[4] = '\0';

               // a valid date is composed only of digits and spaces

        if( !isdigit( dstrl[0] ) && dstrl[0] != ' ' ) { continue; }
        if( !isdigit( dstrl[1] ) && dstrl[1] != ' ' ) { continue; }
        if( !isdigit( dstrl[2] ) && dstrl[2] != ' ' ) { continue; }
        if( !isdigit( dstrl[3] ) && dstrl[3] != ' ' ) { continue; }

        sscanf( dstrl, "%2d%2d", &monthl, &yearl );

        if( monthl <= 0 || monthl > 12 ) { continue; }  // invalid month

        if( year > yearl ||                                 
          ( year == yearl && month > monthl ) ) { continue; } // expired license

              // at this point the license is granted. Now compute the
              // maximum number of processors. If field 3 is present, 
              // the maximum number of processors is decoded, otherwise
              // it is set to the default value of 9999 (i.e., unlimited)

        licp = 9999;

        if( argc > 2 ) {
          field=argvi[2];        // extract third field
          if( isdigit( field[0] ) || field[0] == '-' ) { // field must start with a number
            sscanf( field, "%d", &maxproc );

                // recover maximum number of processors

            mpl = intxor( lock, (void *) &maxproc, sizeof(int) );
            strncpy( mpstrl, (char *) &mpl, 4 );
            mpstrl[4]='\0';
            sscanf( mpstrl, "%d", &maxproc );
            licp=maxproc;
          }

        }
      }
      if( licp > licf ) { licf = licp; }
      
                // if the maximum number of processor is not unlimited, keep reading
                // license lines, just in case a false hit has been generated (might happen
                // in large clusters, with lots of similar MAC addresses)

      if( licf >= 9999 ) { break; }
    }
    fclose( licfile );

  }

  return licf;
}

int lockcode_( const char *product )
{
/*
 * PQS lockcode generator
 *
 * product    product string to generate the lockcode
 *
 * This function returns a lockcode as an integer number
 * obtained by doing a bitwise exlusive or of the 16 characters
 * of the product string (4 bytes at a time) and the
 * 16 charactres of the mac address (4 bytes at a time)
 */

  char hwaddr[17], mac[7];
  int i, ip, xor;
  long macret;

           // get the mac address

  macret = macaddr_( mac );
  if( macret == (long) -1 ) { return 0; }

           // fill the hwaddr string in the proper format (i.e: "0:E0:81:34:5E:CE")

  sprintf( &hwaddr[0], "%1X", mac[0] & 0xff );
  ip = 1;
  for( i=1; i <= 5; i++ ) {
    sprintf( &hwaddr[ip], ":%02X", mac[i] & 0xff );
    ip+=3;
  }
  hwaddr[16] = '\0';

     // xor the  mac address

  xor = intxor( 0, (void *) hwaddr, 16*sizeof(char) );

     // add the product code

  xor = intxor( xor , (void *) product, 16*sizeof(char) );

  return xor;
}

int intxor( int xorini, void *field, int nbytes )
{
/* 
 * This function takes a memory area of length nbytes
 * and starting at the location pointed to by field, treats the
 * area as an integer array, and returns the bitwise exclusive or
 * of all the array locations with the initial value xorini.
 *
 * nbytes should be a multiple of sizeof(int)
 */
  int *intfld = (int *) field;
  int n = nbytes / sizeof( int ), i;
  int retval = xorini;

  for( i=1; i <= n; i++) {
    retval ^= *intfld;
    ++intfld;
  }

  return retval;
}

#ifdef WINDOWS

                                 // windows code for MAC ADDRESS

long macaddr_(unsigned char *MAC) {
IP_ADAPTER_INFO AdapterInfo[16], *pAdapterInfo;
DWORD dwBufLen=sizeof(AdapterInfo);
DWORD dwStatus;

   dwStatus=GetAdaptersInfo(AdapterInfo, &dwBufLen);	
   if (dwStatus!=ERROR_SUCCESS) return -1;    // Verify return value is valid

   pAdapterInfo=AdapterInfo; 
   memcpy(MAC, pAdapterInfo->Address, 6);

   return 0;
}
#else

#ifdef MACINTEL

                                 // MAC Intel code for MAC ADDRESS

// Returns an iterator containing the primary (built-in) Ethernet interface. The caller is responsible for
// releasing the iterator after the caller is done with it.
static kern_return_t FindEthernetInterfaces(io_iterator_t *matchingServices)
{
    kern_return_t    kernResult; 
    CFMutableDictionaryRef  matchingDict;
    CFMutableDictionaryRef  propertyMatchDict;
    
    // Ethernet interfaces are instances of class kIOEthernetInterfaceClass. 
    // IOServiceMatching is a convenience function to create a dictionary with the key kIOProviderClassKey and 
    // the specified value.
    matchingDict = IOServiceMatching(kIOEthernetInterfaceClass);

    // Note that another option here would be:
    // matchingDict = IOBSDMatching("en0");
        
    if (NULL == matchingDict) {
        printf("IOServiceMatching returned a NULL dictionary.\n");
    }
    else {
        // Each IONetworkInterface object has a Boolean property with the key kIOPrimaryInterface. Only the
        // primary (built-in) interface has this property set to TRUE.
        
        // IOServiceGetMatchingServices uses the default matching criteria defined by IOService. This considers
        // only the following properties plus any family-specific matching in this order of precedence 
        // (see IOService::passiveMatch):
        //
        // kIOProviderClassKey (IOServiceMatching)
        // kIONameMatchKey (IOServiceNameMatching)
        // kIOPropertyMatchKey
        // kIOPathMatchKey
        // kIOMatchedServiceCountKey
        // family-specific matching
        // kIOBSDNameKey (IOBSDNameMatching)
        // kIOLocationMatchKey
        
        // The IONetworkingFamily does not define any family-specific matching. This means that in            
        // order to have IOServiceGetMatchingServices consider the kIOPrimaryInterface property, we must
        // add that property to a separate dictionary and then add that to our matching dictionary
        // specifying kIOPropertyMatchKey.
            
        propertyMatchDict = CFDictionaryCreateMutable(kCFAllocatorDefault, 0,
                            &kCFTypeDictionaryKeyCallBacks,
                            &kCFTypeDictionaryValueCallBacks);
    
        if (NULL == propertyMatchDict) {
            printf("CFDictionaryCreateMutable returned a NULL dictionary.\n");
        }
        else {
            // Set the value in the dictionary of the property with the given key, or add the key 
            // to the dictionary if it doesn't exist. This call retains the value object passed in.
            CFDictionarySetValue(propertyMatchDict, CFSTR(kIOPrimaryInterface), kCFBooleanTrue); 
            
            // Now add the dictionary containing the matching value for kIOPrimaryInterface to our main
            // matching dictionary. This call will retain propertyMatchDict, so we can release our reference 
            // on propertyMatchDict after adding it to matchingDict.
            CFDictionarySetValue(matchingDict, CFSTR(kIOPropertyMatchKey), propertyMatchDict);
            CFRelease(propertyMatchDict);
        }
    }
    
    // IOServiceGetMatchingServices retains the returned iterator, so release the iterator when we're done with it.
    // IOServiceGetMatchingServices also consumes a reference on the matching dictionary so we don't need to release
    // the dictionary explicitly.
    kernResult = IOServiceGetMatchingServices(kIOMasterPortDefault, matchingDict, matchingServices);    
    if (KERN_SUCCESS != kernResult) {
        printf("IOServiceGetMatchingServices returned 0x%08x\n", kernResult);
    }
        
    return kernResult;
}
    
// Given an iterator across a set of Ethernet interfaces, return the MAC address of the last one.
// If no interfaces are found the MAC address is set to an empty string.
// In this sample the iterator should contain just the primary interface.
static kern_return_t GetMACAddress(io_iterator_t intfIterator, UInt8 *MACAddress, UInt8 bufferSize)
{
    io_object_t    intfService;
    io_object_t    controllerService;
    kern_return_t  kernResult = KERN_FAILURE;
    
    // Make sure the caller provided enough buffer space. Protect against buffer overflow problems.
  if (bufferSize < kIOEthernetAddressSize) {
    return kernResult;
  }
  
  // Initialize the returned address
    bzero(MACAddress, bufferSize);
    
    // IOIteratorNext retains the returned object, so release it when we're done with it.
    while (intfService = IOIteratorNext(intfIterator))
    {
        CFTypeRef  MACAddressAsCFData;        

        // IONetworkControllers can't be found directly by the IOServiceGetMatchingServices call, 
        // since they are hardware nubs and do not participate in driver matching. In other words,
        // registerService() is never called on them. So we've found the IONetworkInterface and will 
        // get its parent controller by asking for it specifically.
        
        // IORegistryEntryGetParentEntry retains the returned object, so release it when we're done with it.
        kernResult = IORegistryEntryGetParentEntry(intfService,
                           kIOServicePlane,
                           &controllerService);
    
        if (KERN_SUCCESS != kernResult) {
            printf("IORegistryEntryGetParentEntry returned 0x%08x\n", kernResult);
        }
        else {
            // Retrieve the MAC address property from the I/O Registry in the form of a CFData
            MACAddressAsCFData = IORegistryEntryCreateCFProperty(controllerService,
                                 CFSTR(kIOMACAddress),
                                 kCFAllocatorDefault,
                                 0);
            if (MACAddressAsCFData) {
//                CFShow(MACAddressAsCFData); // for display purposes only; output goes to stderr
                
                // Get the raw bytes of the MAC address from the CFData
                CFDataGetBytes(MACAddressAsCFData, CFRangeMake(0, kIOEthernetAddressSize), MACAddress);
                CFRelease(MACAddressAsCFData);
            }
                
            // Done with the parent Ethernet controller object so we release it.
            (void) IOObjectRelease(controllerService);
        }
        
        // Done with the Ethernet interface object so we release it.
        (void) IOObjectRelease(intfService);
    }
        
    return kernResult;
}

long macaddr_( UInt8* MACAddress )
{
    kern_return_t  kernResult = KERN_SUCCESS; // on PowerPC this is an int (4 bytes)
/*
 *  error number layout as follows (see mach/error.h and IOKit/IOReturn.h):
 *
 *  hi                lo
 *  | system(6) | subsystem(12) | code(14) |
 */

    io_iterator_t  intfIterator;
 
    kernResult = FindEthernetInterfaces(&intfIterator);
    
    if (KERN_SUCCESS != kernResult) {
        printf("FindEthernetInterfaces returned 0x%08x\n", kernResult);
    }
    else {
        kernResult = GetMACAddress(intfIterator, MACAddress, kIOEthernetAddressSize );
        
        if (KERN_SUCCESS != kernResult) {
            printf("GetMACAddress returned 0x%08x\n", kernResult);
        }
//        else {
//      printf("This system's built-in MAC address is %02x:%02x:%02x:%02x:%02x:%02x.\n",
//          MACAddress[0], MACAddress[1], MACAddress[2], MACAddress[3], MACAddress[4], MACAddress[5]);
//        }
    }
    
    (void) IOObjectRelease(intfIterator);  // Release the iterator.
        
    return (long) kernResult;
}

#else

                                           // Linux code for MAC Address

long macaddr_(char *MAC) {
    int sd, ret;
    struct ifreq req;
    
    sd=socket(AF_INET,SOCK_DGRAM,0);
    if (sd==-1) 
       {perror("socket()"); return -1;}
    
    sprintf(req.ifr_name, "eth0");
    ret=ioctl(sd,SIOCGIFHWADDR,&req);

    if (ret==-1) { 
       sprintf(req.ifr_name, "eth1");
       ret=ioctl(sd,SIOCGIFHWADDR,&req);
    }
     
    if (ret==-1) { 
       sprintf(req.ifr_name, "eth2");
       ret=ioctl(sd,SIOCGIFHWADDR,&req);
    }
    
    if (ret==-1) {
       fprintf(stderr, "Unable to get a valid machine ID!\n");
       sprintf(MAC, "JUNK");
       close(sd);
       return -1;
    }
    
    //bcopy(req.ifr_hwaddr.sa_data, MAC, 6); 
    memcpy(MAC, req.ifr_hwaddr.sa_data, (size_t) 6); 
    /*sprintf(MAC,"%02x%02x%02x%02x%02x%02x",
                        (int)(req.ifr_hwaddr.sa_data[0] & 0xff),
                        (int)(req.ifr_hwaddr.sa_data[1] & 0xff),
                        (int)(req.ifr_hwaddr.sa_data[2] & 0xff),
                        (int)(req.ifr_hwaddr.sa_data[3] & 0xff),
                        (int)(req.ifr_hwaddr.sa_data[4] & 0xff),
                        (int)(req.ifr_hwaddr.sa_data[5] & 0xff)); 
                        */
    close(sd);
    return 0; 
}
#endif
#endif

#ifdef WININTEL
int LOCKCODE_FILE( const char *product )
#else
int lockcode_file_( const char *product )
#endif
{
  FILE *lockfil;
  char *filnam = "pqs_lockcode";
  char lock_line[256] ="";
  int ierr=0;

            //  write file header, if necessary

  if( ( lockfil = fopen( filnam, "r" ) ) != NULL ){
    fclose( lockfil );
  }
  else{
    if( ( lockfil = fopen( filnam, "w" ) ) != NULL ){
      fprintf( lockfil, 
 "#\n"
 "# PQS lockcode file\n"
 "# -----------------\n"
 "#\n"
 "# This file contains the PQS lockcode(s).\n"
 "#\n"
 "# In order to obtain a valid license, please fill in\n"
 "# the contact information below, then e-mail the file\n"
 "# to licenses@pqs-chem.com.\n"
 "#\n"
 "#  Title          :\n"
 "#  First Name     :\n"
 "#  Last Name      :\n"
 "#\n"
 "#  Organization   :\n"
 "#\n"
 "#  Address 1      :\n"
 "#  Address 2      :\n"
 "#  Address 3      :\n"
 "#  City           :\n"
 "#  State/Province :\n"
 "#  Zip code       :\n"
 "#  Country        :\n"
 "#\n"
 "#  E-mail         :\n"
 "#  Phone          :\n"
 "#  Fax            :\n"
 "#\n"
 "#  License(s) requested (check the corresponding boxes):\n"
 "#\n"
 "#             PQS serial version   [ ]\n"
 "#             PQS parallel version [ ]\n"
 "#             NBO add-on module    [ ]\n"
 "#             SQM                  [ ]\n"
 "#             Convert2             [ ]\n"
 "#             PQSMol (PQS GUI)     [ ]\n"
 "#\n"
 "#             I Accept the License Agreement        [ ]\n"
 "#             I DO NOT Accept the License Agreement [ ]\n"
 "#\n"
 "#---------------------------------------------------------------\n"
 "# LICENSE AGREEMENT FOR SINGLE PROCESSOR VERSION OF PQS & PQSMol\n"
 "#\n"
 "# Parallel Quantum Solutions, LLC (hereafter PQS, LLC) is a supplier\n" 
 "# of parallel hardware/software systems (the QuantumCube(TM)) with a high\n" 
 "# performance/price ratio for quantum chemical and related calculations.\n"
 "#\n"
 "# PQS is a program for performing parallel ab initio computations in \n"
 "# quantum chemistry. We also supply a single-processor version of PQS \n"
 "# for Linux and/or Windows operating systems as well as a graphical user\n"
 "# interface, PQSMol. Both PQS and PQSMol are integral parts of the \n"
 "# Parallel QuantumCube(TM) system supplied by PQS, LLC.\n"
 "#\n"
 "# We work hard to ensure the trouble-free operation of our code. However, \n"
 "# no warranty of any kind as to the accuracy or reliability of the \n"
 "# numbers computed by PQS is given, and PQS, LLC accepts no responsibility \n"
 "# regarding any scientific or other conclusions drawn by the user from the \n"
 "# numbers produced by PQS. We will try to eliminate errors from the code \n"
 "# as we become aware of them but we cannot offer legal guarantees to solve \n"
 "# all possible problems within a certain time period.\n"
 "#\n"
 "# The user warrants that PQS, PQSMol and any other programs supplied by \n"
 "# PQS, LLC will be used only on a single CPU and on no other computer or \n"
 "# CPU, unless PQS, LLC provides written (or email) permission to the \n"
 "# contrary. Note that all programs supplied by PQS, LLC are normally \n"
 "# delivered to the end user only, and may not be redistributed and/or \n"
 "# resold without written (or email) permission to the contrary.\n"
 "#\n"
 "# Users of our Parallel QuantumCube(TM) and PQS and PQSMol programs agree (1) \n"
 "# that they will provide a general informal evaluation of the system \n"
 "# and/or the programs to interested third parties upon request, and (2) \n"
 "# they will cite our Parallel QuantumCube(TM) and/or PQS and PQSMol programs \n"
 "# in all publications and reports in which our hardware/software was used.\n"
 "#---------------------------------------------------------------\n"
 "#\n"
 "# Please wisit our web site at http://www.pqs-chem.com for a\n"
 "# complete list of our software offerings.\n"
 "#\n"
 "#-------------------- PQS lockcode file -------------------\n"
 "#\n"
 "# Please do not modify the file content below this point!\n"
 "#\n"
 "#-------------------- PQS lockcode file -------------------\n"
 "#\n");
      fclose( lockfil );
    }
    else{ ierr = 1; }
  }

  if( ( lockfil = fopen( filnam, "a" ) ) != NULL ){

    lockcode_line(product, lock_line);
    fprintf( lockfil, "%s\n", lock_line );
    fclose( lockfil );
  }
  else{ ierr = 1; }

  return ierr;
}

void lockcode_file_text(const char *product, char *ltext) {
char lock_line[256] ="";
char lock_text[]=
"#  Title          :\n"
"#  First Name     :\n"
"#  Last Name      :\n"
"#  Organization   :\n"
"#  Address 1      :\n"
"#  Address 2      :\n"
"#  Address 3      :\n"
"#  City           :\n"
"#  State/Province :\n"
"#  Zip code       :\n"
"#  Country        :\n"
"#  E-mail         :\n"
"#  Phone          :\n"
"#  Fax            :\n"
"#  License(s) requested (check the corresponding boxes):\n"
"#             PQS serial version   [X]\n"
"#             PQS parallel version [ ]\n"
"#             NBO add-on module    [ ]\n"
"#             SQM                  [ ]\n"
"#             Convert2             [ ]\n"
"#             PQSMol (PQS GUI)     [X]\n"
"#----------------------------------------------------------\n"
"# Please do not modify the file content below this point!\n"
"#----------------------------------------------------------\n";
  
  lockcode_line(product, lock_line);
  sprintf(ltext, "%s%s\n", lock_text, lock_line);
  return;
}


void lockcode_line(const char *product, char *lock_line) {
/* 
 * This function generates a lockcode line
 *
 *  A lockcode line looks like this:
 *
 *  758920707 slater:4 51876 L x86_64 pqsmol 3.3-6-T 11/9/07 YdHjZZEgYbEgYZEg
 *    (1)      (2)    (3) (4)  (5)   (6)     (7)     (8)          (9)
 *
 *   1) lockcode (machine id encoded with the lockcode key)
 *   2) hostname the lockcode was generated on ("Unknown" if not found)
 *      followed by : and number of cores [added by PDW 10/27/2010]
 *   3) kind of random number, so that it will be rather difficult to
 *      generate two identical lockcode lines.
 *   4) operating system L=Linux, W=Windows
 *   5) architecture (x86_64, i386, win32)
 *   6) name of the program generating the lockcode (pqs, pqsmol, sqm)
 *   7) version-revision
 *   8) date of the lockcode generation
 *   9) checksum of fields 1--8 including spaces
 */
char host[80] ="";
char line[256] ="";
int l4;
int numCPU;
struct tm *tptr;
time_t sec;
unsigned int id32 = 0;
int lock;
char chksum[17];
int day, month, year;

#ifdef WINDOWS
  char *pqssys = "W";  // operating system
#else
  #ifdef MACINTEL
    char *pqssys = "M";  // operating system
  #else
    char *pqssys = "L";    // operating system
  #endif
#endif
  char *pqsver = VERSION;
  char *pqsrev = REVISION;
  char *pqsid = IDENT;
  char *pqsarch = ARCHIT;
#ifdef TRIAL
  char *pqstr = "-T";  // trial version
#else
  char *pqstr = "";    // full version
#endif

#ifdef WINDOWS
  SYSTEM_INFO sysinfo;
  char *var;
#else
  int iret;
#endif

              // get host name

#ifdef WINDOWS
  var = getenv("COMPUTERNAME");
  if( !var ) {
    var = getenv("HOSTNAME");
    if( !var ) var = getenv("HOST");
   }
  if( !var ) {
    strcpy( host, "Unknown" );
  }
  else {
    strcpy( host, var);
  }
  GetSystemInfo( &sysinfo );
  numCPU = sysinfo.dwNumberOfProcessors;
#else
    iret = gethostname( host, 79 );
    if( iret < 0 ) strcpy( host, "Unknown" );
    numCPU = sysconf( _SC_NPROCESSORS_ONLN );
#endif

            // get lockcode

   lock = lockcode_( product );

            // get date

   time( &sec );
   tptr = localtime( &sec );
   day   = (int) tptr->tm_mday ;
   month = (int) tptr->tm_mon + 1;
   year  = (int) tptr->tm_year - 100;
   sprintf( line, "%d %s:%d %ld %s %s %s %s-%s%s %d/%d/%.2d", lock, host, numCPU, (long)(sec*sec)%100000, pqssys, pqsarch, pqsid, pqsver, pqsrev, pqstr, month, day, year );

                  // compute the checksum of the longest substring
                  // of line whose length is a multiple of 4
                 
   l4 = strlen(line) / 4;
   id32 = checksum( (void *)line, l4*4, id32 );
   char_encode( ~id32, chksum);
             
                 // write the line computed above and the checksum
                 // on the lockcode line

   sprintf(lock_line, "%s %s\n", line, chksum );

  return;
}


/*
 *  functions for computing checksums and for
 *  obtaining an ascii string based on ckecksum values
 */

unsigned int checksum (
    void *buffer,        /* Input array of bytes to be checksummed */
                         /*  (interpret as 4-byte unsigned ints)   */
    int length,          /* Length of buf array, in bytes          */
                         /*  (must be multiple of 4)               */
    unsigned int sum32) /* 32-bit checksum                        */
{
/*
    Increment the input value of sum32 with the 1's complement sum 
    accumulated over the input buf array.        
*/
    unsigned char *buf = (unsigned char *) buffer;
    unsigned int hi, lo, hicarry, locarry, i;

    /*  Accumulate the sum of the high-order 16 bits and the */
    /*  low-order 16 bits of each 32-bit word, separately.  */
    /*  The first byte in each pair is the most significant. */
    /*  This algorithm works on both big and little endian machines. */

    hi = (sum32 >> 16);
    lo = sum32 & 0xFFFF;
    for (i=0; i < length; i+=4) {
        hi += ((buf[i]   << 8) + buf[i+1]);  
        lo += ((buf[i+2] << 8) + buf[i+3]);
    }

    /* fold carry bits from each 16 bit sum into the other sum */
    hicarry = hi >> 16; 
    locarry = lo >> 16;
    while (hicarry || locarry) {
        hi = (hi & 0xFFFF) + locarry;
        lo = (lo & 0xFFFF) + hicarry;
        hicarry = hi >> 16;
        locarry = lo >> 16;
    }

    /* concatenate the full 32-bit value from the 2 halves */
    sum32 = (hi << 16) + lo;

    return sum32;
}

void char_encode (
     unsigned int value,  /* 1's complement of the checksum value */
                          /*     to be encoded                    */
    char *ascii)          /* Output 16-character encoded string   */
{
    unsigned int exclude[13] = { 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40,
                                 0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60 };
    
    int offset = 0x30;                          /* ASCII 0 (zero) */
    unsigned long mask[4] = { 0xff000000, 0xff0000, 0xff00, 0xff  };
    int byte, quotient, remainder, ch[4], check, i, j, k;
    char asc[32];

    for (i=0; i < 4; i++) {
        /* each byte becomes four */
        byte = (value & mask[i]) >> ((3 - i) * 8);
        quotient = byte / 4 + offset;
        remainder = byte % 4;
        for (j=0; j < 4; j++)
            ch[j] = quotient;

        ch[0] += remainder;

        for (check=1; check;)           /* avoid ASCII punctuation */
            for (check=0, k=0; k < 13; k++)
                for (j=0; j < 4; j+=2)
                    if (ch[j]==exclude[k] || ch[j+1]==exclude[k]) {
                        ch[j]++;
                        ch[j+1]--;
                        check++;
                    }

        for (j=0; j < 4; j++)           /* assign the bytes */
            asc[4*j+i] = ch[j];
    }

    for (i=0; i < 16; i++)              /* permute the bytes for FITS */
        ascii[i] = asc[(i+15)%16];

    ascii[16] = 0;                      /* terminate the string */
}

// very sofisticated pieces of code needed to fool the linker into linking
// statically. I cannot beleive they actually work.
// Actually, I think there is an official name for this kind of shenanigans:
// they are called "weak symbols". 
#ifdef WININTEL
// This one is needed to link statically with acml on Windows
int _imp__putenv() {
  fprintf( stderr, "\nOh my gosh! I am being called!\nThis might end in tears ...\n");
  return 0;
}
// This one is needed so the executable does not need libfcoremd.dll
int for_errsns_load() {
  fprintf( stderr, "\nOh my gosh! I am being called!\nThis might end in tears ...\n");
  return 0;
}
// This one is needed so the executable does not need a certan DLL (do not remeber which one)
int _imp__libc_fpieee_version() {
  fprintf( stderr, "\nOh my gosh! I am being called!\nThis might end in tears ...\n");
  return 0;
}
#endif
// This is for linking statically with the OPENMPI libraries on Linux
int __vdso_clock_gettime() {
  fprintf( stderr, "\nOh my gosh! I am being called!\nThis might end in tears ...\n");
  return 0;
}
