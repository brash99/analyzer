/////////////////////////////////////////////////////////////////////
//
//  THaCodaFile
//  File of CODA data
//
//  CODA data file, and facilities to open, close, read,
//  write, filter CODA data to disk files, etc.
//
//  This is largely a wrapper around the JLAB EVIO library.
//
//  author  Robert Michaels (rom@jlab.org)
//
/////////////////////////////////////////////////////////////////////

#include "THaCodaFile.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include "evio.h"

using namespace std;

namespace Decoder {

//Constructors

  THaCodaFile::THaCodaFile()
    : ffirst(0), max_to_filt(0), maxflist(0), maxftype(0)
  {
    // Default constructor. Do nothing (must open file separately).
  }

  THaCodaFile::THaCodaFile(const char* fname, const char* readwrite)
    : ffirst(0), max_to_filt(0), maxflist(0), maxftype(0)
  {
    // Standard constructor. Pass read or write flag
    if( codaOpen(fname,readwrite) != CODA_OK )
      fIsGood = false;
  }

  THaCodaFile::~THaCodaFile () {
    //Destructor
    codaClose();
  };

  Int_t THaCodaFile::codaOpen(const char* fname, Int_t mode )
  {
    // Open CODA file 'fname' in read-only mode
    return codaOpen( fname, "r", mode );
  }

  Int_t THaCodaFile::codaOpen(const char* fname, const char* readwrite,
                              Int_t /* mode */ )
  {
    // Open CODA file 'fname' with 'readwrite' access
    init(fname);
    // evOpen really wants char*, so we need to do this safely. (The string
    // _might_ be modified internally ...) Silly, really.
    char *d_fname = strdup(fname), *d_flags = strdup(readwrite);
    Int_t status = evOpen(d_fname,d_flags,&handle);
    fIsGood = (status == S_SUCCESS);
    staterr("open",status);
    free(d_fname); free(d_flags);
    return ReturnCode(status);
  };

  Int_t THaCodaFile::codaClose() {
// Close the file. Do nothing if file not opened.
    if( handle ) {
      Int_t status = evClose(handle);
      handle = 0;
      fIsGood = (status == S_SUCCESS);
      staterr("close",status);
      return ReturnCode(status);
    }
    return ReturnCode(S_SUCCESS);
  }


  Int_t THaCodaFile::codaRead() {
// codaRead: Reads data from file, stored in evbuffer.
// Must be called once per event.
    Int_t status;
    if( handle ) {
      status = evRead(handle, evbuffer, MAXEVLEN);
      fIsGood = (status == S_SUCCESS || status == EOF );
      staterr("read",status);
    } else {
      if(CODA_VERBOSE) {
	 cout << "codaRead ERROR: tried to access a file with handle = 0" << endl;
	 cout << "You need to call codaOpen(filename)" << endl;
	 cout << "or use the constructor with (filename) arg" << endl;
      }
      status = S_EVFILE_BADHANDLE;
    }
    return ReturnCode(status);
  };


  Int_t THaCodaFile::codaWrite(const UInt_t* evbuf) {
// codaWrite: Writes data from 'evbuf' to file
     Int_t status;
     if ( handle ) {
       status = evWrite(handle, evbuf);
       fIsGood = (status == S_SUCCESS);
       staterr("write",status);
     } else {
       cout << "codaWrite ERROR: tried to access file with handle = 0" << endl;
       status = S_EVFILE_BADHANDLE;
     }
     return ReturnCode(status);
   };

  bool THaCodaFile::isOpen() const {
    return (handle!=0);
  }

  Int_t THaCodaFile::filterToFile(const char* output_file) {
// A call to filterToFile filters from present file to output_file
// using filter criteria defined by evtypes, evlist, and max_to_filt
// which are loaded by public methods of this class.  If no conditions
// were loaded, it makes a copy of the input file (i.e. no filtering).

       Int_t i;
       if(filename == output_file) {
	 if(CODA_VERBOSE) {
	   cout << "filterToFile: ERROR: ";
	   cout << "Input and output files cannot be same " << endl;
	   cout << "This is to protect you against overwriting data" << endl;
	 }
	 return CODA_ERROR;
       }
       FILE *fp;
       if ((fp = fopen(output_file,"r")) != NULL) {
	  if(CODA_VERBOSE) {
	    cout << "filterToFile:  ERROR:  ";
	    cout << "Output file `" << output_file << "' exists " << endl;
	    cout << "You must remove it by hand first. " << endl;
	    cout << "This forces you to think and not overwrite data." << endl;
	  }
	  fclose(fp);
	  fIsGood = false;
	  return CODA_FATAL;
       }
       THaCodaFile* fout = new THaCodaFile(output_file,"w");
       if( !fout || !fout->isGood() ) {
         delete fout;
         fIsGood = false;
         return CODA_FATAL;
       }
       Int_t nfilt = 0;

       Int_t status, fout_status = CODA_OK;
       while( (status = codaRead()) == CODA_OK ) {
	   UInt_t* rawbuff = getEvBuffer();
	   Int_t evtype = rawbuff[1]>>16;
	   Int_t evnum = rawbuff[4];
	   Int_t oktofilt = 1;
	   if (CODA_DEBUG) {
	     cout << "Input evtype " << dec << evtype;
	     cout << "  evnum " << evnum << endl;
	     cout << "max_to_filt = " << max_to_filt << endl;
	     cout << "evtype size = " << evtypes[0] << endl;
	     cout << "evlist size = " << evlist[0] << endl;
	   }
	   if ( evtypes[0] > 0 ) {
	       oktofilt = 0;
	       for (i=1; i<=evtypes[0]; i++) {
		   if (evtype == evtypes[i]) {
		       oktofilt = 1;
		       break;
		   }
	       }
	   }
	   if ( evlist[0] > 0 ) {
	       oktofilt = 0;
	       for (i=1; i<=evlist[0]; i++) {
		   if (evnum == evlist[i]) {
		       oktofilt = 1;
		       break;
		   }
	       }
	   }
	   if (oktofilt) {
	     nfilt++;
	     if (CODA_DEBUG) {
	       cout << "Filtering event, nfilt " << dec << nfilt << endl;
	     }
	     fout_status = fout->codaWrite(getEvBuffer());
	     if (fout_status != CODA_OK) {
	       if (CODA_VERBOSE) {
		 cout << "Error in filterToFile ! " << endl;
		 cout << "codaWrite returned status " << fout_status << endl;
	       }
	       break;
	     }
	     if (max_to_filt > 0) {
		if (nfilt == max_to_filt) {
		  break;
		}
	     }
	   }
       }
       if( status == CODA_EOF ) // EOF is normal
         status = CODA_OK;
       fIsGood = (status == CODA_OK);

       fout_status = fout->codaClose();
       delete fout;

       return fIsGood ? fout_status : status;
  };



  void THaCodaFile::addEvTypeFilt(Int_t evtype_to_filt)
// Function to set up filtering by event type
  {
     initFilter();
     if (evtypes[0] >= maxftype-1) {
	TArrayI temp = evtypes;
	maxftype = maxftype + 100;
	evtypes.Set(maxftype);
	for (Int_t i=0; i<=temp[0]; i++) evtypes[i]=temp[i];
	temp.~TArrayI();
     }
     evtypes[0] = evtypes[0] + 1;  // 0th element = num elements in list
     Int_t n = evtypes[0];
     evtypes[n] = evtype_to_filt;
     return;
  };


  void THaCodaFile::addEvListFilt(Int_t event_num_to_filt)
// Function to set up filtering by list of event numbers
  {
     initFilter();
     if (evlist[0] >= maxflist-1) {
	TArrayI temp = evlist;
	maxflist = maxflist + 100;
	evlist.Set(maxflist);
	for (Int_t i=0; i<=temp[0]; i++) evlist[i]=temp[i];
     }
     evlist[0] = evlist[0] + 1;  // 0th element = num elements in list
     Int_t n = evlist[0];
     evlist[n] = event_num_to_filt;
     return;
  };

  void THaCodaFile::setMaxEvFilt(Int_t max_event)
// Function to set up the max number of events to filter
  {
     max_to_filt = max_event;
     return;
  };

  void THaCodaFile::init(const char* fname) {
    if( filename != fname ) {
      codaClose();
      filename = fname;
    }
    handle = 0;
  };

  void THaCodaFile::initFilter() {
    if (!ffirst) {
       ffirst = 1;
       maxflist = 100;
       maxftype = 100;
       evlist.Set(maxflist);
       evtypes.Set(maxftype);
       evlist[0] = 0;
       evtypes[0] = 0;
    }
  };

}

ClassImp(Decoder::THaCodaFile)
