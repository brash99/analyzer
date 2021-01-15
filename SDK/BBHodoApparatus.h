#ifndef ROOT_BBHodoApparatus
#define ROOT_BBHodoApparatus

//////////////////////////////////////////////////////////////////////////
//
// BBHodoApparatus
//
//////////////////////////////////////////////////////////////////////////

#include "THaApparatus.h"

class BBHodoApparatus : public THaApparatus {
  
public:
  BBHodoApparatus( const char* name="U", 
		 const char* description="BigBite Hodoscope Apparatus" );
  virtual ~BBHodoApparatus();

  virtual void  Clear( Option_t* opt="");
  virtual Int_t Reconstruct();

protected:
  
  Int_t  fNtotal;  // Total number of hits in all detectors

  virtual Int_t DefineVariables( EMode mode = kDefine );

  ClassDef(BBHodoApparatus,0) // Hadron Arm
};

#endif

