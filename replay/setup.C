#include "TString.h"
#include <cstring>
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaGlobals.h"

#include "../SDK/CdetApparatus.h"

int setup()
{
  //  R. Michaels, Jan 2017
  //  Steering script for Hall A analyzer

  //gSystem->Load("./libSBS.so");

  //THaApparatus* SBSFb = new FbusApparatus("S","SBS Fastbus Test Stand");
  //gHaApps->Add( SBSFb );

  //THaApparatus* HadArm = new HcalApparatus("H","SBS Hadron Arm HCAL");
  //gHaApps->Add( HadArm );
  
  THaApparatus* Cdet = new CdetApparatus("C","Cdet Detector");
  gHaApps->Add( Cdet );
  
  int run_number = 1089;
  TString experiment = "scint";

  // Get the run number
  cout << "Run number (hint: 1089, 1090, ...)? ";
  cin >> run_number;
  
  // Check if the input file exists
  TString data_dir = gSystem->Getenv("DATA_DIR");
  if( data_dir.IsNull() )
    data_dir = "../../data";
  TString run_file = data_dir + "/" + experiment + Form("_%d.dat",run_number);
  if( gSystem->AccessPathName(run_file) ) {
  	Error("setup.C", "Input file does not exist: %s", run_file.Data() );
  	return 1;
  }
  
  THaRun* run = new THaRun( run_file, "CDet Cosmics Data" );
  run->SetDataRequired( THaRunBase::kDate );

  int nev = -1;
  cout << "Number of events to replay (-1=all)? ";
  cin >> nev;
  if( nev > 0 )
    run->SetLastEvent(nev);

  THaEvent* event = new THaEvent;
  
  THaAnalyzer* analyzer = new THaAnalyzer; 

  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = "./rootfiles";
  TString out_file = out_dir + "/" + experiment + Form("_%d.root", run_number);

  analyzer->SetOutFile( out_file );
 
  analyzer->EnableBenchmarks();
  analyzer->SetEvent( event );
  analyzer->SetCutFile("onlana.cuts");
  analyzer->SetOdefFile("output.def");

  analyzer->Process(run);

  analyzer->Close();
  delete analyzer;
  gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();

  TFile* rootfile = new TFile( out_file, "READ" );

  return 0;
}
