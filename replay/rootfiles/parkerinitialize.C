//
// =================================================================================
// PROGRAM INITIALIZE:
// PURPOSE: THIS PROGRAM RETRIEVES DATA FROM THE ROOT FILE AND SETS THEIR RESPECTIVE
//          BRANCH ADDRESSES. (ADC-TDC..ETC)
//     THIS INITIALZATION IS FOR TEST PURPOSES, WRITTEN BY PETER MONAGHAN AND STUDENTS
// =================================================================================
//

#include "TROOT.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TGraph.h"
#include "parkeranalyze.C"
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
// ok lets find out what is going on!
static const Int_t NUMPMT = 14;
static const Int_t NUMPIXEL = 16;
static const Int_t NUMPADDLE = 14;
static const Int_t NUMPADDLES = NUMPMT*NUMPIXEL; // straightforward

static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const Int_t xcanvas = 800; // width of canvases
static const Int_t ycanvas = 800; // height of canvases

TTree *T; // tree from root file, raw adc/tdc/hit datas
// local arrays to hold data per event from tree
Double_t adc[NUMPADDLES];
Double_t adc_c[NUMPADDLES];
Double_t tdcl[NUMPADDLES];
Double_t tdcl_c[NUMPADDLES];
Double_t tdct[NUMPADDLES];
Double_t tdct_c[NUMPADDLES];
Double_t nahit;
Double_t nthit;
Double_t nhit;

//Latest map
Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13}; // excluded pixel 1's
Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16}; // excluded pixel 2's
//Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 6,13, 4,13,13,13, 5,13};
//Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};
//Original map of Ralph/Tommy
//Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13};
//Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};
//Map to test for missing pixels
//Int_t pixel1[NUMPMT]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//Int_t pixel2[NUMPMT]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t paddleindex[NUMPADDLES];

Int_t run; // run number used in titles of plots
Int_t n_events_to_analyze; // number of events to analyze ... -1 = all.
TStyle *MyStyle = new TStyle("MyStyle","MyStyle");

void parkerinitialize(Int_t runno, Int_t n_events=-1){
 
  TString filename;
  filename.Form("scint_%d.root",runno);
  TFile *_file0 = TFile::Open(filename);

  run=runno; // we lazy now 
// schematics for graph styles, not important right now
  MyStyle->SetTitleFontSize(0.08);
  MyStyle->SetTitleX(0.15);
  MyStyle->SetTitleY(0.99);
  MyStyle->SetStatW(0.9);
  MyStyle->SetMarkerStyle(6);
  gStyle->SetCanvasDefH(xcanvas);
  gStyle->SetCanvasDefW(ycanvas); // this seems backwards? (height and width)
  gStyle->SetPalette(1);
  gROOT->SetStyle("MyStyle");

  T = (TTree *)_file0->Get("T");
  T->SetBranchAddress("C.cdetm1r.adc",&adc);
  T->SetBranchAddress("C.cdetm1r.adc_c",&adc_c);
  T->SetBranchAddress("C.cdetm1r.tdcl",&tdcl);
  T->SetBranchAddress("C.cdetm1r.tdct",&tdct);
  T->SetBranchAddress("C.cdetm1r.tdcl_c",&tdcl_c);
  T->SetBranchAddress("C.cdetm1r.tdct_c",&tdct_c);
  T->SetBranchAddress("C.cdetm1r.nhit",&nhit);
  T->SetBranchAddress("C.cdetm1r.nahit",&nahit);
  T->SetBranchAddress("C.cdetm1r.nthit",&nthit);

  Int_t n_entries = T->GetEntries();  // USEFUL
  cout << "Found " << n_entries << " events"<<endl;
  if (n_events==-1){
	n_events_to_analyze = n_entries;
	cout << "Analyzing all events." << endl;
  } else {
	if (n_events < n_entries) {
		n_events_to_analyze = n_events;
		cout << "Analyzing " << n_events << " events." << endl;
  	}else{
		n_events_to_analyze = n_entries;
		cout << "Analyzing " << n_entries << " events." << endl;
	}
  }

}

void setPaddleIndices(){
	
          for(Int_t pmt=0; pmt<NUMPMT;pmt++){
                Int_t ipaddle = (pmt+1)*NUMPADDLE+1;
                for(Int_t pixel=0; pixel<NUMPIXEL;pixel++){
                        Int_t index = pmt*NUMPIXEL+pixel;
                        if (pixel!=pixel1[pmt]-1&&pixel!=pixel2[pmt]-1){
                                ipaddle--;
				paddleindex[ipaddle]=index;
                        }
                }
           }
	   return;
}

Int_t getPaddleIndex(Int_t pmt, Int_t ipaddle){

	setPaddleIndices();

	Int_t localpaddle = (pmt-1)*NUMPADDLE+ipaddle;

	return paddleindex[localpaddle];
}

Int_t getPaddlePixel(Int_t pmt, Int_t ipaddle){

	setPaddleIndices();

	Int_t localpaddle = (pmt-1)*NUMPADDLE+ipaddle;

	return paddleindex[localpaddle]-(pmt-1)*NUMPIXEL+1;
}

void print_event(Int_t adc_cut=50, Int_t start_event=1, Int_t num_events=1, Int_t tdc_min=750, Int_t tdc_width=300){

        for (Int_t id=start_event;id<start_event+num_events;id++){
          T->GetEntry(id);
	  cout << "Event " << id << endl;
	  for (Int_t pmt=1; pmt<=NUMPMT;pmt++){
	    for (Int_t index=1; index<=NUMPIXEL; index++){
	      Int_t ipixel = (pmt-1)*NUMPIXEL+index-1;
	      if(tdcl[ipixel]>tdc_min&&tdcl[ipixel]<tdc_min+tdc_width){
		cout << "Leading Edge TDC Hit on global pixel " << ipixel << " = " << tdcl[ipixel] << endl;
	      }
	      if(tdct[ipixel]>tdc_min&&tdct[ipixel]<tdc_min+tdc_width){
		cout << "Trailing Edge TDC Hit on global pixel " << ipixel << " = " << tdct[ipixel] << endl;
	      }
	      if(adc_c[ipixel]>adc_cut){
		cout << "ADC Hit on global pixel " << ipixel << " = " << adc_c[ipixel] << endl;
	      }	
	    }
	  }
	}
	
	return;
}
//
//
//    OK: Modifying this portion to do my analysis later on. Will copy the useful parts out
//    and annotate as I go along
//=======================================================
//  "METHOD" PLOT_ADC
//   Written by Peter Monaghan and CNU students
//   This program takes the ADC data and stores the pedestal and
//   signal separately in htmpa and htmpb respectively
//=======================================================
/*TCanvas *plot_adc(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADC",run);
        TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=600;  // arbitrary bins, causing problems sometimes
        Int_t min=-100, max=500;
        for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=n_events_to_analyze; // I think this part is what actually fills the histograms
// There is a distinction made here with the ADC's being made with their neighboring TDC's
        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id); // grabs a single entry
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(adc_c[ipaddle]);
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpb[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}
// ======================================================================
//
//  Analysis of tail end of adc.
//  Notes: The ADC signal has already been determined by the tdc analysis
//  

// cout<<"Max bin: "<< htmpa1->GetMaximumBin();


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//  This is where the loop will begin
//  note that this needs to loop from inputs adc(1-14), it may be reasonable to store
//  data in a file for each iteration and then move onto the next.
//  so ideally
//   loop of adc()
//   
//
//
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

 TF1 *myfit = new TF1("myfit", "gaus"); // line fit



   myfit->SetParName(0,"Height of Gaus");
    myfit->SetParName(1,"not used");
    myfit->SetParName(2,"sig"); // this is the sigma parameter which we need
//
//
 cout << " we made it this far ";


 // Double_t binmax = htmpa[1]->GetMaximumBin(); // need to work on this binmax, 
htmpa[1]->Fit("myfit","","", -50,50);
//		function = htmpb[1]->GetFunction("gaus");
//		function->SetLineColor(1);

//		Double_t constants = function->GetParameter(0);
//		Double_t means = function->GetParameter(1);
		Double_t sigmas = myfit->GetParameter(2);
                
 cout << "This is sigma of the gaussian yo " << sigmas;





//
// ======================================================================
/*	cADC->Clear();
        cADC->Divide(4,4);

        //plot histos
        Int_t icount=0;
        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cADC->cd( icount + 1 );
             gPad->SetLogy();

            cADC->Update();
                
            Int_t entries = htmpa[i]->GetEntries();
            float mean = htmpa[i]->GetMean(1);
            float RMS = htmpa[i]->GetRMS(1);
            //cout << entries <<" "<< mean <<" "<< RMS <<endl;
                
            htmpa[i]->SetStats(0);
            //current->Modified();
            htmpa[i]->Draw();
            htmpb[i]->Draw("same");

            icount++;
          
	  }
        }

        title.Form("run_%d_ADC_pmt_%d_tdc_min_%d_max_%d.png",
                   run,pmt,tdc_min,tdc_min+tdc_width);
        cADC->Print(title);
        cADC->cd(0); */
   //     return cADC;
//}

