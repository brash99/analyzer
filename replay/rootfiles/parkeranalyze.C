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

#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

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
double all[15][15];
 ofstream myfile;
  myfile.open ("Integral.txt");
TCanvas *check_signal(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADC",run);
        TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);


        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=200;  // arbitrary bins, causing problems sometimes
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
  double sigmas[15];
// double 3sig[15];
double integral[15];
double total [15];
   myfit->SetParName(0,"Height of Gaus");
    myfit->SetParName(1,"not used");
    myfit->SetParName(2,"sig"); // this is the sigma parameter which we need
//
//
 cout << " we made it this far ";

   for (int i = 1; i<15; i++){
 // Double_t binmax = htmpa[1]->GetMaximumBin(); // need to work on this binmax, 
htmpa[i]->Fit("myfit","","", -100,100);
//		function = htmpb[1]->GetFunction("gaus");
//		function->SetLineColor(1);

//		Double_t constants = function->GetParameter(0);
//		Double_t means = function->GetParameter(1);
		 sigmas[i] = myfit->GetParameter(2);
                 sigmas[i]  = 3*sigmas[i];
                
 cout << "This is sigma of the gaussian " << sigmas[i] << endl;


}
// integrate over the entire thing first to check event
   for (int i = 1; i<15; i++){
  TAxis *axis = htmpa[i]->GetXaxis();
  int bmin = axis->FindBin(-100); //in your case xmin=-1.5
  int bmax = axis->FindBin(500); //in your case xmax=0.8
   total[i] = htmpa[i]->Integral(bmin,bmax);

}

// now we integrate the rest from 3 sigma
   for (int i = 1; i<15; i++){
  TAxis *axis = htmpa[i]->GetXaxis();
  int bmin = axis->FindBin(sigmas[i]); 
  int bmax = axis->FindBin(500); 
   integral[i] = htmpa[i]->Integral(bmin,bmax);
cout << "PMT: " <<pmt <<" PADDLE: " << i << "\t 3 SIGMA:  "<< setprecision(3) << sigmas[i] <<" \t EVENTS AFTER 3 SIGMA: "<< setprecision(5) << integral[i] <<"\t TOTAL EVENTS: "<< setprecision(7) << total[i] << endl;

// check for suffiecient number of "non pedestal" events
   	if(integral[i]/total[i] < 0.001){ cout << "paddle " << i << " does not exceed 0.1 percent of events "<< endl;
}
// check for too many "non pedestal" events
	if(integral[i]/total[i] > 0.05) {
cout << "paddle " << i << " exceeds 5 percent of events! " << endl;
// add to total array
  all[pmt][i] = integral[i];
  myfile <<  all[pmt][i] << endl;

//	gr3 = new TGraph(14,pmt,integral[i]);
//	gr3->SetLineColor(3);
//	gr3->SetLineWidth(4);

//	gr->Draw("AP");
   if (pmt = 14){myfile.close();}
}

}
// ======================================================================

// 
	cADC->Clear();
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
        cADC->cd(0); 
        return cADC;
}
TCanvas *cross_pixel(Int_t pmt=1, Int_t pix1 = 1 , Int_t tdc_min=750, Int_t tdc_width=300){
//
//
//
//
//
//
//
//
// ----------------------------------------------------------------------
//   TO DO:
// Now need to write a script that compares the two adc values of pmt 9 to one another, to determine
// whether there is a clear relation between it's gradual increase/decrease of pixel width and pickup
// from the giant adc signal "solenoids" located in the test lab.
// 
        TString cut, draw, draw1, title;
        title.Form("run_%d_ADC",run);
        TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);

	TH1D *htmpb[NUMPIXEL];// 
        TH1D *htmpa[NUMPIXEL];// 
	TH2D *cross1[NUMPIXEL];// this will be used to check pixel correlations


      TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=200;  // arbitrary bins, causing problems sometimes
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
   cross1[icounter - 1] = new TH2D("cross","cross",nbin,min,max,nbin,min,max); // 2d histogram initalization,
        }

        Int_t nentries=n_events_to_analyze; // This part is what actually fills the histograms
// There is a distinction made here with the ADC's being made with their neighboring TDC's
        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id); // grabs a single entry. ok, now what
          
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(adc_c[ipaddle]);
            cross1[pix1-1] -> Fill(adc_c[ipaddle],adc_c[ipaddle+1])
            
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpb[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}
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


            icount++; }}// questionable close braces
// test gaussian
   for (int n = 0; n<10000; n++){
cross[pix1]->Fill( gRandom->Gaus(0.2,0.2), gRandom->Gaus(1.4,0.2) );}
            cross1[pix1]->Draw();
        title.Form("run_%d_ADC_pmt_%d_tdc_min_%d_max_%d.png",
                   run,pmt,tdc_min,tdc_min+tdc_width);
        cADC->Print(title);
        cADC->cd(0); 
        return cADC;
//    Take cross1[i] and cross2[i+1] to compare each event. 
//
//
//    each event has an idex from 0 - NEvents. so we grab each event and compare
//    amplitudes, and plot against one another. 
//
//    if there is a pickup relation we should see a correlation on the 2d plot.
//
//   perform a linear fit and look for the coffecient of regression.
//
// ----------------------------------------------------------------------
//
TCanvas *ADC_Cut_Fit(Int_t pmt=1, Int_t ADClmax=10, Int_t ADCrmax=10){
Int_t indx = (pmt-1)*NUMPIXEL // sets the starting point
        TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);
        cADC -> Divide(4,4);
	TH1D *ADC[NUMPIXEL];// fix this later.
	 

return cADC
}


}
