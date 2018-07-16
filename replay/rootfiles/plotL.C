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

using namespace std;

static const Int_t NUMPMT = 14;
static const Int_t NUMPIXEL = 16;
static const Int_t NUMPADDLE = 14;
static const Int_t NUMPADDLES = NUMPMT*NUMPIXEL;

static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const Int_t xcanvas = 800; // width of canvases
static const Int_t ycanvas = 800; // height of canvases

TTree *T; // tree from root file, raw adc/tdc/hit data
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

//Latest map set for M1-L
Int_t pixel1[NUMPMT]={1, 2, 13,4, 6, 12, 1,1, 1,5,13,1, 1,8};
Int_t pixel2[NUMPMT]={14, 16,16,6,13,13,13,3,6,9,15,4,16,12};
//Map for M1-R
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

void plotL(Int_t runno, Int_t n_events=-1){
 
  TString filename;
  filename.Form("scint_%d.root",runno);
  TFile *_file0 = TFile::Open(filename);

  run=runno;
// Set Canvas preferences
  MyStyle->SetTitleFontSize(0.08);
  MyStyle->SetTitleX(0.15);
  MyStyle->SetTitleY(0.99);
  MyStyle->SetStatW(0.9);
  MyStyle->SetMarkerStyle(6);
  gStyle->SetCanvasDefH(xcanvas);
  gStyle->SetCanvasDefW(ycanvas);
  gStyle->SetPalette(1);
  gROOT->SetStyle("MyStyle");

  T = (TTree *)_file0->Get("T"); 
// Setting addresses so we can call them later
  T->SetBranchAddress("C.cdetm1r.adc",&adc);
  T->SetBranchAddress("C.cdetm1r.adc_c",&adc_c);
  T->SetBranchAddress("C.cdetm1r.tdcl",&tdcl);
  T->SetBranchAddress("C.cdetm1r.tdct",&tdct);
  T->SetBranchAddress("C.cdetm1r.tdcl_c",&tdcl_c);
  T->SetBranchAddress("C.cdetm1r.tdct_c",&tdct_c);
  T->SetBranchAddress("C.cdetm1r.nhit",&nhit);
  T->SetBranchAddress("C.cdetm1r.nahit",&nahit);
  T->SetBranchAddress("C.cdetm1r.nthit",&nthit);

  Int_t n_entries = T->GetEntries(); // Checks how many total entries in T
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
	// Setting the location of paddles, taking into account the missing pixels in each PMT
          for(Int_t pmt=0; pmt<NUMPMT;pmt++){
	    Int_t ipaddle = (pmt)*NUMPADDLE;
                for(Int_t pixel=0; pixel<NUMPIXEL;pixel++){
		  Int_t index = (pmt*NUMPIXEL)+pixel;
                        if (pixel!=pixel1[pmt]-1&&pixel!=pixel2[pmt]-1){
                                ipaddle++;
				paddleindex[ipaddle]=index;
                        }
                }
           }
	   return;
}



Int_t getPaddleIndex(Int_t pmt, Int_t ipaddle){
  //pmt numbers are from 0 to 13 as are ipaddle numbers and the code accounts for removed pixels
	setPaddleIndices();

	Int_t localpaddle = (pmt+1)*NUMPADDLE-ipaddle;

	return paddleindex[localpaddle];
}

Int_t getPaddlePixel(Int_t pmt, Int_t ipaddle){
  //gives the pixel number for the given paddle number of the pmt, again paddle and pmt numbers are from 0 to 13
	setPaddleIndices();

	Int_t localpaddle = (pmt+1)*NUMPADDLE-ipaddle;

	return paddleindex[localpaddle]-(pmt)*NUMPIXEL;
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

TCanvas *plot_xtalk(Int_t pmt=1, Int_t highv=800, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_TDC",run);
        TCanvas *cXTalk= new TCanvas("cXTalk",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

	Int_t adc_cut = (highv-650)/(10);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=tdc_width;
        Int_t min=tdc_min, max=tdc_min+tdc_width;
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

        Int_t nentries = n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
           T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    if(tdcl[ipaddle]>500)
	      {
		htmpa[index-1]->Fill(tdcl[ipaddle]);
		if(adc_c[ipaddle]<adc_cut)
		  {
		    htmpb[index-1]->Fill(tdcl[ipaddle]);
		  }
	      }
	  }
	}

	cXTalk->Clear();
        cXTalk->Divide(4,4);

        //plot histos
        Int_t icount=0;
        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cXTalk->cd( icount + 1 );
            //gPad->SetLogy();

            cXTalk->Update();
                
            Int_t entriesa = htmpa[i]->GetEntries();
            float meana = htmpa[i]->GetMean(1);
            float RMSa = htmpa[i]->GetRMS(1);

	    Int_t entriesb = htmpb[i]->GetEntries();
            float meanb = htmpb[i]->GetMean(1);
            float RMSb = htmpb[i]->GetRMS(1);

            cout <<"Pixel number: "<< i+1<<"   \t entries = "<< entriesa <<"   \t entries = "<< entriesb<<endl;
                
            htmpa[i]->SetStats(0);
            //current->Modified();
            htmpa[i]->Draw();
            htmpb[i]->Draw("same");

            icount++;
          
	  }
        }

        title.Form("run_%d_TDC_pmt_%d_adc_cut_%d.png",
                   run,pmt,adc_cut);
        cXTalk->Print(title);
        cXTalk->cd(0);
        return cXTalk;
}

TCanvas *plot_adc_fit(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=400, Int_t adc_neighbor_cut=20, Int_t adc_cut=50){

        TString cut, draw, draw1, title, grtitle;
        title.Form("run_%d_ADC_Fit",run);
        TCanvas *cADCFit = new TCanvas("cADCFit",title,xcanvas,ycanvas);
	title.Form("run_%d_ADC_Mean_Fit",run);
	TCanvas *cADCMeanFit= new TCanvas("cADCMeanFit",title,xcanvas,ycanvas);

	// Setting up different arrays and objects used to create the two canvases.

	TH1D *htmp[NUMPIXEL];
	TF1 *function;
	TF1 *function1;
	Int_t histmax_cut[16], gaus_cut_plus[16];
	Float_t gaus_cut_minus[16];

	Double_t upixel[NUMPIXEL]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	Double_t epixel[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t errors[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t constants[16], means[16], sigmas[16];

	TString tmpentry;
	MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatY(0.4);

	setPaddleIndices(); //Setting the geometric paddle locations

	Int_t nbin=90;
	Int_t min=-100, max=900; // change max to 1500 for PMT 9

	for(Int_t i = 1; i <= NUMPIXEL; i++)
	  {
	    tmpentry.Form("htmp%d",i);
	    htmp[i-1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
	    htmp[i-1]->SetLineColor(kRed);
	    title.Form("Run %d ADC Fit pmt %d, paddle %d: %d < tdc < %d",run,pmt,i,tdc_min,tdc_min+tdc_width);
	    htmp[i-1]->SetTitle(title);
	  }

	// Filling the histograms using only adc data with good tdc and good adc neighbor cut.

	Int_t nentries=n_events_to_analyze;

       	for (Int_t id=1;id<=nentries;id++)
	  {
	    T->GetEntry(id);
	    Int_t ipaddle = (pmt)*NUMPADDLE-NUMPADDLE;
	    for (Int_t pixel=0; pixel < NUMPIXEL; pixel++)
	    {
		Int_t index = (pmt-1)*NUMPIXEL+pixel;
		if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1)
		{   
		  ipaddle++;
		    if (ipaddle < 2) 
		    {
			if (adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut)
			{
			  htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		    else if (ipaddle > NUMPMT*NUMPADDLE)
		    {
			if (adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut)
			{
			    htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		    else
		    {
			if (adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut && adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut)
			{
			  htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		}
	    }
	  }


	// Creating the canvas of adc data with good tdc and fitting the adc with a gaussian or "landau" function.

	cADCFit->Clear();
	cADCMeanFit->Clear();
	cADCFit->Divide(4,4);

        Int_t count = 0;
        for (Int_t i=0; i<NUMPIXEL; i++)
	  {
	    if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) 
	      {
         //       means[i]=0.0;
	//	int check = 0;
		cADCFit->cd(count+1);
		 gPad->SetLogy();
		cADCFit->Update();
                
		Int_t entries = htmp[i]->GetEntries();
		//float mean = htmp[i]->GetMean(1);
		//float RMS = htmp[i]->GetRMS(1);
                
		htmp[i]->SetStats(0);
		htmp[i]->Draw();

		Double_t bc =0.0;
		Double_t bn = 0.0;

	    	for (Int_t j=15; j<100 ;j++){
		  bc = htmp[i]->GetBinContent(j);
		  if( bc == 0.0 || (bc <htmp[i]->GetBinContent(j+1) && bc < htmp[i]->GetBinContent(j+2) && bc <htmp[i]->GetBinContent(j-1) && bc <htmp[i]->GetBinContent(j-2) && bc <htmp[i]->GetBinContent(j-3) && bc <htmp[i]->GetBinContent(j-4))  ){
		    bn = htmp[i]->GetBinCenter(j);
		    htmp[i]->Fit("gaus","","",bn,bn+600);
		    break;
		  }
		}       

		//htmp[i]->Fit("landau","","", 0, 250);
		function = htmp[i]->GetFunction("gaus");
		function->SetLineColor(1);

		constants[i] = function->GetParameter(0);
		means[i] = function->GetParameter(1);
		sigmas[i] = function->GetParameter(2);
              
	count++;
	      
	      }
	  }
	  

	// Printing out means and sigmas to screen with corresponding pixel number.

	for (Int_t i=0; i<NUMPIXEL; i++)
	  {
	    cout << "Pixel number: " << i+1 << "   \t mean = " << means[i] << "    \t sigma = " << sigmas[i] << endl;
	  }

        title.Form("run_%d_ADC_pmt_%d_tdc_min_%d_max_%d.png",run,pmt,tdc_min,tdc_min+tdc_width);
        cADCFit->Print(title);
        cADCFit->cd(0);

	// End of work on canvas with adc data fitted with a function

	// Creating the canvas of mean adc data with good tdc data.

	cADCMeanFit->cd();

	Double_t mmean[2]={0,0};
	Double_t merror[2]={1,1};
	Double_t mpixel[2]={0,0};
	Double_t mepixel[2]={0,0};
	mpixel[0]=pixel1[pmt-1];
	mpixel[1]=pixel2[pmt-1];

	Double_t yline[2]={100,100};
	Double_t xline[2]={0,17};

	// Calculating errors on the mean adc data and printing out number of entries, means, and errors
	// to screen with corresponding pixel number.

	count = 0;
	for (Int_t i = 0; i < NUMPIXEL; i++)
	{
	    if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) 
	    {
		Int_t entries = htmp[i]->GetEntries();
		errors[i] = sigmas[i]/sqrt(entries);
		if(means[i] > 1)
		{
		    cout << "Pixel number: " <<  i+1 << "  \t entries: " << entries << "\t mean: " 
			 << means[i] << "\t error: " << errors[i] << endl;
		}
               	count++;
	    }
	}

	// Creating three different graphs for used pixels, missing pixels, and general formatting.

	TGraphErrors *gr = new TGraphErrors(NUMPIXEL,upixel,means,epixel,errors);
	gr->SetMarkerStyle(21);
	gr->GetXaxis()->SetTitle("Pixel Number");
	gr->GetYaxis()->SetTitle("Mean ADC Fit (Good TDC)");
	gr->GetYaxis()->SetTitleOffset(1.4);
	grtitle.Form("run_%d_pmt_%d_adc_mean",run,pmt);
	gr->SetTitle(grtitle);
	TGraphErrors *gr2 = new TGraphErrors(2,mpixel,mmean,mepixel,merror);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(2);
	gr2->SetTitle("");
	TGraph *gr3 = new TGraph(2,xline,yline);
	gr3->SetLineColor(2);
	gr3->SetLineWidth(2);
	gr3->SetLineStyle(2);
	gr3->SetTitle("");

	// Drawing the three different graphs to the adc mean values canvas.

	gr->Draw("AP");
	gr2->Draw("P");
	gr3->Draw("L");
            
        cADCMeanFit->Update();
	cADCMeanFit->cd(0);

	// End of work on canvas with adc mean data.

	// Creating a file to output parameters of function fitted to the ADC data with good TDC.
        TString filename;
        filename.Form("run_%d_pmt_%d_adc_fitresults.txt", run, pmt);
	ofstream outfile(filename);

	outfile << "Run #: " << run << "\t PMT #: " << pmt << " \t  Paramters of fitted ADC data with good TDC. \n";
	outfile << "Pixel # \t Mean Value \t Error \t Sigma \n" << endl;
	for(Int_t i = 0; i < NUMPIXEL; i++)
	  {
	    outfile << i+1 << "\t" << means[i] << "\t" << errors[i] << "\t" << sigmas[i] << endl;
	  }
	cout << "Created output file for parameters." << endl;
	outfile.close();
	
	return cADCMeanFit;
	return cADCFit;
}

TCanvas *plot_tdc_fit(Int_t pmt=1, Int_t pixel=1, Int_t highv=800, Int_t tdc_cut=500, Int_t adc_neighbor_cut=20, 
		      Int_t tdc_neighbor_cut=200){

        TString cut, draw, draw1, title;
        title.Form("run_%d_TDC_Fit",run);
        TCanvas *cTDCFit= new TCanvas("cTDCFit",title,xcanvas,ycanvas);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

	TF1 *f1 = new TF1("f1","pol1",700,1100);
	f1->SetLineColor(2);
	f1->SetLineWidth(2);
	f1->SetLineStyle(9);

	Int_t nentries = n_events_to_analyze;
	Int_t adc_cut = (highv-650)/(10);

	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleSize(0.05,"x");
	gStyle->SetTitleSize(0.05,"a");
	gStyle->SetTitleOffset(1.5,"a");
	gStyle->SetTitleOffset(0.9,"x");
	gStyle->SetTitleOffset(1.5,"y");
	gStyle->SetPadColor(0);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(1);
	gStyle->SetStatColor(10); 
	gStyle->SetStatH(0.16);
	gStyle->SetStatW(0.16);
	gStyle->UseCurrentStyle();

	TGraphErrors *g1 = new TGraphErrors();
	g1->SetName("g1");
	g1->SetMarkerStyle(20);       
	g1->SetMarkerColor(1);   
	g1->SetLineColor(kRed);
	g1->SetMarkerSize(0.5);

	TGraphErrors *g2 = new TGraphErrors();
	g2->SetName("g2");
	g2->SetMarkerStyle(20);       
	g2->SetMarkerColor(1);   
	g2->SetLineColor(kRed);
	g2->SetMarkerSize(0.5);

	TGraphErrors *g3 = new TGraphErrors();
	g3->SetName("g3");
	g3->SetMarkerStyle(20);       
	g3->SetMarkerColor(3);   
	g3->SetLineColor(kRed);
	g3->SetMarkerSize(0.5);

	cTDCFit->Clear();
     	cTDCFit->cd(0);

	Int_t shift = -23;
	Int_t paddle = (pmt-1)*NUMPIXEL+pixel-1;
	Int_t count1 = 0;
	Int_t count2 = 0;
	Int_t count3 = 0;
	Double_t diff = 0;
	Double_t slope = 0; 
	Double_t ymin = 0;
	Double_t xmax = 0;
	Double_t tdc_width_cut = 0;
	
	for (Int_t id=1;id<=nentries;id++){
	  T->GetEntry(id);
	  if(tdcl[paddle]>500)
	    {
	     g1->SetPoint(count1,tdct[paddle],tdcl[paddle]);
	     g1->SetPointError(count1,0.0,0.0);
	     if(tdcl[paddle]<880 || tdct[paddle]>900)
		//&& adc_c[paddle+1]<adc_neighbor_cut && adc_c[paddle-1]<adc_neighbor_cut && tdcl[paddle+1]<200
		//&& tdcl[paddle-1]<200 && tdcl[paddle+2]<200 && tdcl[paddle+3]<200 && tdcl[paddle+4]<200
		//&& tdcl[paddle+5]<200 && tdcl[paddle+6]<200 && tdcl[paddle+7]<200 && tdcl[paddle+8]<200
		//&& tdcl[paddle+9]<200 && tdcl[paddle+10]<200 && tdcl[paddle+11]<200 && tdcl[paddle+12]<200
		//&& tdcl[paddle+13]<200 && tdcl[paddle+14]<200 && tdcl[paddle+15]<200)
	       {
		 g2->SetPoint(count2,tdct[paddle],tdcl[paddle]);
		 g2->SetPointError(count2,0.0,0.0);
		 count2++;
	       }
	     count1++;
	    }
	}

	TString histo_title;
	histo_title.Form("Comparison of TDCL vs. TDCT for PMT %d Pixel %d",pmt,pixel);
	TH1F *mdf = new TH1F("mdf",histo_title,1050,750,1050);
	mdf->SetMaximum(1050.0);
	mdf->SetMinimum(750.0);
	mdf->SetStats(0);

	mdf->GetYaxis()->SetTitle("TDC Leading Value");
	mdf->GetXaxis()->SetTitle("TDC Trailing Value");
	mdf->SetTitleOffset(1.6,"y");
	mdf->Draw();

	g1->Draw("same P");
	
	if(count2>5){
	  g2->Draw("same P");
	  g2->Fit("f1","B");
	  Double_t diff = g2->GetFunction("f1")->GetParameter(0);
	  Double_t slope = g2->GetFunction("f1")->GetParameter(1);

	  Double_t ymin = 750+diff;
	  Double_t xmax = 750+((1020-ymin)/(slope));
	  Double_t tdc_width_cut = diff+shift;

	  TLine *lg = new TLine(750.0,ymin+shift,xmax,1020+shift);
	  lg->SetLineColor(8);
	  lg->SetLineStyle(2);
	  lg->SetLineWidth(3);
	  lg->Draw();

	  for (Int_t id=1;id<=nentries;id++){
	    T->GetEntry(id);
	    if(tdcl[paddle]>500 && (tdcl[paddle]-tdct[paddle])>tdc_width_cut)
	      {
		g3->SetPoint(count3,tdct[paddle],tdcl[paddle]);
		g3->SetPointError(count3,0.0,0.0);
		count3++;
	      }
	  }
	
	  g3->Draw("same P");
	}

	Int_t tdccount = 0;
	Int_t tdcwcount = 0;

	for (Int_t id=1;id<=nentries;id++){
	  T->GetEntry(id);
	  if(tdcl[paddle]>500)
	    {
	      if(adc_c[paddle]<adc_cut && (tdcl[paddle]-tdct[paddle])>tdc_width_cut)
	       {
		 tdcwcount++;
	       }
	     tdccount++;
	    }
	}

	cout << "PMT " << pmt << " Pixel " << pixel << "\t " << tdcwcount << "\t" << tdccount << "\t" <<
	  "TDC Width Cut: " << tdc_width_cut << endl;

        title.Form("run_%d_TDC_pmt_%d.png",run,pmt);
        cTDCFit->Print(title);
        return cTDCFit;
}
  

TCanvas *plot_adc(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADC",run);
        TCanvas *cADC= new TCanvas("cADC",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=600;
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

        Int_t nentries=n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(adc_c[ipaddle]);
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpb[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}

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

TCanvas *plot_tdc(Int_t pmt=1, Int_t adc_cut=40, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_TDC",run);
        TCanvas *cTDC= new TCanvas("cTDC",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=tdc_width;
        Int_t min=tdc_min, max=tdc_min+tdc_width;
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

        Int_t nentries = n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(tdcl[ipaddle]);
	    if(adc_c[ipaddle]>adc_cut){
		htmpb[index-1]->Fill(tdcl[ipaddle]);
	    }
	  }
	}

	cTDC->Clear();
        cTDC->Divide(4,4);

        //plot histos
        Int_t icount=0;
        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cTDC->cd( icount + 1 );
            //gPad->SetLogy();

            cTDC->Update();
                
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

        title.Form("run_%d_TDC_pmt_%d_adc_cut_%d.png",
                   run,pmt,adc_cut);
        cTDC->Print(title);
        cTDC->cd(0);
        return cTDC;
}

TCanvas *plot_tdc_adc(Int_t pmt=1, Int_t tdc_min=1, Int_t tdc_width=100){

        TString cut, draw, draw1, title;
        title.Form("run_%d_TDC",run);
        TCanvas *cTDCA= new TCanvas("cTDCA",title,xcanvas,ycanvas);

        TH2D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbinx=tdc_width;
        Int_t minx=tdc_min, maxx=tdc_min+tdc_width;
        Int_t nbiny=350;
        Int_t miny=0, maxy=350;
        for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH2D(tmpentry,tmpentry,nbinx,minx,maxx,nbiny,miny,maxy);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
        }

        Int_t nentries = n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(tdcl[ipaddle]-tdct[ipaddle],adc_c[ipaddle]);
	  }
	}

	cTDCA->Clear();
        cTDCA->Divide(4,4);

        //plot histos
        Int_t icount=0;
        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cTDCA->cd( icount + 1 );
            //gPad->SetLogy();

            cTDCA->Update();
                
            Int_t entries = htmpa[i]->GetEntries();
            float mean = htmpa[i]->GetMean(1);
            float RMS = htmpa[i]->GetRMS(1);
            //cout << entries <<" "<< mean <<" "<< RMS <<endl;
                
            htmpa[i]->SetStats(0);
            //current->Modified();
            htmpa[i]->Draw("COLZ");

            icount++;
          
	  }
        }

        title.Form("run_%d_TDC_pmt_%d.png",
                   run,pmt);
        cTDCA->Print(title);
        cTDCA->cd(0);
        return cTDCA;
}

TCanvas *plot_tdc_diff(Int_t pmt=1, Int_t adc_cut=30, Int_t tdc_min=1, Int_t tdc_width=100){

        TString cut, draw, draw1, title;
        title.Form("run_%d_TDC",run);
        TCanvas *cTDCT= new TCanvas("cTDCT",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=tdc_width;
        Int_t min=tdc_min, max=tdc_min+tdc_width;
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

        Int_t nentries = n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(tdcl[ipaddle]-tdct[ipaddle]);
	    if(adc_c[ipaddle]>adc_cut){
		htmpb[index-1]->Fill(tdcl[ipaddle]-tdct[ipaddle]);
	    }
	  }
	}

	cTDCT->Clear();
        cTDCT->Divide(4,4);

        //plot histos
        Int_t icount=0;
        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cTDCT->cd( icount + 1 );
            //gPad->SetLogy();

            cTDCT->Update();
                
            Int_t entriesa = htmpa[i]->GetEntries();
            float meana = htmpa[i]->GetMean(1);
            float RMSa = htmpa[i]->GetRMS(1);
            //cout << "No ADC Cut:" << "\t" <<  entriesa << "\t" << meana << "\t" << RMSa <<endl;
                
            htmpa[i]->SetStats(0);
            //current->Modified();
            htmpa[i]->Draw();
            htmpb[i]->Draw("same");

	    Int_t entriesb = htmpb[i]->GetEntries();
	    float meanb = htmpb[i]->GetMean(1);
	    float RMSb = htmpb[i]->GetRMS(1);
	    //cout << "With ADC Cut:" << entriesb <<" "<< meanb <<" "<< RMSb <<endl;

            icount++;
          
	  }
        }

        title.Form("run_%d_TDC_pmt_%d_adc_cut_%d.png",
                   run,pmt,adc_cut);
        cTDCT->Print(title);
        cTDCT->cd(0);
        return cTDCT;
}

TCanvas *plot_occupancy(Int_t adc_cut=40, Int_t adc_neighbor_cut=10000, Int_t multiplicity_cut=12, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_Occupancy",run);

        Int_t nbin=196;
        Int_t min=1, max=197;
        Int_t nbinm=11;
        Int_t minm=-1, maxm=10;
        TH1D *hoccupancy = new TH1D("hoccupancy","hoccupancy",nbin,min,max);
        TH1D *hmultiplicity = new TH1D("hmultiplicity","hmultiplicity",nbinm,minm,maxm);
        TH2D *heatmap = new TH2D("heatmap","heatmap",nbin,min,max,nbinm,minm,maxm);

        TCanvas *cOCCUPANCY= new TCanvas("cOCCUPANCY",title,xcanvas,ycanvas);

        Int_t nentries=n_events_to_analyze;
        
	TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.9);
        MyStyle->SetStatW(0.2);

	setPaddleIndices();

        for (Int_t id=1;id<=nentries;id++){
	  T->GetEntry(id);
	 
          Int_t nmultiplicity=0;
          Int_t good_paddle[500];
          for(Int_t icount=0;icount<500;icount++){good_paddle[icount]=-1;}

          for(Int_t pmt=0; pmt<NUMPMT;pmt++){
		Int_t ipaddle = (pmt)*NUMPADDLE;
		for(Int_t pixel=0; pixel<NUMPIXEL;pixel++){
			Int_t index = pmt*NUMPIXEL+pixel;
			if (pixel!=pixel1[pmt]-1&&pixel!=pixel2[pmt]-1){
				ipaddle++;
	    			if (tdcl[index]>tdc_min&&tdcl[index]<tdc_min+tdc_width){
	      			  nmultiplicity++;
				  if (ipaddle == 1){
	      			   if (adc_c[index]>adc_cut && adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut ){
					good_paddle[nmultiplicity-1]=ipaddle;
	      			   }
				  }else if (ipaddle == NUMPMT*NUMPADDLE+1){
	      			   if (adc_c[index]>adc_cut && adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut ){
					good_paddle[nmultiplicity-1]=ipaddle;
	      			   }
				  }else{
	      			   if (adc_c[index]>adc_cut && adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut && adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut ){
					good_paddle[nmultiplicity-1]=ipaddle;
	      			   }
				  }
	    			}
           		}
		}
	   }

           if(nmultiplicity>0&&nmultiplicity<=multiplicity_cut){
             for(Int_t icount=0;icount<nmultiplicity;icount++){
		     hoccupancy->Fill(good_paddle[icount]);
	     	     heatmap->Fill(good_paddle[icount],nmultiplicity);
	     }
           }
           hmultiplicity->Fill(nmultiplicity);
        }

        cOCCUPANCY->Clear();
        cOCCUPANCY->Divide(1,2) ;

        title.Form("run_%d_OCCUPANCY_tdc_min_%d_max_%d.png",
                   run,tdc_min,tdc_min+tdc_width);
        cOCCUPANCY->Print(title);
        cOCCUPANCY->cd(1);
	//gPad->SetLogy();
        hoccupancy->Draw();
        hoccupancy->GetXaxis()->SetNdivisions(NUMPMT,NUMPMT,0,0);
        hoccupancy->SetLineColor(kBlue);
        gPad->SetGridx();
        TVirtualPad *c1_1 = cOCCUPANCY->cd(2);
	c1_1->Divide(2,1);
	c1_1->cd(1);
        hmultiplicity->Draw();
        hmultiplicity->SetLineColor(kBlue);
	c1_1->cd(2);
        heatmap->Draw("COLZ");
        heatmap->GetXaxis()->SetNdivisions(NUMPMT,NUMPMT,0,0);
        gPad->SetGridx();
        return cOCCUPANCY;
}

TCanvas *plot_occupancy_single(Int_t event_id = 1, Int_t adc_cut=40, Int_t multiplicity_cut=12, Int_t tdc_min=750, Int_t tdc_width=200){

        TString cut, draw, draw1, title;
        title.Form("run_%d_Occupancy",run);

        Int_t nbin=196;
        Int_t min=1, max=197;
        Int_t nbinm=11;
        Int_t minm=-1, maxm=10;
        TH1D *hoccupancy = new TH1D("hoccupancy","hoccupancy",nbin,min,max);
        TH1D *hmultiplicity = new TH1D("hmultiplicity","hmultiplicity",nbinm,minm,maxm);

        TCanvas *cOCCUPANCYE= new TCanvas("cOCCUPANCYE",title,xcanvas,ycanvas);

        Int_t nentries=n_events_to_analyze;
        
	TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.9);
        MyStyle->SetStatW(0.2);

        //for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(event_id);
	  
          Int_t nmultiplicity=0;
          Int_t good_paddle[500];
          for(Int_t icount=0;icount<500;icount++){good_paddle[icount]=-1;}

          for(Int_t pmt=0; pmt<NUMPMT;pmt++){
		Int_t ipaddle = (pmt+1)*NUMPADDLE+1;
		for(Int_t pixel=0; pixel<NUMPIXEL;pixel++){
			Int_t index = pmt*NUMPIXEL+pixel;
			if (pixel!=pixel1[pmt]-1&&pixel!=pixel2[pmt]-1){
				ipaddle--;
	    			if (tdcl[index]>tdc_min&&tdcl[index]<tdc_min+tdc_width){
	      			  nmultiplicity++;
	      			  if (adc_c[index]>adc_cut){
					good_paddle[nmultiplicity-1]=ipaddle;
	      			  }
	    			}
           		}
		}
	   }

           if(nmultiplicity>0&&nmultiplicity<=multiplicity_cut){
             for(Int_t icount=0;icount<nmultiplicity;icount++){hoccupancy->Fill(good_paddle[icount]);}
           }
           hmultiplicity->Fill(nmultiplicity);
        //}

        cOCCUPANCYE->Clear();
        cOCCUPANCYE->Divide(1,2) ;

        title.Form("run_%d_OCCUPANCYE_tdc_min_%d_max_%d.png",
                   run,tdc_min,tdc_min+tdc_width);
        cOCCUPANCYE->Print(title);
        cOCCUPANCYE->cd(1);
        hoccupancy->Draw();
        hoccupancy->GetXaxis()->SetNdivisions(NUMPMT,NUMPMT,0,0);
        hoccupancy->SetLineColor(kBlue);
        gPad->SetGridx();
        cOCCUPANCYE->cd(2);
        hmultiplicity->Draw();
        hmultiplicity->SetLineColor(kBlue);
        return cOCCUPANCYE;
}

TCanvas *plot_ratio(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=300){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADCRATIO",run);
        TCanvas *cADCRATIO= new TCanvas("cADCRATIO",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);
        TH1D *htmpc[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=100;
        Int_t min=0, max=500;
        for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpc%d", icounter);
                htmpc[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
                htmpc[icounter - 1]->SetLineColor(kRed);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title);
                htmpc[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(adc_c[ipaddle]);
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpb[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}

	cADCRATIO->Clear();
        cADCRATIO->Divide(4,4);

        //plot histos
        Int_t icount=0;
        TF1 *myfit = new TF1("myfit","1.0-0.5*ROOT::Math::erfc((x-[0])/[1])",0,1);
        myfit->SetParName(0,"Mean");
        myfit->SetParName(1,"Width");

        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cADCRATIO->cd( icount + 1 );
            //gPad->SetLogy();

            //cADC->Update();
                
            Int_t entries = htmpa[i]->GetEntries();
            float mean = htmpa[i]->GetMean(1);
            float RMS = htmpa[i]->GetRMS(1);
            //cout << entries <<" "<< mean <<" "<< RMS <<endl;
                
            htmpb[i]->SetStats(0);
            //current->Modified();

            myfit->SetParameter(0,40.0);
            myfit->SetParameter(1,10.0);
            htmpc[i] = (TH1D*)htmpb[i]->Clone();
            htmpc[i]->Divide(htmpa[i]);
            htmpc[i]->Fit("myfit");

	    // write the mean and sigma of the fit 
	    // on the canvas itself

	    float fitmean = myfit->GetParameter(0);
	    float fitsig = myfit->GetParameter(1);

	    // create strings to write on plot
	    TString lmean = Form("mean = %1.1f",fitmean);
	    TString lsig = Form("sigma = %1.1f",fitsig);
	    
	    // create the latex object to write
	    TLatex *tl;
	    tl->SetTextSize(0.07);
	    tl->SetTextAngle(0);
	    tl->SetTextColor(1);

	    // NOTE: first two numbers are the x and y VALUES
	    // on the individual graph/histogram
	    tl->DrawLatex(30,0.85,lmean);
	    tl->DrawLatex(30,0.75,lsig);
	    
            icount++;
          
	  }
        }

        title.Form("run_%d_ADCRATIO_pmt_%d_tdc_min_%d_max_%d.png",
                   run,pmt,tdc_min,tdc_min+tdc_width);
        cADCRATIO->Print(title);
        cADCRATIO->cd(0);
        return cADCRATIO;
}

TCanvas *plot_mean_ratio(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=200){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADCMEANRATIO",run);
        TCanvas *cADCMEANRATIO= new TCanvas("cADCMEANRATIO",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);
        TH1D *htmpc[NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=25;
        Int_t min=0, max=100;
        for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb%d", icounter);
                htmpb[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpc%d", icounter);
                htmpc[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[icounter - 1]->SetLineColor(kBlue);
                htmpb[icounter - 1]->SetLineColor(kRed);
                htmpc[icounter - 1]->SetLineColor(kRed);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
                htmpb[icounter - 1]->SetTitle(title);
                htmpc[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
	  T->GetEntry(id);
	 
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[index-1]->Fill(adc_c[ipaddle]);
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpb[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}

	cADCMEANRATIO->Clear();

        //plot histos
        Int_t icount=0;
	cADCMEANRATIO->cd();

        Double_t mean[NUMPADDLE]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        Double_t sigma[NUMPADDLE]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t paddle[NUMPADDLE]={1,2,3,4,5,6,7,8,9,10,11,12,13,14};
        Double_t epaddle[NUMPADDLE]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        TF1 *myfit = new TF1("myfit","1.0-0.5*ROOT::Math::erfc((x-[0])/[1])",0,1);
        myfit->SetParName(0,"Mean");
        myfit->SetParName(1,"Width");

        for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            htmpb[i]->SetStats(0);

            myfit->SetParameter(0,40.0);
            myfit->SetParameter(1,10.0);
            htmpc[i] = (TH1D*)htmpb[i]->Clone();
            htmpc[i]->Divide(htmpa[i]);
            htmpc[i]->Fit("myfit");

	    Int_t numentries = htmpc[i]->GetEntries();
	    mean[icount]=myfit->GetParameter(0);
	    //sigma[icount]=myfit->GetParameter(1)/sqrt(numentries);
	    sigma[icount]=myfit->GetParameter(1);

            icount++;
          
	  }
        }

        TGraphErrors *gr = new TGraphErrors(NUMPADDLE,paddle,mean,epaddle,sigma);
        gr->SetMarkerStyle(21);
        gr->GetXaxis()->SetTitle("Paddle Number");
        gr->GetYaxis()->SetTitle("50% Threshold (Good TDC)");
	gr->GetHistogram()->SetMaximum(100);
	gr->GetHistogram()->SetMinimum(0);

        gr->Draw("AP");

        cADCMEANRATIO->Update();

        title.Form("run_%d_ADCMEANRATIO_pmt_%d_tdc_min_%d_max_%d.png",
                   run,pmt,tdc_min,tdc_min+tdc_width);
        cADCMEANRATIO->Print(title);
        cADCMEANRATIO->cd();
        return cADCMEANRATIO;

}

TCanvas *plot_mean_adc(Int_t pmt=1, Int_t tdc_min=750, Int_t tdc_width=200){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADCMean",run);
        TCanvas *cADCMean= new TCanvas("cADCMean",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=600;
        Int_t min=-100, max=1000;
        for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa%d", icounter);
                htmpa[icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[icounter - 1]->SetTitle(title);
        }

        Int_t nentries=n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    if(tdcl[ipaddle]>tdc_min&&tdcl[ipaddle]<tdc_min+tdc_width){
		htmpa[index-1]->Fill(adc_c[ipaddle]);
	    }
	  }
	}

	cADCMean->Clear();

        //plot histos
        Int_t icount=0;
        cADCMean->cd();

	Double_t mean[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t sigma[NUMPIXEL]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t pixel[NUMPIXEL]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	Double_t epixel[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t mmean[2]={0,0};
	Double_t msigma[2]={1,1};
	Double_t mpixel[2]={0,0};
	Double_t mepixel[2]={0,0};
	mpixel[0]=pixel1[pmt-1];
	mpixel[1]=pixel2[pmt-1];

	Double_t yline[2]={100,100};
	Double_t xline[2]={0,17};
        
	for (Int_t i=0; i<NUMPIXEL; i++){

          if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            Int_t entries = htmpa[i]->GetEntries();
            mean[i] = htmpa[i]->GetMean(1);
            sigma[i] = htmpa[i]->GetRMS(1)/sqrt(entries);
            cout << i+1 << " "<< entries <<" "<< mean[i] <<" "<< sigma[i] <<endl;
                
            icount++;
          
	  }
        }

	TGraphErrors *gr = new TGraphErrors(NUMPIXEL,pixel,mean,epixel,sigma);
	gr->SetMarkerStyle(21);
	gr->GetXaxis()->SetTitle("Pixel Number");
	gr->GetYaxis()->SetTitle("Mean ADC (Good TDC)");
	TGraphErrors *gr2 = new TGraphErrors(2,mpixel,mmean,mepixel,msigma);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(3);
	TGraph *gr3 = new TGraph(2,xline,yline);
	gr3->SetLineColor(3);
	gr3->SetLineWidth(4);

	gr->Draw("AP");
	gr2->Draw("P");
	gr3->Draw("L");
            
        cADCMean->Update();

        title.Form("run_%d_ADCMean_pmt_%d_tdc_min_%d_max_%d.png",
                   run,pmt,tdc_min,tdc_min+tdc_width);
        cADCMean->Print(title);
        cADCMean->cd(0);
        return cADCMean;
}

TCanvas *calibrate_adc(Int_t ped = 1, Int_t tdc_min=750, Int_t tdc_width=200){

        TString cut, draw, draw1, title;
        title.Form("run_%d_ADCCalibrate",run);
        TCanvas *cADCCalib= new TCanvas("cADCCalib",title,xcanvas,ycanvas);

        TH1D *htmpa[NUMPMT][NUMPIXEL];//=new TH1D("htmp","htmp",nbin,min,max);
        TH1D *htmpb[NUMPMT][NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);
        TH1D *htmpc[NUMPMT][NUMPIXEL];//=new TH1D("htmp1","htmp1",nbin,min,max);

        TString tmpentry;
        MyStyle->SetStatX(0.9);
        MyStyle->SetStatY(0.6);
        MyStyle->SetStatW(0.4);

        Int_t nbin=600;
        Int_t min=-100, max=1000;
	for (Int_t pmt=1; pmt<=NUMPMT; pmt++){
         for(Int_t icounter=1;icounter<=NUMPIXEL;icounter++){
                tmpentry.Form("htmpa_%d_%d", pmt,icounter);
                htmpa[pmt-1][icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpb_%d_%d", pmt,icounter);
                htmpb[pmt-1][icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                tmpentry.Form("htmpc_%d_%d", pmt,icounter);
                htmpc[pmt-1][icounter - 1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
                htmpa[pmt-1][icounter - 1]->SetLineColor(kBlue);
                htmpb[pmt-1][icounter - 1]->SetLineColor(kRed);
                htmpc[pmt-1][icounter - 1]->SetLineColor(kGreen);
                title.Form("Run %d ADC pmt %d, paddle %d: %d < tdc < %d",run,pmt,icounter,tdc_min,tdc_min+tdc_width);
                htmpa[pmt-1][icounter - 1]->SetTitle(title);
                htmpb[pmt-1][icounter - 1]->SetTitle(title);
                htmpc[pmt-1][icounter - 1]->SetTitle(title);
         }
	}

        Int_t nentries=n_events_to_analyze;

        for (Int_t id=1;id<=nentries;id++){
          T->GetEntry(id);
	  for (Int_t pmt=1; pmt<=NUMPMT; pmt++){
	   for (Int_t index=1; index<=NUMPIXEL; index++){
	    Int_t ipaddle = (pmt-1)*NUMPIXEL+index-1;
	    htmpa[pmt-1][index-1]->Fill(adc[ipaddle]);
	    if(tdcl[ipaddle]<tdc_min||tdcl[ipaddle]>tdc_min+tdc_width){
		htmpb[pmt-1][index-1]->Fill(adc[ipaddle]);
	    }else{
		htmpc[pmt-1][index-1]->Fill(adc[ipaddle]);
	    }
	   }
	  }
	}

	cADCCalib->Clear();
        cADCCalib->Divide(4,4);

        //plot histos
	for (Int_t pmt=1; pmt<=NUMPMT; pmt++){
         Int_t icount=0;
         for (Int_t i=0; i<NUMPIXEL; i++){

          //if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) {

            //cout<<"into loop 2, i = " << i << endl;

            cADCCalib->cd( icount + 1 );
            gPad->SetLogy();

            cADCCalib->Update();
                
            Int_t entries = htmpb[pmt-1][i]->GetEntries();
            float mean = htmpb[pmt-1][i]->GetMean(1);
            float RMS = htmpb[pmt-1][i]->GetRMS(1);
	    Int_t imean = mean;
            float mean2 = htmpc[pmt-1][i]->GetMean(1);
   	    float gain;
	    if (mean2-mean > 2.0){
	    	gain = 100.0/(mean2-mean);
	    }else{
		gain = 1.0;
	    }
            //cout << entries <<" "<< mean <<" "<< RMS <<endl;
            if (ped == 1){
		if(pmt==1&&i==0) cout << "C.cdetm1r.ped = ";
	    	cout << imean << " ";
	    }else{
		if(pmt==1&&i==0) cout << "C.cdetm1r.gain = ";
	    	cout << gain << " ";
	    }
		
                
            htmpa[pmt-1][i]->SetStats(0);
            //current->Modified();
            htmpa[pmt-1][i]->Draw();
            htmpb[pmt-1][i]->Draw("same");
            htmpc[pmt-1][i]->Draw("same");

            icount++;
          
	  //}
         }
	cout << endl;
	}

        title.Form("run_%d_ADCCalib_tdc_min_%d_max_%d.png",
                   run,tdc_min,tdc_min+tdc_width);
        cADCCalib->Print(title);
        cADCCalib->cd(0);
        return cADCCalib;
}
