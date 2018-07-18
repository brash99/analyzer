// 25th May 2018
//
// Script to make a comparison plot of the efficiency ratio
// values as a function of high voltage settings for PMT 7 
// between runs 1317 - 1320, 1324 with resistors installed.
// The threshold numbers for each setting have to be read in.


#ifndef __CINT__
#include <iostream>
#include <fstream>
#endif
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TFormula.h"
#include "TTree.h"
#include "TPad.h"
#include "TArrow.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "TPostScript.h"

using namespace std;

static const Int_t NUMPMT = 14;
static const Int_t NUMPIXEL = 16;
static const Int_t NUMPADDLE = 14;
static const Int_t NUMPADDLES = NUMPMT*NUMPIXEL;
static const Int_t NUMRUNS = 5;
static const Int_t NUMBINS = NUMRUNS*NUMPADDLE;

double_t adc[NUMPADDLES];
double_t adc_c[NUMPADDLES];
double_t tdct[NUMPADDLES];
double_t tdcl[NUMPADDLES];

//Latest map set for M1-L
//Int_t pixel1[NUMPMT]={1, 2, 13, 4, 6, 12, 4,1, 1,5,13,1, 1,8};
//Int_t pixel2[NUMPMT]={14, 16,16,6,13,13,16,3,6,9,15,4,16,12};
//Map set for M1-R
Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13};
Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};

void effplot_ejb_bool(){

  // Set some of the styles first
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleOffset(1.5,"a");
  gStyle->SetTitleOffset(0.9,"x");
  // gStyle->SetTitleOffset(1.7,"y");
  gStyle->SetPadColor(10);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(1);
  gStyle->SetStatColor(10); 
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.16);
  gStyle->UseCurrentStyle();
  gStyle->SetMarkerStyle(22);
  gStyle->SetMarkerColor(kBlack);
  gStyle->SetLineColor(kBlack);
  gStyle->SetMarkerSize(1.2);

  Int_t paddleindex[NUMPADDLES];

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

  // ------------------------------------------------

  // Read in the data from a text file.
  // First line are run numbers corresponding to runs at 
  // HV settings of 700, 725, 750, 775 and 800 V respectively

  // define arrays for threshold values from each run.
 
  Int_t thr[NUMRUNS][NUMPIXEL];

  //Int_t* thr1 = new Int_t[NUMPIXEL];     // 700 V
  //Int_t* thr2 = new Int_t[NUMPIXEL];     // 725 V
  //Int_t* thr3 = new Int_t[NUMPIXEL];     // 750 V
  //Int_t* thr4 = new Int_t[NUMPIXEL];     // 775 V
  //Int_t* thr5 = new Int_t[NUMPIXEL];     // 800 V

  Int_t* pixel = new Int_t[NUMPIXEL];   // pixel number in PMT
  Int_t pmt = 0;   // tag for PMT number
  Int_t rn[NUMRUNS];  // run numbers for root files
  Int_t hv[NUMRUNS];  // HV settings
  Int_t mp1 = 0;   // missing pixel numbers
  Int_t mp2 = 0; 
  Double_t* bf = new Double_t[NUMBINS];   // counts before threshold cut
  Double_t* af = new Double_t[NUMBINS];   // counts after threshold cut
  Double_t* eff = new Double_t[NUMBINS];  // efficiency numbers
  string tx,px;

  // array of TDC width cuts for different HV
  Int_t tdcw[NUMRUNS][NUMPIXEL] = 
  {24, 32, 32, 0, 32, 26, 26, 32, 0, 24, 32, 32, 26, 26, 30, 0,
  Int_t tdcw[NUMRUNS][NUMPIXEL] = 
  {24, 32, 32, 0, 32, 26, 26, 32, 0, 24, 32, 32, 26, 26, 30, 0,
   30, 33, 33, 0, 34, 30, 30, 34, 0, 30, 34, 34, 30, 30, 33, 0,
   32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0,
   32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0,
   32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0};

  // read in the input file of threshold values
  ifstream tt; 
  tt.open("testtext.txt",ios::in);  
  if(tt.is_open()){
    tt >> tx >> rn[0] >> rn[1] >> rn[2] >> rn[3] >> rn[4];
    tt >> px >> pmt >> hv[0] >> hv[1] >> hv[2] >> hv[3] >> hv[4];
    cout << px << "\t" << pmt << endl;
    cout << tx << "\t" << rn[0] << "\t" << rn[1] << "\t" << rn[2] << "\t" << rn[3] << "\t" << rn[4] << endl;
    cout << "HV:" << "\t" << hv[0] << "\t" << hv[1] << "\t" << hv[2] << "\t" << hv[3] << "\t" << hv[4] << endl;
    for(Int_t i=0; i<NUMPIXEL; i++){
      tt >> pixel[i]>> thr[0][i] >> thr[1][i] >> thr[2][i] >> thr[3][i] >> thr[4][i];
      cout<<pixel[i]<<"\t"<<thr[0][i]<<"\t"<<thr[1][i]<<"\t"<<thr[2][i]<<"\t"<<thr[3][i]<<"\t"<<thr[4][i]<<endl;
    }
  }
  tt.close();
 
  for(Int_t i=0; i<NUMPIXEL; i++){
    if(thr[4][i] == 0 && mp1 == 0){
      mp1 = i+1;
    }
    if(thr[4][i] == 0 && mp1>0){
      mp2 = i+1;
    }
  }

  cout << "missing pixels are: " << mp1 << "\t" << mp2 << endl;

  // create the string for the root filenames
  
  TString fname[NUMRUNS];
  TFile *f[NUMRUNS];
  TTree *t[NUMRUNS];
  Int_t numentries[NUMRUNS];

  for (Int_t index = 0; index<NUMRUNS; index++){
  	fname[index] = Form("../rootfiles/scint_%d.root",rn[index]);
  	f[index] = TFile::Open(fname[index],"OPEN");
	f[index]->cd();
	t[index] = (TTree*)f[index]->Get("T");
	t[index]->SetBranchAddress("C.cdetm1r.adc_c",&adc_c);
	t[index]->SetBranchAddress("C.cdetm1r.adc",&adc_c);
	t[index]->SetBranchAddress("C.cdetm1r.tdcl",&tdcl);
	t[index]->SetBranchAddress("C.cdetm1r.tdct",&tdct);
  	numentries[index] = t[index]->GetEntries();
  }

  // Figure out the root indices for all the pixels

  Int_t pc = ((pmt-1)*NUMPIXEL);

  // Create all fourteen canvases and all the histograms now

  TCanvas *sc[NUMPIXEL-1];

  TH1D *hh[NUMBINS]; // = new TH1D[NUMBINS];   // adc spectrum with TDC & ADC neighbors cut
  TH1D *ha[NUMBINS]; // = new TH1D[NUMBINS];   // adc with TDC crosstalk cut added to previous cuts
  TH1D *hb[NUMBINS]; // = new TH1D[NUMBINS];   // adc with TDC width cut included also
  TH1D *hc[NUMBINS]; // = new TH1D[NUMBINS];   // adc with threshold cut added to other cuts

  TString tmpentry;
  for(Int_t counter=0;counter<NUMBINS;counter++){
	  tmpentry.Form("hh%d",counter);
	  hh[counter] = new TH1D(tmpentry,tmpentry,1000,0,1000); 
	  tmpentry.Form("ha%d",counter);
	  ha[counter] = new TH1D(tmpentry,tmpentry,1000,0,1000); 
	  tmpentry.Form("hb%d",counter);
	  hb[counter] = new TH1D(tmpentry,tmpentry,1000,0,1000); 
	  tmpentry.Form("hc%d",counter);
	  hc[counter] = new TH1D(tmpentry,tmpentry,1000,0,1000); 
  }

  Bool_t temph1,temph2,temph3;
  Bool_t* hcuts = new Bool_t[NUMPADDLE];
  Bool_t* acuts = new Bool_t[NUMPADDLE];
  Bool_t* bcuts = new Bool_t[NUMBINS];
  Bool_t* ccuts = new Bool_t[NUMBINS];

  Int_t ipaddle = (pmt)*NUMPADDLE+1;
  Int_t hcutindex = 0;
	    
  for (Int_t pixel=0; pixel < NUMPIXEL; pixel++){
    Int_t index = (pmt-1)*NUMPIXEL+pixel;
    if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1){
      ipaddle--;
      if (ipaddle < 2){
	temph2 = Form("C.cdetm1r.adc_c[%d]<20",paddleindex[ipaddle+1]);
	temph3 = Form("C.cdetm1r.tdcl[%d]>500",paddleindex[ipaddle]);
	hcuts[hcutindex] = temph2 && temph3;
	//cout << hcuts[hcutindex].GetTitle() << endl;
	hcutindex++;
      }
      else if (ipaddle > NUMPMT*NUMPADDLE){
	temph1 = Form("C.cdetm1r.adc_c[%d]<20",paddleindex[ipaddle-1]);
	temph3 = Form("C.cdetm1r.tdcl[%d]>500",paddleindex[ipaddle]);
	hcuts[hcutindex] = temph1 && temph3;
	//cout << hcuts[hcutindex].GetTitle() << endl;
	hcutindex++;
      }
      else{
	temph1 = Form("C.cdetm1r.adc_c[%d]<20",paddleindex[ipaddle-1]);
	temph2 = Form("C.cdetm1r.adc_c[%d]<20",paddleindex[ipaddle+1]);
	temph3 = Form("C.cdetm1r.tdcl[%d]>500",paddleindex[ipaddle]);
	hcuts[hcutindex] = temph1 && temph2 && temph3;
	//cout << hcuts[hcutindex].GetTitle() << endl;
	hcutindex++;
      }
    }
  }

  Bool_t tempa;
  Int_t ipad = (pmt)*NUMPADDLE+1;
  Int_t acutindex = 0;
 
  for(Int_t pix=0;pix<NUMPIXEL;pix++){
    Int_t ind = (pmt)*NUMPADDLE;
    if (pix!=pixel1[pmt-1]-1 && pix!=pixel2[pmt-1]-1){
      ipad--;
      if (ipad < 2){
	for(Int_t s=0;s<NUMPADDLE;s++){
	  if(ind-s != ipad){
	    acuts[acutindex] += Form("C.cdetm1r.tdcl[%d]<200",paddleindex[ind-s]);
	  }
	}
	tempa = hcuts[acutindex];
	acuts[acutindex] = acuts[acutindex] && tempa;
	//cout << acuts[acutindex].GetTitle() << endl;
	acutindex++;
      }
      else if(ipad > NUMPMT*NUMPADDLE){
	for(Int_t s=0;s<NUMPADDLE;s++){
	  if(ind-s != ipad){
	    acuts[acutindex] += Form("C.cdetm1r.tdcl[%d]<200",paddleindex[ind-s]);
	  }
	}
	tempa = hcuts[acutindex];
	acuts[acutindex] = acuts[acutindex] && tempa;
	//cout << acuts[acutindex].GetTitle() << endl;
	acutindex++;
      }
      else{
	for(Int_t s=0;s<NUMPADDLE;s++){
	  if(ind-s != ipad){
	    acuts[acutindex] += Form("C.cdetm1r.tdcl[%d]<200",paddleindex[ind-s]);
	  }
	}
	tempa = hcuts[acutindex];
	acuts[acutindex] = acuts[acutindex] && tempa;
	//cout << acuts[acutindex].GetTitle() << endl;
	acutindex++;
      }
    }
    if(acutindex == NUMPADDLE){
      break;
    }
  }
  
  Bool_t tempb;
  Int_t hacutindex;
  Int_t bcutindex;
  Int_t tdcwindex;

  for (Int_t hv=1; hv < 6; hv++){
    hacutindex = 0;
    tdcwindex = 0;
    bcutindex = 0;
    for(Int_t pxl=0;pxl<NUMPIXEL;pxl++){
      Int_t idx = (pmt)*NUMPADDLE+pxl-2;
      if (pxl!=pixel1[pmt-1]-1 && pxl!=pixel2[pmt-1]-1){
	tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw[hv-1][tdcwindex]);
	bcuts[bcutindex+(hv-1)*NUMPADDLE] = acuts[hacutindex] && tempb;
	//cout << "HV = " << hv << " Cut # " << bcutindex+(hv-1)*NUMPADDLE << bcuts[bcutindex+(hv-1)*NUMPADDLE].GetTitle() << endl;
	bcutindex++;
	hacutindex++;
      }
      tdcwindex++;
    }
  }
 
  Bool_t tempc;
  Int_t habcutindex;
  Int_t ccutindex;
  Int_t thrindex;

  for(Int_t HV=1; HV < 6; HV++){
    habcutindex = 0;
    thrindex = 0;
    ccutindex = 0;
    for(Int_t pxel=0; pxel < NUMPIXEL; pxel++){
      Int_t inx = (pmt)*NUMPADDLE+pxel-2;
      if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr[HV-1][thrindex]);
	  ccuts[ccutindex+(HV-1)*NUMPADDLE] = bcuts[ccutindex+(HV-1)*NUMPADDLE] && tempc;
	  //cout << "HV = " << HV << " Cut # " << ccutindex+(HV-1)*NUMPADDLE << ccuts[ccutindex+(HV-1)*NUMPADDLE].GetTitle() << endl;
	  ccutindex++;
	  habcutindex++;
      }
      thrindex++;
    }
  }

  cout << "Starting tree looping now!!!" << endl;


  Int_t canindexd=0;
  for (Int_t canvasd=0;canvasd<NUMPIXEL;canvasd++){
   if (canvasd!=pixel1[pmt-1]-1 && canvasd!=pixel2[pmt-1]-1){
    TString cname = Form("sc%d",canindexd+1);
    TString ctitle = Form("Pixel %d",canvasd+1);
    sc[canindexd] = new TCanvas(cname,ctitle,1300,400);
    sc[canindexd]->Divide(5,1);
    sc[canindexd]->Draw();
   }
   canindexd++;
  }

  for (Int_t counter=0;counter<NUMBINS;counter++){
	  af[counter]=0;
	  bf[counter]=0;
	  eff[counter]=0;
  }

  TString temptitle;
  for(Int_t hvolt=1; hvolt < 6; hvolt++){
    if (id%1000==0) cout << "Event = " << id << endl;
    Int_t canindex = 0;
    for(Int_t canvas=0; canvas < NUMPIXEL; canvas++){
     if(canvas!=pixel1[pmt-1]-1 && canvas!=pixel2[pmt-1]-1){
        paddle--;
        temptitle.Form("C.cdetm1r.adc_c[%d]",paddleindex[paddle]);

        //sc[canindex]->cd(hvolt);
	
	//if(hcuts[canindex]) hh[canindex+(hvolt-1)*NUMPADDLE]->Fill(temptitle);
	if(hcuts[canindex]) hh[canindex+(hvolt-1)*NUMPADDLE]->Fill(C.cdetm1r.adc_c[paddleindex[paddle]]);
	
	if(acuts[canindex]) ha[canindex+(hvolt-1)*NUMPADDLE]->Fill(temptitle);
	
	if(bcuts[canindex+(hvolt-1)*NUMPADDLE]) {
		hb[canindex+(hvolt-1)*NUMPADDLE]->Fill(temptitle);
		bf[canindex+(hvolt-1)*NUMPADDLE]++;
	}
	
	if(ccuts[canindex+(hvolt-1)*NUMPADDLE]) {
		hc[canindex+(hvolt-1)*NUMPADDLE]->Fill(temptitle);
		af[canindex+(hvolt-1)*NUMPADDLE]++;
	}
	
	canindex++;
       }
      }
    }
  }

  cout << "Finished tree looping ... " << endl;
  
  for (Int_t counter=0;counter<NUMBINS;counter++){
	  eff[counter]=af[counter]/bf[counter];
  }

  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(3);

  ofstream outf;  // output file
  TString outftitle;
  outftitle.Form("PMT_%d_Efficiency.txt",pmt);
  outf.open(outftitle,ios::out);

  outf.setf(ios::fixed);
  outf.setf(ios::showpoint);
  outf.precision(3);

  outf << "PMT " << pmt << " Efficiency" << endl;
  outf << tx << "\t" << rn[0] << "\t" << rn[1] << "\t" << rn[2] << "\t" << rn[3] << "\t" << rn[4] << endl;
  outf << "HV:" << "\t" << hv[0] << "\t" << hv[1] << "\t" << hv[2] << "\t" << hv[3] << "\t" << hv[4] << endl;
  cout << "PMT " << pmt << " Efficiency" << endl;
  cout << tx << "\t" << rn[0] << "\t" << rn[1] << "\t" << rn[2] << "\t" << rn[3] << "\t" << rn[4] << endl;
  cout << "HV:" << "\t" << hv[0] << "\t" << hv[1] << "\t" << hv[2] << "\t" << hv[3] << "\t" << hv[4] << endl;

  Int_t effindex = 0;

  for(Int_t pid=0; pid < NUMPIXEL; pid++){
    if(pid!=pixel1[pmt-1]-1 && pid!=pixel2[pmt-1]-1){
      outf << pid+1 << "\t" << eff[effindex] << "\t" << eff[effindex+NUMPADDLE] << "\t" << eff[effindex+2*NUMPADDLE] << "\t" <<
	eff[effindex+3*NUMPADDLE] << "\t" << eff[effindex+4*NUMPADDLE] << endl;
      cout << pid+1 << "\t" << eff[effindex] << "\t" << eff[effindex+NUMPADDLE] << "\t" << eff[effindex+2*NUMPADDLE] << "\t" <<
	eff[effindex+3*NUMPADDLE] << "\t" << eff[effindex+4*NUMPADDLE] << endl;
      effindex++;
    }
  }

  outf.close();

  TString gname;
  TCanvas *ec[NUMPADDLE];
  TGraph *gc[NUMPADDLE];
  TF1 *fun1 = new TF1("fun1","pol3",700,800);
  fun1->SetLineColor(1);
  fun1->SetLineWidth(2);
  fun1->SetLineStyle(7);

  Int_t ecindex = 0;
  for(Int_t pind=0; pind < NUMPIXEL; pind++){
    if(pind!=pixel1[pmt-1]-1 && pind!=pixel2[pmt-1]-1){
      TString ecname = Form("ec%d",effindex+1);
      TString ectitle = Form("Pixel %d Eff",pind+1);
      ec[ecindex] = new TCanvas(ecname,ectitle,800,800);
      ec[ecindex]->cd();
      if(ecindex==NUMPADDLE){
	break;
      }
      else {
	//gc[ecindex]->SetPoint(0,700,eff[ecindex]);
	//gc[ecindex]->SetPoint(1,725,eff[ecindex+NUMPADDLE]);
	//gc[ecindex]->SetPoint(2,750,eff[ecindex+2*NUMPADDLE]);
	//gc[ecindex]->SetPoint(3,775,eff[ecindex+3*NUMPADDLE]);
	//gc[ecindex]->SetPoint(4,800,eff[ecindex+4*NUMPADDLE]);
	//gc[ecindex]->Draw();
	//gc[ecindex]->Fit("fun1","B");
      }
      ecindex++;
    }
  }
  
  outf.close();

}

