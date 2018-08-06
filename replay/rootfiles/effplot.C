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

//Latest map set for M1-L
//Int_t pixel1[NUMPMT]={1, 2, 13, 4, 6, 12, 4,1, 1,5,13,1, 1,8};
//Int_t pixel2[NUMPMT]={14, 16,16,6,13,13,16,3,6,9,15,4,16,12};
//Map set for M1-R
Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13};
Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};

void effplot(){

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

	cout << "tdc daq index = " << index << "\t paddleindex = " << ipaddle << endl;
      }
    }
  }

  // ------------------------------------------------

  // Read in the data from a text file.
  // First line are run numbers corresponding to runs at 
  // HV settings of 700, 725, 750, 775 and 800 V respectively

  // define arrays for threshold values from each run.
  Int_t* thr1 = new Int_t[16];     // 700 V
  Int_t* thr2 = new Int_t[16];     // 725 V
  Int_t* thr3 = new Int_t[16];     // 750 V
  Int_t* thr4 = new Int_t[16];     // 775 V
  Int_t* thr5 = new Int_t[16];     // 800 V

  Int_t* pixel = new Int_t[16];   // pixel number in PMT
  Int_t pmt = 0;   // tag for PMT number
  Int_t rn1, rn2, rn3, rn4, rn5;  // run numbers for root files
  Int_t hv1, hv2, hv3, hv4, hv5;  // HV settings
  Int_t mp1 = 0;   // missing pixel numbers
  Int_t mp2 = 0; 
  Double_t* bf = new Double_t[70];   // counts before threshold cut
  Double_t* af = new Double_t[70];   // counts after threshold cut
  Double_t* eff = new Double_t[70];  // efficiency numbers
  string tx,px;

  // array of TDC width cuts for different HV
  Int_t tdcw1[16] = {24, 32, 32, 0, 32, 26, 26, 32, 0, 24, 32, 32, 26, 26, 30, 0};
  Int_t tdcw2[16] = {30, 33, 33, 0, 34, 30, 30, 34, 0, 30, 34, 34, 30, 30, 33, 0};
  Int_t tdcw3[16] = {32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0};
  Int_t tdcw4[16] = {32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0};
  Int_t tdcw5[16] = {32, 34, 34, 0, 34, 32, 32, 35, 0, 32, 34, 34, 32, 32, 34, 0};

  // read in the input file of threshold values
  ifstream tt; 
  tt.open("M1R_threshold_pmt3.dat",ios::in);  
  if(tt.is_open()){
    tt >> tx >> rn1 >> rn2 >> rn3 >> rn4 >> rn5;
    tt >> px >> pmt >> hv1 >> hv2 >> hv3 >> hv4 >> hv5;
    cout << px << "\t" << pmt << endl;
    cout << tx << "\t" << rn1 << "\t" << rn2 << "\t" << rn3 << "\t" << rn4 << "\t" << rn5 << endl;
    cout << "HV:" << "\t" << hv1 << "\t" << hv2 << "\t" << hv3 << "\t" << hv4 << "\t" << hv5 << endl;
    for(Int_t i=0; i<16; i++){
      tt >> pixel[i]>> thr1[i] >> thr2[i] >> thr3[i] >> thr4[i] >> thr5[i];
      cout<<pixel[i]<<"\t"<<thr1[i]<<"\t"<<thr2[i]<<"\t"<<thr3[i]<<"\t"<<thr4[i]<<"\t"<<thr5[i]<<endl;
    }
  }
  tt.close();
 
  for(Int_t i=0; i<16; i++){
    if(thr5[i] == 0 && mp1 == 0){
      mp1 = i+1;
    }
    if(thr5[i] == 0 && mp1>0){
      mp2 = i+1;
    }
  }

  cout << "missing pixels are: " << mp1 << "\t" << mp2 << endl;

  // create the string for the root filenames
  
  TString fname1 = Form("../rootfiles/scint_%d.root",rn1);
  TString fname2 = Form("../rootfiles/scint_%d.root",rn2);
  TString fname3 = Form("../rootfiles/scint_%d.root",rn3);
  TString fname4 = Form("../rootfiles/scint_%d.root",rn4);
  TString fname5 = Form("../rootfiles/scint_%d.root",rn5);

  // Now open the root files

  TFile *f1 = TFile::Open(fname1,"OPEN");
  TFile *f2 = TFile::Open(fname2,"OPEN");
  TFile *f3 = TFile::Open(fname3,"OPEN");
  TFile *f4 = TFile::Open(fname4,"OPEN");
  TFile *f5 = TFile::Open(fname5,"OPEN");

  // Figure out the root indices for all the pixels

  Int_t pc = ((pmt-1)*16);

  // Get access to the tree in each root file

  f1->cd();
  TTree *t1 = (TTree*)f1->Get("T");

  f2->cd();
  TTree *t2 = (TTree*)f2->Get("T");

  f3->cd();
  TTree *t3 = (TTree*)f3->Get("T");

  f4->cd();
  TTree *t4 = (TTree*)f4->Get("T");

  f5->cd();
  TTree *t5 = (TTree*)f5->Get("T");

  // Create all fourteen canvases and all the histograms now

  TCanvas *sc[15];

  TH1D *hh = new TH1D[70];   // adc spectrum with TDC & ADC neighbors cut
  TH1D *ha = new TH1D[70];   // adc with TDC crosstalk cut added to previous cuts
  TH1D *hb = new TH1D[70];   // adc with TDC width cut included also
  TH1D *hc = new TH1D[70];   // adc with threshold cut added to other cuts

  Int_t numentries1 = t1->GetEntries();
  Int_t numentries2 = t2->GetEntries();
  Int_t numentries3 = t3->GetEntries();
  Int_t numentries4 = t4->GetEntries();
  Int_t numentries5 = t5->GetEntries();

  TCut temph1,temph2,temph3;
  TCut* hcuts = new TCut[14];
  TCut* acuts = new TCut[14];
  TCut* bcuts = new TCut[70];
  TCut* ccuts = new TCut[70];

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
	cout << hcuts[hcutindex].GetTitle() << endl;
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

  TCut tempa;
  Int_t ipad = (pmt)*NUMPADDLE+1;
  Int_t acutindex = 0;
 
  for(Int_t pix=0;pix<NUMPIXEL;pix++){
    Int_t ind = (pmt)*NUMPADDLE;
    if (pix!=pixel1[pmt-1]-1 && pix!=pixel2[pmt-1]-1){
      ipad--;
      if (ipad < 2){
	for(Int_t s=0;s<14;s++){
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
	for(Int_t s=0;s<14;s++){
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
	for(Int_t s=0;s<14;s++){
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
    if(acutindex == 14){
      break;
    }
  }
  
  TCut tempb;
  Int_t hacutindex;
  Int_t bcutindex = 0;
  Int_t tdcwindex = 0;

  for (Int_t hv=1; hv < 6; hv++){
    hacutindex = 0;
    tdcwindex = 0;
    for(Int_t pxl=0;pxl<NUMPIXEL;pxl++){
      Int_t idx = (pmt)*NUMPADDLE+pxl-2;
      if (pxl!=pixel1[pmt-1]-1 && pxl!=pixel2[pmt-1]-1){
	if(hv==1){
	  tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw1[tdcwindex]);
	  bcuts[bcutindex] = acuts[hacutindex] && tempb;
	  cout << "Cut # " << bcutindex << " " << bcuts[bcutindex].GetTitle() << endl;
	  bcutindex++;
	}
	else if(hv==2){
	  tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw2[tdcwindex]);
	  bcuts[bcutindex] = acuts[hacutindex] && tempb;
	  cout << "Cut # " << bcutindex << " " << bcuts[bcutindex].GetTitle() << endl;
	  bcutindex++;
	}
	else if(hv==3){
	  tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw3[tdcwindex]);
	  bcuts[bcutindex] = acuts[hacutindex] && tempb;
	  cout << "Cut # " << bcutindex << " " << bcuts[bcutindex].GetTitle() << endl;
	  bcutindex++;
	}
	else if(hv==4){
	  tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw4[tdcwindex]);
	  bcuts[bcutindex] = acuts[hacutindex] && tempb;
	  cout << "Cut # " << bcutindex << " " << bcuts[bcutindex].GetTitle() << endl;
	  bcutindex++;
	}
	else{
	  tempb = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d])>%d",idx,idx,tdcw5[tdcwindex]);
	  bcuts[bcutindex] = acuts[hacutindex] && tempb;
	  cout << "Cut # " << bcutindex << " " << bcuts[bcutindex].GetTitle() << endl;
	  bcutindex++;
	}
	hacutindex++;
      }
      tdcwindex++;
    }
  }
 
  TCut tempc;
  Int_t ccutindex = 0;
  Int_t habcutindex = 0;
  Int_t thrindex1 = 0;
  Int_t thrindex2 = 0;
  Int_t thrindex3 = 0;
  Int_t thrindex4 = 0;
  Int_t thrindex5 = 0;

  for(Int_t HV=1; HV < 6; HV++){
    for(Int_t pxel=0; pxel < NUMPIXEL; pxel++){
      Int_t inx = (pmt)*NUMPADDLE+pxel-2;
      if(HV==1){
	if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr1[thrindex1]);
	  ccuts[ccutindex] = bcuts[ccutindex] && tempc;
	  //cout << ccuts[ccutindex].GetTitle() << endl;
	  ccutindex++;
	}
	thrindex1++;
      }
      else if(HV==2){
	if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr2[thrindex2]);
	  ccuts[ccutindex] = bcuts[ccutindex] && tempc;
	  //cout << ccuts[ccutindex].GetTitle() << endl;
	  ccutindex++;
	}
	thrindex2++;
      }
      else if(HV==3){
	if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr3[thrindex3]);
	  ccuts[ccutindex] = bcuts[ccutindex] && tempc;
	  //cout << ccuts[ccutindex].GetTitle() << endl;
	  ccutindex++;
	}
	thrindex3++;
      }
      else if(HV==4){
	if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr4[thrindex4]);
	  ccuts[ccutindex] = bcuts[ccutindex] && tempc;
	  //cout << ccuts[ccutindex].GetTitle() << endl;
	  ccutindex++;
	}
	thrindex4++;
      }
      else{
	if(pxel!=pixel1[pmt-1]-1 && pxel!=pixel2[pmt-1]-1){
	  tempc = Form("C.cdetm1r.adc_c[%d]>%d", inx, thr5[thrindex5]);
	  ccuts[ccutindex] = bcuts[ccutindex] && tempc;
	  //cout << ccuts[ccutindex].GetTitle() << endl;
	  ccutindex++;
	}
	thrindex5++;
      }
    }
  }

  TString temptitle1, temptitle2, temptitle3, temptitle4;
  Int_t paddle = (pmt)*NUMPADDLE+1;
  Int_t canindex = 0;
  Int_t bcindex = 0;

  for(Int_t canvas=0; canvas < NUMPIXEL; canvas++){
    if(canvas!=pixel1[pmt-1]-1 && canvas!=pixel2[pmt-1]-1){
      paddle--;

      TString cname = Form("sc%d",canindex+1);
      TString ctitle = Form("Pixel %d",canvas+1);
      sc[canindex] = new TCanvas(cname,ctitle,1300,400);
      sc[canindex]->Divide(5,1);
      sc[canindex]->Draw();

      for(Int_t hvolt=1; hvolt < 6; hvolt++){
	sc[canindex]->cd(hvolt);
	if(hvolt==1){
	  t1->SetLineColor(1);
	  temptitle1.Form("C.cdetm1r.adc_c[%d] >> hh[%d]",paddleindex[paddle],canindex);
	  t1->Draw(temptitle1, hcuts[canindex]);
	  t1->SetLineColor(4);
	  temptitle2.Form("C.cdetm1r.adc_c[%d] >> ha[%d]",paddleindex[paddle],canindex);
	  t1->Draw(temptitle2, acuts[canindex], "same");
	  t1->SetLineColor(2);
	  temptitle3.Form("C.cdetm1r.adc_c[%d] >> hb[%d]",paddleindex[paddle],canindex);
	  t1->Draw(temptitle3, bcuts[canindex], "same");
	  t1->SetLineColor(3);
	  t1->SetFillColor(3);
	  temptitle4.Form("C.cdetm1r.adc_c[%d] >> hc[%d]",paddleindex[paddle],canindex);
	  t1->Draw(temptitle4, ccuts[canindex], "same");
	  t1->SetFillColor(0);
	  bf[canindex] = t1->GetEntries(bcuts[canindex].GetTitle());
	  af[canindex] = t1->GetEntries(ccuts[canindex].GetTitle());
	  eff[canindex] = (af[canindex])/(bf[canindex]);
	}   
	else if(hvolt==2){
	  t2->SetLineColor(1);
	  temptitle1.Form("C.cdetm1r.adc_c[%d] >> hh[%d]",paddleindex[paddle],canindex+14);
	  t2->Draw(temptitle1, hcuts[canindex]);
	  t2->SetLineColor(4);
	  temptitle2.Form("C.cdetm1r.adc_c[%d] >> ha[%d]",paddleindex[paddle],canindex+14);
	  t2->Draw(temptitle2, acuts[canindex], "same");
	  t2->SetLineColor(2);
	  temptitle3.Form("C.cdetm1r.adc_c[%d] >> hb[%d]",paddleindex[paddle],canindex+14);
	  t2->Draw(temptitle3, bcuts[canindex+14], "same");
	  t2->SetLineColor(3);
	  t2->SetFillColor(3);
	  temptitle4.Form("C.cdetm1r.adc_c[%d] >> hc[%d]",paddleindex[paddle],canindex+14);
	  t2->Draw(temptitle4, ccuts[canindex+14], "same");
	  t2->SetFillColor(0);
	  bf[canindex+14] = t2->GetEntries(bcuts[canindex+14].GetTitle());
	  af[canindex+14] = t2->GetEntries(ccuts[canindex+14].GetTitle());
	  eff[canindex+14] = (af[canindex+14])/(bf[canindex+14]);
	} 
	else if(hvolt==3){
	  t3->SetLineColor(1);
	  temptitle1.Form("C.cdetm1r.adc_c[%d] >> hh[%d]",paddleindex[paddle],canindex+28);
	  t3->Draw(temptitle1, hcuts[canindex]);
	  t3->SetLineColor(4);
	  temptitle2.Form("C.cdetm1r.adc_c[%d] >> ha[%d]",paddleindex[paddle],canindex+28);
	  t3->Draw(temptitle2, acuts[canindex], "same");
	  t3->SetLineColor(2);
	  temptitle3.Form("C.cdetm1r.adc_c[%d] >> hb[%d]",paddleindex[paddle],canindex+28);
	  t3->Draw(temptitle3, bcuts[canindex+28], "same");
	  t3->SetLineColor(3);
	  t3->SetFillColor(3);
	  temptitle4.Form("C.cdetm1r.adc_c[%d] >> hc[%d]",paddleindex[paddle],canindex+28);
	  t3->Draw(temptitle4, ccuts[canindex+28], "same");
	  t3->SetFillColor(0);
	  bf[canindex+28] = t3->GetEntries(bcuts[canindex+28].GetTitle());
	  af[canindex+28] = t3->GetEntries(ccuts[canindex+28].GetTitle());
	  eff[canindex+28] = (af[canindex+28])/(bf[canindex+28]);
	}
	else if(hvolt==4){
	  t4->SetLineColor(1);
	  temptitle1.Form("C.cdetm1r.adc_c[%d] >> hh[%d]",paddleindex[paddle],canindex+42);
	  t4->Draw(temptitle1, hcuts[canindex]);
	  t4->SetLineColor(4);
	  temptitle2.Form("C.cdetm1r.adc_c[%d] >> ha[%d]",paddleindex[paddle],canindex+42);
	  t4->Draw(temptitle2, acuts[canindex], "same");
	  t4->SetLineColor(2);
	  temptitle3.Form("C.cdetm1r.adc_c[%d] >> hb[%d]",paddleindex[paddle],canindex+42);
	  t4->Draw(temptitle3, bcuts[canindex+42], "same");
	  t4->SetLineColor(3);
	  t4->SetFillColor(3);
	  temptitle4.Form("C.cdetm1r.adc_c[%d] >> hc[%d]",paddleindex[paddle],canindex+42);
	  t4->Draw(temptitle4, ccuts[canindex+42], "same");
	  t4->SetFillColor(0);
	  bf[canindex+42] = t4->GetEntries(bcuts[canindex+42].GetTitle());
	  af[canindex+42] = t4->GetEntries(ccuts[canindex+42].GetTitle());
	  eff[canindex+42] = (af[canindex+42])/(bf[canindex+42]);
	}
	else{
	  t5->SetLineColor(1);
	  temptitle1.Form("C.cdetm1r.adc_c[%d] >> hh[%d]",paddleindex[paddle],canindex+56);
	  t5->Draw(temptitle1, hcuts[canindex]);
	  t5->SetLineColor(4);
	  temptitle2.Form("C.cdetm1r.adc_c[%d] >> ha[%d]",paddleindex[paddle],canindex+56);
	  t5->Draw(temptitle2, acuts[canindex], "same");
	  t5->SetLineColor(2);
	  temptitle3.Form("C.cdetm1r.adc_c[%d] >> hb[%d]",paddleindex[paddle],canindex+56);
	  t5->Draw(temptitle3, bcuts[canindex+56], "same");
	  t5->SetLineColor(3);
	  t5->SetFillColor(3);
	  temptitle4.Form("C.cdetm1r.adc_c[%d] >> hc[%d]",paddleindex[paddle],canindex+56);
	  t5->Draw(temptitle4, ccuts[canindex+56], "same");
	  t5->SetFillColor(0);
	  bf[canindex+56] = t5->GetEntries(bcuts[canindex+56].GetTitle());
	  af[canindex+56] = t5->GetEntries(ccuts[canindex+56].GetTitle());
	  eff[canindex+56] = (af[canindex+56])/(bf[canindex+56]);
	} 
      }
      canindex++;
    }
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
  outf << tx << "\t" << rn1 << "\t" << rn2 << "\t" << rn3 << "\t" << rn4 << "\t" << rn5 << endl;
  outf << "HV:" << "\t" << hv1 << "\t" << hv2 << "\t" << hv3 << "\t" << hv4 << "\t" << hv5 << endl;
  cout << "PMT " << pmt << " Efficiency" << endl;
  cout << tx << "\t" << rn1 << "\t" << rn2 << "\t" << rn3 << "\t" << rn4 << "\t" << rn5 << endl;
  cout << "HV:" << "\t" << hv1 << "\t" << hv2 << "\t" << hv3 << "\t" << hv4 << "\t" << hv5 << endl;

  Int_t effindex = 0;

  for(Int_t pid=0; pid < NUMPIXEL; pid++){
    if(pid!=pixel1[pmt-1]-1 && pid!=pixel2[pmt-1]-1){
      outf << pid+1 << "\t" << eff[effindex] << "\t" << eff[effindex+14] << "\t" << eff[effindex+28] << "\t" <<
	eff[effindex+42] << "\t" << eff[effindex+56] << endl;
      cout << pid+1 << "\t" << eff[effindex] << "\t" << eff[effindex+14] << "\t" << eff[effindex+28] << "\t" <<
	eff[effindex+42] << "\t" << eff[effindex+56] << endl;
      effindex++;
    }
  }

  outf.close();

  TString gname;
  TCanvas *ec[14];
  TGraph *gc1;
  TGraph *gc2;
  TGraph *gc3;
  TGraph *gc4;
  TGraph *gc5;
  TGraph *gc6;
  TGraph *gc7;
  TGraph *gc8;
  TGraph *gc9;
  TGraph *gc10;
  TGraph *gc11;
  TGraph *gc12;
  TGraph *gc13;
  TGraph *gc14;
  Int_t ecindex = 0;
  TF1 *fun1 = new TF1("fun1","pol3",700,800);
  fun1->SetLineColor(1);
  fun1->SetLineWidth(2);
  fun1->SetLineStyle(7);

  for(Int_t pind=0; pind < NUMPIXEL; pind++){
    if(pind!=pixel1[pmt-1]-1 && pind!=pixel2[pmt-1]-1){
      TString ecname = Form("ec%d",effindex+1);
      TString ectitle = Form("Pixel %d Eff",pind+1);
      ec[ecindex] = new TCanvas(ecname,ectitle,800,800);
      ec[ecindex]->cd();
      if(ecindex==14){
	break;
      }
      else if(ecindex==0){
	gc1->SetPoint(0,700,eff[ecindex]);
	gc1->SetPoint(1,725,eff[ecindex+14]);
	gc1->SetPoint(2,750,eff[ecindex+28]);
	gc1->SetPoint(3,775,eff[ecindex+42]);
	gc1->SetPoint(4,800,eff[ecindex+56]);
	gc1->Draw();
	gc1->Fit("fun1","B");
      }
      else if(ecindex==1){
	gc2->SetPoint(0,700,eff[ecindex]);
	gc2->SetPoint(1,725,eff[ecindex+14]);
	gc2->SetPoint(2,750,eff[ecindex+28]);
	gc2->SetPoint(3,775,eff[ecindex+42]);
	gc2->SetPoint(4,800,eff[ecindex+56]);
	gc2->Draw();
	gc2->Fit("fun1","B");
      }
      else if(ecindex==2){
	gc3->SetPoint(0,700,eff[ecindex]);
	gc3->SetPoint(1,725,eff[ecindex+14]);
	gc3->SetPoint(2,750,eff[ecindex+28]);
	gc3->SetPoint(3,775,eff[ecindex+42]);
	gc3->SetPoint(4,800,eff[ecindex+56]);
	gc3->Draw();
	gc3->Fit("fun1","B");
      }
      else if(ecindex==3){
	gc4->SetPoint(0,700,eff[ecindex]);
	gc4->SetPoint(1,725,eff[ecindex+14]);
	gc4->SetPoint(2,750,eff[ecindex+28]);
	gc4->SetPoint(3,775,eff[ecindex+42]);
	gc4->SetPoint(4,800,eff[ecindex+56]);
	gc4->Draw();
	gc4->Fit("fun1","B");
      }
      else if(ecindex==4){
	gc5->SetPoint(0,700,eff[ecindex]);
	gc5->SetPoint(1,725,eff[ecindex+14]);
	gc5->SetPoint(2,750,eff[ecindex+28]);
	gc5->SetPoint(3,775,eff[ecindex+42]);
	gc5->SetPoint(4,800,eff[ecindex+56]);
	gc5->Draw();
	gc5->Fit("fun1","B");
      }
      else if(ecindex==5){
	gc6->SetPoint(0,700,eff[ecindex]);
	gc6->SetPoint(1,725,eff[ecindex+14]);
	gc6->SetPoint(2,750,eff[ecindex+28]);
	gc6->SetPoint(3,775,eff[ecindex+42]);
	gc6->SetPoint(4,800,eff[ecindex+56]);
	gc6->Draw();
	gc6->Fit("fun1","B");
      }
      else if(ecindex==6){
	gc7->SetPoint(0,700,eff[ecindex]);
	gc7->SetPoint(1,725,eff[ecindex+14]);
	gc7->SetPoint(2,750,eff[ecindex+28]);
	gc7->SetPoint(3,775,eff[ecindex+42]);
	gc7->SetPoint(4,800,eff[ecindex+56]);
	gc7->Draw();
	gc7->Fit("fun1","B");
      }
      else if(ecindex==7){
	gc8->SetPoint(0,700,eff[ecindex]);
	gc8->SetPoint(1,725,eff[ecindex+14]);
	gc8->SetPoint(2,750,eff[ecindex+28]);
	gc8->SetPoint(3,775,eff[ecindex+42]);
	gc8->SetPoint(4,800,eff[ecindex+56]);
	gc8->Draw();
	gc8->Fit("fun1","B");
      }
      else if(ecindex==8){
	gc9->SetPoint(0,700,eff[ecindex]);
	gc9->SetPoint(1,725,eff[ecindex+14]);
	gc9->SetPoint(2,750,eff[ecindex+28]);
	gc9->SetPoint(3,775,eff[ecindex+42]);
	gc9->SetPoint(4,800,eff[ecindex+56]);
	gc9->Draw();
	gc9->Fit("fun1","B");
      }
      else if(ecindex==9){
	gc10->SetPoint(0,700,eff[ecindex]);
	gc10->SetPoint(1,725,eff[ecindex+14]);
	gc10->SetPoint(2,750,eff[ecindex+28]);
	gc10->SetPoint(3,775,eff[ecindex+42]);
	gc10->SetPoint(4,800,eff[ecindex+56]);
	gc10->Draw();
	gc10->Fit("fun1","B");
      }
      else if(ecindex==10){
	gc11->SetPoint(0,700,eff[ecindex]);
	gc11->SetPoint(1,725,eff[ecindex+14]);
	gc11->SetPoint(2,750,eff[ecindex+28]);
	gc11->SetPoint(3,775,eff[ecindex+42]);
	gc11->SetPoint(4,800,eff[ecindex+56]);
	gc11->Draw();
	gc11->Fit("fun1","B");
      }
      else if(ecindex==11){
	gc12->SetPoint(0,700,eff[ecindex]);
	gc12->SetPoint(1,725,eff[ecindex+14]);
	gc12->SetPoint(2,750,eff[ecindex+28]);
	gc12->SetPoint(3,775,eff[ecindex+42]);
	gc12->SetPoint(4,800,eff[ecindex+56]);
	gc12->Draw();
	gc12->Fit("fun1","B");
      }
      else if(ecindex==12){
	gc13->SetPoint(0,700,eff[ecindex]);
	gc13->SetPoint(1,725,eff[ecindex+14]);
	gc13->SetPoint(2,750,eff[ecindex+28]);
	gc13->SetPoint(3,775,eff[ecindex+42]);
	gc13->SetPoint(4,800,eff[ecindex+56]);
	gc13->Draw();
	gc13->Fit("fun1","B");
      }
      else if(ecindex==13){
	gc14->SetPoint(0,700,eff[ecindex]);
	gc14->SetPoint(1,725,eff[ecindex+14]);
	gc14->SetPoint(2,750,eff[ecindex+28]);
	gc14->SetPoint(3,775,eff[ecindex+42]);
	gc14->SetPoint(4,800,eff[ecindex+56]);
	gc14->Draw();
	gc14->Fit("fun1","B");
      }
      ecindex++;
    }
  }
  
  outf.close();
  /*
  TString rtitle;
  rtitle.Form("M1-R_Eff_PMT_%d.root",pmt);
  TFile *outroot = new TFile(rtitle,"NEW");

  for(Int_t r=0; r < NUMPADDLE; r++){
    outroot->WriteTObject(sc[r]);
    outroot->WriteTObject(ec[r]);
  }
  outroot->Write(outf);
  */
}

/*
for(Int_t r = 1; r<6; r++){
if (r == 1){
for (Int_t ent=1;ent<=numentries1;ent++){
           T->GetEntry(ent);
	  for (Int_t j=0; j<NUMPIXEL; j++){
	    Int_t k = (pmt-1)*NUMPIXEL+k;
	    if(bcuts[k])
	      {
		hb[j]->Fill(adc_c[k]);
		if(ccuts[k])
		  {
		    hc[j]->Fill(adc_c[k]);
		  }
	      }
  bf[j] = hb[j]->GetEntries();
  af[j] = hc[j]->GetEntries();

  eff[j] = af[j]/bf[j];

  cout << "# events before = " << bf[j] << endl;
  cout << "# events after = " << af[j] << endl;
  cout << " efficiency = " << eff[j] << endl;
	  }
	}
  bf[j] = hb[j]->GetEntries();
  af[j] = hc[j]->GetEntries();

  eff[j] = af[j]/bf[j];

  cout << "# events before = " << bf[j] << endl;
  cout << "# events after = " << af[j] << endl;
  cout << " efficiency = " << eff[j] << endl;
}
  else if (r == 2){

  }
  else if (r == 3){

  }
  else if (r == 4){

  }
  else{

  }
}
*/

  /*
  for(Int_t i=0; i<6; i++){

    if(i != mp1-1 && i != mp2-1){

      //  Int_t paddle = i+1

      Int_t indexup = pc+i+1;
      Int_t indexlw = pc+i-1;

      TCut tdccutpx = Form("C.cdetm1r.tdcl[%d]>500",pc+i);
      TCut tdccutlw = Form("C.cdetm1r.tdcl[%d]<500",pc+i-1);
      TCut tdccutup = Form("C.cdetm1r.tdcl[%d]<500",pc+i+1);

      TCut adccutlw = Form("C.cdetm1r.adc_c[%d]<20",indexlw);
      TCut adccutup = Form("C.cdetm1r.adc_c[%d]<20",indexup);    

      cout << i <<  "\t pixel = " << pc+i << "\t cut indices: " << pc+i-1 << "  " << pc+i+1 << endl;

      cout << "lower ADC cut: " << adccutlw.GetTitle() << endl;

      if(i+1 == mp1-1 || i+1 == mp2-1){
	TCut tdccutup = Form("C.cdetm1r.tdcl[%d]<500",pc+i+2);
	TCut adccutup = Form("C.cdetm1r.adc_c[%d]<20",pc+i+2);
      
	cout << "upper: " << adccutup.GetTitle() << endl;
      }

      else if(i-1 == mp1-1 || i-1 == mp2-1){
	TCut tdccutlw = Form("C.cdetm1r.tdcl[%d]<500",pc+i-2);
	TCut adccutlw = Form("C.cdetm1r.adc_c[%d]<20",pc+i-2);

	cout << "lower: " << adccutlw.GetTitle() << endl;
      }
    }
  }

  TCut tcxcut, ttcut;

  for(Int_t j=0; j<16; j++){

      if(j != i ){
	ttcut = Form("C.cdetm1r.tdcl[%d]<500",pc+j);
	tcxcut += ttcut;

	//	cout << j << "\t" << tcxcut << endl;
      }

      //      cout << tcxcut << endl;
  }

  TCut wdtdc = Form("(C.cdetm1r.tdcl[%d]-C.cdetm1r.tdct[%d]) > %f",pc+i,pc+i,tdcw[i]);

  cout << wdtdc.GetTitle() << endl;

  TCut cutthrs = Form("C.cdetm1r.adc_c[%d]>%f",pc+i,thr5[i]);

  cout << cutthrs << endl;

  TString cname = Form("sc%d",i+1);
  TString ctitle = Form("Pixel %d",i+1);
  sc[i+1] = new TCanvas(cname,ctitle,1200,700);
  sc[i+1]->Divide(5,2);
  sc[i+1]->Draw();

  TString hhname = Form("hh[%d]",i);
  TString hhtitle = Form("Pixel %d at 700 V",i+1);
  hh[i] = new TH1F(hhname,hhtitle,100,-100,400);

  TString haname = Form("ha[%d]",i);
  ha[i] = new TH1F(haname,"",100,-100,400);
  ha[i]->SetLineColor(2);

  TString hbname = Form("hb[%d]",i);
  hb[i] = new TH1F(hbname,"",100,-100,400);
  hb[i]->SetLineColor(3);

  TString hcname = Form("hc[%d]",i);
  hc[i] = new TH1F(hcname,"",100,-100,400);
  hc[i]->SetLineColor(3);
  hc[i]->SetFillColor(3);

  TString adc = Form("C.cdetm1r.adc_c[%d] >> hh[%d]",pc+i,i);
  TString adccx = Form("C.cdetm1r.adc_c[%d] >> ha[%d]",pc+i,i);
  TString adcwd = Form("C.cdetm1r.adc_c[%d] >> hb[%d]",pc+i,i);
  TString adcthr = Form("C.cdetm1r.adc_c[%d] >> hc[%d]",pc+i,i);
    
  sc[i+1]->cd(5);

  t5->Draw(adc,tdccutpx && tdccutlw && tdccutup && adccutlw && adccutup);

  t5->Draw(adccx,tdccutpx && tdccutlw && tdccutup && adccutlw && adccutup && tcxcut,"same");

  t5->Draw(adcwd,tdccutpx && tdccutlw && tdccutup && adccutlw && adccutup && wdtdc && tcxcut,"same");

  t5->Draw(adcthr,tdccutpx && tdccutlw && tdccutup && adccutlw && adccutup && wdtdc && tcxcut && cutthrs,"same");
    
  bf[i] = hb[i]->GetEntries();
  af[i] = hc[i]->GetEntries();

  eff[i] = af[i]/bf[i];

  cout << "# events before = " << bf[i] << endl;
  cout << "# events after = " << af[i] << endl;
  cout << " efficiency = " << eff[i] << endl;

    // write the efficiency values to an output file

  outf << i+1 << "\t" << eff[i] << endl;

  }
  // end of 'i' for loop

  // Create all the cuts for each pixel and make the adc plots

  //cout << "end of script here" << endl;
  
  // going to create the plots for each pixel in turn

  TCanvas *s1 = new TCanvas("s1","pixel 1 - 700 V",700,900);
  s1->cd();
  s1->Divide(1,2);

  s1->cd(1);

  // t1->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]","C.cdetm1r.tdcl[96]>500");
  // tg->SetMarkerStyle(20);
  // tg->SetMarkerSize(0.7);

  TCut tdccut = Form("C.cdetm1r.tdcl[%d]>500 && C.cdetm1r.tdcl[%d]<500",pc,pc-1);

  TCut adccut = Form("C.cdetm1r.adc_c[%d]<20 && C.cdetm1r.adc_c[%d]<20",pc-1, pc+1);

  // cout << adccut << "\n" << endl;
  // cout << tdccut << "\n" << endl;
  
  TCut totalcut = adccut && tdccut;

  // cout << totalcut << endl;

  
  //  cout << thrscut << endl;

  // // s1->Update();

  // s1->cd(2);

  //  t1->Draw("C.cdetm1r.adc_c[96]",adccut);
  TH1F *hh = new TH1F("hh","ADC PMT 7, pixel 1;ADC channels; Counts",100,-100,400);
  TH1F *ha = new TH1F("ha","",100,-100,400);
  ha->SetLineColor(2);
  TH1F *hb = new TH1F("hb","",100,-100,400);
  hb->SetLineColor(3);

  t1->Draw("C.cdetm1r.adc_c[96] >> hh",totalcut);
  //  hh->SetName("aa");
  // apply the tdc width cut

  t1->SetLineColor(2);
  t1->Draw("C.cdetm1r.adc_c[96] >> ha",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  cout << ha->GetEntries() << endl;

  t1->SetLineColor(3);
  t1->Draw("C.cdetm1r.adc_c[96] >> hb",(totalcut && thrscut1)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");  // width, crosstalk and threshold cut
  cout << hb->GetEntries() << endl;
  //  htemp->SetName("bb"):
  
  // Float_t eff = 0.0;
  // cout << eff = hb->GetEntries()/ha->GetEntries() << endl;

  // Now plot the leading versus trailing tdc
  s1->cd(2);
  t1->SetMarkerStyle(21);
  t1->SetMarkerSize(0.7);
  t1->SetMarkerColor(1);

  t1->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96] >> h1",totalcut);

  t1->SetMarkerColor(2);
  t1->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  

  t1->SetMarkerColor(3);
  t1->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",(totalcut && thrscut1)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");


  //------------------------------------------------------------------

  // Now repeat only for this pixel but for the other HV settings
  // HV = 725 V

  TCanvas *s2 = new TCanvas("s2","pixel 1 - 725 V",700,900);
  s2->cd();
  s2->Divide(1,2);

  s2->cd(1);



  //  t1->Draw("C.cdetm1r.adc_c[96]",adccut);
  TH1F *jh = new TH1F("jh","ADC PMT 7, pixel 1;ADC channels; Counts",100,-100,400);
  TH1F *ja = new TH1F("ja","",100,-100,400);
  ja->SetLineColor(2);
  TH1F *jb = new TH1F("jb","",100,-100,400);
  jb->SetLineColor(3);

  t2->Draw("C.cdetm1r.adc_c[96] >> jh",totalcut);
  //  hh->SetName("aa");
  // apply the tdc width cut

  t2->SetLineColor(2);
  t2->Draw("C.cdetm1r.adc_c[96] >> ja",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  cout << ja->GetEntries() << endl;

  t2->SetLineColor(3);
  t2->Draw("C.cdetm1r.adc_c[96] >> jb",(totalcut && thrscut2)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");  // width, crosstalk and threshold cut
  cout << jb->GetEntries() << endl;

  // Now plot the leading versus trailing tdc
  s2->cd(2);
  t2->SetMarkerStyle(21);
  t2->SetMarkerSize(0.7);
  t2->SetMarkerColor(1);

  t2->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96] >> h1",totalcut);

  t2->SetMarkerColor(2);
  t2->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  

  t2->SetMarkerColor(3);
  t2->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",(totalcut && thrscut2)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");

  
  //--------------------------------------------------------------
  // HV = 750 V

  TCanvas *s3 = new TCanvas("s3","pixel 1 - 750 V",700,900);
  s3->cd();
  s3->Divide(1,2);

  s3->cd(1);


  TH1F *kh = new TH1F("kh","ADC PMT 7, pixel 1;ADC channels; Counts",100,-100,400);
  TH1F *ka = new TH1F("ka","",100,-100,400);
  ka->SetLineColor(2);
  TH1F *kb = new TH1F("kb","",100,-100,400);
  kb->SetLineColor(3);

  t3->Draw("C.cdetm1r.adc_c[96] >> kh",totalcut);
  //  hh->SetName("aa");
  // apply the tdc width cut

  t3->SetLineColor(2);
  t3->Draw("C.cdetm1r.adc_c[96] >> ka",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  cout << ka->GetEntries() << endl;

  t3->SetLineColor(3);
  t3->Draw("C.cdetm1r.adc_c[96] >> kb",(totalcut && thrscut3)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");  // width, crosstalk and threshold cut
  cout << kb->GetEntries() << endl;

  s3->cd(2);
  t3->SetMarkerStyle(21);
  t3->SetMarkerSize(0.7);
  t3->SetMarkerColor(1);

  t3->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96] >> h1",totalcut);

  t3->SetMarkerColor(2);
  t3->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  

  t3->SetMarkerColor(3);
  t3->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",(totalcut && thrscut3)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");


  //---------------------------------------------------
  // HV = 775 V

  TCanvas *s4 = new TCanvas("s4","pixel 1 - 775 V",700,900);
  s4->cd();
  s4->Divide(1,2);

  s4->cd(1);



  //  t1->Draw("C.cdetm1r.adc_c[96]",adccut);
  TH1F *lh = new TH1F("lh","ADC PMT 7, pixel 1;ADC channels; Counts",100,-100,400);
  TH1F *la = new TH1F("la","",100,-100,400);
  la->SetLineColor(2);
  TH1F *lb = new TH1F("lb","",100,-100,400);
  lb->SetLineColor(3);

  t4->Draw("C.cdetm1r.adc_c[96] >> lh",totalcut);
  //  hh->SetName("aa");
  // apply the tdc width cut

  t4->SetLineColor(2);
  t4->Draw("C.cdetm1r.adc_c[96] >> la",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  cout << la->GetEntries() << endl;

  t4->SetLineColor(3);
  t4->Draw("C.cdetm1r.adc_c[96] >> lb",(totalcut && thrscut4)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");  // width, crosstalk and threshold cut
  cout << lb->GetEntries() << endl;

  s4->cd(2);
  t4->SetMarkerStyle(21);
  t4->SetMarkerSize(0.7);
  t4->SetMarkerColor(1);

  t4->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96] >> h1",totalcut);

  t4->SetMarkerColor(2);
  t4->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  

  t4->SetMarkerColor(3);
  t4->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",(totalcut && thrscut4)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");

  //---------------------------------------------------------
  // HV = 800 V

  TCanvas *s5 = new TCanvas("s5","pixel 1 - 800 V",700,900);
  s5->cd();
  s5->Divide(1,2);

  s5->cd(1);



  //  t1->Draw("C.cdetm1r.adc_c[96]",adccut);
  TH1F *mh = new TH1F("mh","ADC PMT 7, pixel 1;ADC channels; Counts",100,-100,400);
  TH1F *ma = new TH1F("ma","",100,-100,400);
  ma->SetLineColor(2);
  TH1F *mb = new TH1F("mb","",100,-100,400);
  mb->SetLineColor(3);

  t5->Draw("C.cdetm1r.adc_c[96] >> mh",totalcut);
  //  hh->SetName("aa");
  // apply the tdc width cut

  t5->SetLineColor(2);
  t5->Draw("C.cdetm1r.adc_c[96] >> ma",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>50","same");
  cout << ma->GetEntries() << endl;

  t5->SetLineColor(3);
  t5->Draw("C.cdetm1r.adc_c[96] >> mb",(totalcut && thrscut5)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>50","same");  // width, crosstalk and threshold cut
  cout << mb->GetEntries() << endl;
  
  s5->cd(2);
  t5->SetMarkerStyle(21);
  t5->SetMarkerSize(0.7);
  t5->SetMarkerColor(1);

  t5->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96] >> h1",totalcut);

  t5->SetMarkerColor(2);
  t5->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",totalcut*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");
  

  t5->SetMarkerColor(3);
  t5->Draw("C.cdetm1r.tdcl[96]:C.cdetm1r.tdct[96]",(totalcut && thrscut5)*"(C.cdetm1r.tdcl[96]-C.cdetm1r.tdct[96])>40","same");


  //-------------------------------------------------------

  // Now create the efficiency plot and check

  TCanvas *s6 = new TCanvas("s6","efficiency",800,600);
  s6->cd();

  TH2F *ef1 = new TH2F("ef2","; High Voltage [V]; Efficiency [%]",100,690,810,100,90.0,100.0);
  // ef1->SetMarkerStyle(22);
  // ef1->SetMarkerColor(2);
  // ef1->SetMarkerSize(0.9);
  ef1->Draw();

  TGraph *gr1 = new TGraph();
  gr1->SetName("gr1");
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(7);
  //  gr1->SetMarkerSize(0.9);

  // fill the arrays with the numbers of events

  bf[0] = ha->GetEntries();
  bf[1] = ja->GetEntries();
  bf[2] = ka->GetEntries();
  bf[3] = la->GetEntries();
  bf[4] = ma->GetEntries();

  af[0] = hb->GetEntries();
  af[1] = jb->GetEntries();
  af[2] = kb->GetEntries();
  af[3] = lb->GetEntries();
  af[4] = mb->GetEntries();

  for(Int_t i=0; i<5; i++){
    gr1->SetPoint(i,700+25*i,100*af[i]/bf[i]);

    cout << i << "\t" << af[i] << "\t" << 100*af[i]/bf[i] << endl;
  } 
  
  gr1->Draw("same P");

  s6->Print("eff-test.pdf");
























 /*
 
  TH2F *h1 = new TH2F("h1","TDC comparison;Leading Edge; Trailing Edge",100,800,950,100,800,950);
  h1->SetMarkerStyle(21);
  h1->SetMarkerColor(1);
  h1->SetMarkerSize(0.7);


  //  cout << ff->GetEntries() << endl;




  /*
  Float_t bfm = 0.0;
  Float_t afm = 0.0;  // to calculate the mean of each set

  // Define a graph for each set of values
  
  TGraphErrors *g1 = new TGraphErrors();  // before
  g1->SetName("g1");
  g1->SetMarkerStyle(22);       // up triangles
  g1->SetMarkerColor(kBlack);   // color = black
  g1->SetLineColor(kBlack);
  g1->SetMarkerSize(1.2);
  
  TGraphErrors *g2 = new TGraphErrors();  // after
  g2->SetName("g2");
  g2->SetMarkerStyle(20);       // circles
  g2->SetMarkerColor(kGreen+2);   // color = green
  g2->SetLineColor(kGreen+2);
  g2->SetMarkerSize(1.2);

  // Use a loop to set the values in each graph

  for(Int_t i=0; i<14; i++){
    
    g1->SetPoint(i,i+1,abef[i]);
    g1->SetPointError(i,0.0,0.0);

    g2->SetPoint(i,i+1,aaft[i]);
    g2->SetPointError(i,0.0,erraft[i]);

    // calculate the means

    bfm += abef[i];
    afm += aaft[i];

  }

  // ------------------------------------------------

  // Create the canvas and histogram and draw the graphs

  TCanvas *s1 = new TCanvas("s1","adc-compare",900,700);
  s1->cd();

  TH1F *mdf = new TH1F("mdf","Comparison of Efficiency Ratio - PMT 7",15,0,15);
  mdf->SetMaximum(150.0);
  mdf->SetMinimum(0.0);
  mdf->SetStats(0);

  mdf->GetYaxis()->SetTitle("50% ADC Value");
  mdf->GetXaxis()->SetTitle("Paddle Number");
  mdf->SetTitleOffset(1.55,"y");
  mdf->Draw();

  // now draw the two graphs on the same plot

  g1->Draw("same P");

  g2->Draw("same P");

  // add in a legend for this plot

  lgto = new TLegend(0.63,0.7,0.88,0.85,"");
  lgto->AddEntry(g1, "CAPACITORS", "P");
  lgto->AddEntry(g2, "NO CAPACITORS", "P");
  lgto->SetTextFont(42);
  lgto->SetTextSize(0.025);
  lgto->SetFillStyle(0);
  lgto->SetBorderSize(1);
  lgto->Draw();

  // add in two lines as eye-guides
  // use the mean of each distribution for numbers

  /*
  TLine *lg1 = new TLine(0.0,bfm/14 ,15.0, bfm/14);
  lg1->SetLineColor(kRed+2);
  lg1->SetLineStyle(2);
  lg1->SetLineWidth(3);
  lg1->Draw();

  TLine *lg2 = new TLine(0.0,afm/14 ,15.0, afm/14);
  lg2->SetLineColor(kGreen+1);
  lg2->SetLineStyle(2);
  lg2->SetLineWidth(3);
  lg2->Draw();
  */



  // end of script here

//}
