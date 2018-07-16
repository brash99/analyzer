#define carbon_chain_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

#include "carbon_chain.h"

using namespace std;

void carbon_chain::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L carbon_chain.C
//      Root > carbon_chain t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("CT.ct_1by2",1);
   fChain->SetBranchStatus("CT.ct_2by1",1);
   fChain->SetBranchStatus("CT.ncomb",1);
   fChain->SetBranchStatus("CT.ntr1",1);
   fChain->SetBranchStatus("CT.ntr2",1);
 
   fChain->SetBranchStatus("fEvtHdr.fEvtType",1);

   fChain->SetBranchStatus("EK_L.Q2",1);
   fChain->SetBranchStatus("EK_L.W2",1);
   fChain->SetBranchStatus("EK_L.angle",1);
   fChain->SetBranchStatus("EK_L.epsilon",1);
   fChain->SetBranchStatus("EK_L.nu",1);
   fChain->SetBranchStatus("EK_L.omega",1);
   fChain->SetBranchStatus("EK_L.ph_q",1);
   fChain->SetBranchStatus("EK_L.q3m",1);
   fChain->SetBranchStatus("EK_L.th_q",1);
   fChain->SetBranchStatus("EK_L.x_bj",1);

   fChain->SetBranchStatus("L.gold.dp",1);
   fChain->SetBranchStatus("L.gold.index",1);
   fChain->SetBranchStatus("L.gold.ok",1);
   fChain->SetBranchStatus("L.gold.p",1);
   fChain->SetBranchStatus("L.gold.ph",1);
   fChain->SetBranchStatus("L.gold.px",1);
   fChain->SetBranchStatus("L.gold.py",1);
   fChain->SetBranchStatus("L.gold.pz",1);
   fChain->SetBranchStatus("L.gold.th",1);
   fChain->SetBranchStatus("L.gold.x",1);
   fChain->SetBranchStatus("L.gold.y",1);

   fChain->SetBranchStatus("NA.nd.p1.la",1);
   fChain->SetBranchStatus("NA.nd.p1.la_c",1);
   fChain->SetBranchStatus("NA.nd.p1.lt",1);
   fChain->SetBranchStatus("NA.nd.p1.lt_c",1);
   fChain->SetBranchStatus("NA.nd.p1.ra",1);
   fChain->SetBranchStatus("NA.nd.p1.ra_c",1);
   fChain->SetBranchStatus("NA.nd.p1.rt",1);
   fChain->SetBranchStatus("NA.nd.p1.rt_c",1);
  
   fChain->SetBranchStatus("NA.nd.p2.la",1);
   fChain->SetBranchStatus("NA.nd.p2.la_c",1);
   fChain->SetBranchStatus("NA.nd.p2.lt",1);
   fChain->SetBranchStatus("NA.nd.p2.lt_c",1);
   fChain->SetBranchStatus("NA.nd.p2.ra",1);
   fChain->SetBranchStatus("NA.nd.p2.ra_c",1);
   fChain->SetBranchStatus("NA.nd.p2.rt",1);
   fChain->SetBranchStatus("NA.nd.p2.rt_c",1);

   fChain->SetBranchStatus("NA.nd.p3.la",1);
   fChain->SetBranchStatus("NA.nd.p3.la_c",1);
   fChain->SetBranchStatus("NA.nd.p3.lt",1);
   fChain->SetBranchStatus("NA.nd.p3.lt_c",1);
   fChain->SetBranchStatus("NA.nd.p3.ra",1);
   fChain->SetBranchStatus("NA.nd.p3.ra_c",1);
   fChain->SetBranchStatus("NA.nd.p3.rt",1);
   fChain->SetBranchStatus("NA.nd.p3.rt_c",1);

   fChain->SetBranchStatus("NA.nd.p4.la",1);
   fChain->SetBranchStatus("NA.nd.p4.la_c",1);
   fChain->SetBranchStatus("NA.nd.p4.lt",1);
   fChain->SetBranchStatus("NA.nd.p4.lt_c",1);
   fChain->SetBranchStatus("NA.nd.p4.ra",1);
   fChain->SetBranchStatus("NA.nd.p4.ra_c",1);
   fChain->SetBranchStatus("NA.nd.p4.rt",1);
   fChain->SetBranchStatus("NA.nd.p4.rt_c",1);

   fChain->SetBranchStatus("NA.nd.ntracks",1);
   fChain->SetBranchStatus("NA.nd.p1.npadhit",1);
   fChain->SetBranchStatus("NA.nd.p2.npadhit",1);
   fChain->SetBranchStatus("NA.nd.p3.npadhit",1);
   fChain->SetBranchStatus("NA.nd.p4.npadhit",1);
   fChain->SetBranchStatus("NA.ntracks",1);
   fChain->SetBranchStatus("NA.veto.npadhit",1);
   fChain->SetBranchStatus("NA.wov.ntracks",1);
   fChain->SetBranchStatus("NA.wv.ntracks",1);
 
   fChain->SetBranchStatus("SK.Erecoil",1);
   fChain->SetBranchStatus("SK.MandelS",1);
   fChain->SetBranchStatus("SK.MandelT",1);
   fChain->SetBranchStatus("SK.MandelU",1);
   fChain->SetBranchStatus("SK.Mrecoil",1);
   fChain->SetBranchStatus("SK.emiss",1);
   fChain->SetBranchStatus("SK.ph_bq",1);
   fChain->SetBranchStatus("SK.ph_xq",1);
   fChain->SetBranchStatus("SK.phb_cm",1);
   fChain->SetBranchStatus("SK.phx_cm",1);
   fChain->SetBranchStatus("SK.pmiss",1);
   fChain->SetBranchStatus("SK.pmiss_x",1);
   fChain->SetBranchStatus("SK.pmiss_y",1);
   fChain->SetBranchStatus("SK.pmiss_z",1);
   fChain->SetBranchStatus("SK.px_cm",1);
   fChain->SetBranchStatus("SK.t_tot_cm",1);
   fChain->SetBranchStatus("SK.tb",1);
   fChain->SetBranchStatus("SK.tb_cm",1);
   fChain->SetBranchStatus("SK.th_bq",1);
   fChain->SetBranchStatus("SK.th_xq",1);
   fChain->SetBranchStatus("SK.thb_cm",1);
   fChain->SetBranchStatus("SK.thx_cm",1);
   fChain->SetBranchStatus("SK.tx",1);
   fChain->SetBranchStatus("SK.tx_cm",1);
   fChain->SetBranchStatus("SK.xangle",1);

                                                                             
   TFile *out = new TFile("output_chain_carbon_T3.root","RECREATE");
   TTree t1("t1","T3 events only tree");
                                                                               
   // Only worrying about the first elements in the coincidence-time
   // arrays, corresponding to the 'first' track in each spectrometer.
                                                                               
   t1.Branch("CT.ncomb",&CT_ncomb,"CT.ncomb/D");
   t1.Branch("CT.ntr1",&CT_ntr1,"CT.ntr1/D");
   t1.Branch("CT.ntr2",&CT_ntr2,"CT.ntr2/D");
   t1.Branch("CT.ct_1by2",CT_ct_1by2,"CT.ct_1by2/D");
   t1.Branch("CT.ct_2by1",CT_ct_2by1,"CT.ct_2by1/D");

   t1.Branch("fEvtHdr.fEvtType",&fEvtHdr_fEvtType,"fEvtHdr.fEvtType/I");

   t1.Branch("L.gold.dp",&L_gold_dp,"L.gold.dp/D");
   t1.Branch("L.gold.index",&L_gold_index,"L.gold.index/D");
   t1.Branch("L.gold.ok",&L_gold_ok,"L.gold.ok/D");
   t1.Branch("L.gold.p",&L_gold_p,"L.gold.p/D");
   t1.Branch("L.gold.ph",&L_gold_ph,"L.gold.ph/D");
   t1.Branch("L.gold.px",&L_gold_px,"L.gold.px/D");
   t1.Branch("L.gold.py",&L_gold_py,"L.gold.py/D");
   t1.Branch("L.gold.pz",&L_gold_pz,"L.gold.pz/D");
   t1.Branch("L.gold.th",&L_gold_th,"L.gold.th/D");
   t1.Branch("L.gold.x",&L_gold_x,"L.gold.x/D");
   t1.Branch("L.gold.y",&L_gold_y,"L.gold.y/D");
 
   t1.Branch("EK_L.Q2",&EK_L_Q2,"EK_L.Q2/D");
   t1.Branch("EK_L.W2",&EK_L_W2,"EK_L.W2/D");
   t1.Branch("EK_L.angle",&EK_L_angle,"EK_L.angle/D");
   t1.Branch("EK_L.epsilon",&EK_L_epsilon,"EK_L.epsilon/D");
   t1.Branch("EK_L.nu",&EK_L_nu,"EK_L.nu/D");
   t1.Branch("EK_L.omega",&EK_L_omega,"EK_L.omega/D");
   t1.Branch("EK_L.ph_q",&EK_L_ph_q,"EK_L.ph_q/D");
   t1.Branch("EK_L.q3m",&EK_L_q3m,"EK_L.q3m/D");
   t1.Branch("EK_L.th_q",&EK_L_th_q,"EK_L.th_q/D");
   t1.Branch("EK_L.x_bj",&EK_L_x_bj,"EK_L.x_bj/D");

   t1.Branch("NA.nd.p1.la",&NA_nd_p1_la,"NA.nd.p1.la/D");
   t1.Branch("NA.nd.p1.la_c",&NA_nd_p1_la_c,"NA.nd.p1.la_c/D");
   t1.Branch("NA.nd.p1.lt",&NA_nd_p1_lt,"NA.nd.p1.lt/D");
   t1.Branch("NA.nd.p1.lt_c",&NA_nd_p1_lt_c,"NA.nd.p1.lt_c/D");
   t1.Branch("NA.nd.p1.ra",&NA_nd_p1_ra,"NA.nd.p1.ra/D");
   t1.Branch("NA.nd.p1.ra_c",&NA_nd_p1_ra_c,"NA.nd.p1.ra_c/D");
   t1.Branch("NA.nd.p1.rt",&NA_nd_p1_rt,"NA.nd.p1.rt/D");
   t1.Branch("NA.nd.p1.rt_c",&NA_nd_p1_rt_c,"NA.nd.p1.rt_c/D");
  
   t1.Branch("NA.nd.p2.la",&NA_nd_p2_la,"NA.nd.p2.la/D");
   t1.Branch("NA.nd.p2.la_c",&NA_nd_p2_la_c,"NA.nd.p2.la_c/D");
   t1.Branch("NA.nd.p2.lt",&NA_nd_p2_lt,"NA.nd.p2.lt/D");
   t1.Branch("NA.nd.p2.lt_c",&NA_nd_p2_lt_c,"NA.nd.p2.lt_c/D");
   t1.Branch("NA.nd.p2.ra",&NA_nd_p2_ra,"NA.nd.p2.ra/D");
   t1.Branch("NA.nd.p2.ra_c",&NA_nd_p2_ra_c,"NA.nd.p2.ra_c/D");
   t1.Branch("NA.nd.p2.rt",&NA_nd_p2_rt,"NA.nd.p2.rt/D");
   t1.Branch("NA.nd.p2.rt_c",&NA_nd_p2_rt_c,"NA.nd.p2.rt_c/D");

   t1.Branch("NA.nd.p3.la",&NA_nd_p3_la,"NA.nd.p3.la/D");
   t1.Branch("NA.nd.p3.la_c",&NA_nd_p3_la_c,"NA.nd.p3.la_c/D");
   t1.Branch("NA.nd.p3.lt",&NA_nd_p3_lt,"NA.nd.p3.lt/D");
   t1.Branch("NA.nd.p3.lt_c",&NA_nd_p3_lt_c,"NA.nd.p3.lt_c/D");
   t1.Branch("NA.nd.p3.ra",&NA_nd_p3_ra,"NA.nd.p3.ra/D");
   t1.Branch("NA.nd.p3.ra_c",&NA_nd_p3_ra_c,"NA.nd.p3.ra_c/D");
   t1.Branch("NA.nd.p3.rt",&NA_nd_p3_rt,"NA.nd.p3.rt/D");
   t1.Branch("NA.nd.p3.rt_c",&NA_nd_p3_rt_c,"NA.nd.p3.rt_c/D");

   t1.Branch("NA.nd.p4.la",&NA_nd_p4_la,"NA.nd.p4.la/D");
   t1.Branch("NA.nd.p4.la_c",&NA_nd_p4_la_c,"NA.nd.p4.la_c/D");
   t1.Branch("NA.nd.p4.lt",&NA_nd_p4_lt,"NA.nd.p4.lt/D");
   t1.Branch("NA.nd.p4.lt_c",&NA_nd_p4_lt_c,"NA.nd.p4.lt_c/D");
   t1.Branch("NA.nd.p4.ra",&NA_nd_p4_ra,"NA.nd.p4.ra/D");
   t1.Branch("NA.nd.p4.ra_c",&NA_nd_p4_ra_c,"NA.nd.p4.ra_c/D");
   t1.Branch("NA.nd.p4.rt",&NA_nd_p4_rt,"NA.nd.p4.rt/D");
   t1.Branch("NA.nd.p4.rt_c",&NA_nd_p4_rt_c,"NA.nd.p4.rt_c/D");

   /*   t1.Branch("NA.nd.ntracks",&NA.nd.ntracks,"NA.nd.ntracks/D");
   t1.Branch("NA.nd.p1.npadhit",&NA.nd.p1.npadhit,"NA.nd.p1.npadhit/D");
   t1.Branch("NA.nd.p2.npadhit",&NA.nd.p2.npadhit,"NA.nd.p2.npadhit/D");
   t1.Branch("NA.nd.p3.npadhit",&NA.nd.p3.npadhit,"NA.nd.p3.npadhit/D");
   t1.Branch("NA.nd.p4.npadhit",&NA.nd.p4.npadhit,"NA.nd.p4.npadhit/D");
   t1.Branch("NA.ntracks",&NA.ntracks,"NA.ntracks/D");
   t1.Branch("NA.veto.npadhit",&NA.veto.npadhit,"NA.veto.npadhit/D");
   t1.Branch("NA.wov.ntracks",&NA.wov.ntracks,"NA.wov.ntracks/D");
   t1.Branch("NA.wv.ntracks",&NA.wv.ntracks,"NA.wv.ntracks/D");
   */

   t1.Branch("SK.Erecoil",&SK_Erecoil,"SK.Erecoil/D");
   t1.Branch("SK.MandelS",&SK_MandelS,"SK.MandelS/D");
   t1.Branch("SK.MandelT",&SK_MandelT,"SK.MandelT/D");
   t1.Branch("SK.MandelU",&SK_MandelU,"SK.MandelU/D");
   t1.Branch("SK.Mrecoil",&SK_Mrecoil,"SK.Mrecoil/D");
   t1.Branch("SK.emiss",&SK_emiss,"SK.emiss/D");
   t1.Branch("SK.ph_bq",&SK_ph_bq,"SK.ph_bq/D");
   t1.Branch("SK.ph_xq",&SK_ph_xq,"SK.ph_xq/D");
   t1.Branch("SK.phb_cm",&SK_phb_cm,"SK.phb_cm/D");
   t1.Branch("SK.phx_cm",&SK_phx_cm,"SK.phx_cm/D");
   t1.Branch("SK.pmiss",&SK_pmiss,"SK.pmiss/D");
   t1.Branch("SK.pmiss_x",&SK_pmiss_x,"SK.pmiss_x/D");
   t1.Branch("SK.pmiss_y",&SK_pmiss_y,"SK.pmiss_y/D");
   t1.Branch("SK.pmiss_z",&SK_pmiss_z,"SK.pmiss_z/D");
   t1.Branch("SK.px_cm",&SK_px_cm,"SK.px_cm/D");
   t1.Branch("SK.t_tot_cm",&SK_t_tot_cm,"SK.t_tot_cm/D");
   t1.Branch("SK.tb",&SK_tb,"SK.tb/D");
   t1.Branch("SK.tb_cm",&SK_tb_cm,"SK.tb_cm/D");
   t1.Branch("SK.th_bq",&SK_th_bq,"SK.th_bq/D");
   t1.Branch("SK.th_xq",&SK_th_xq,"SK.th_xq/D");
   t1.Branch("SK.thb_cm",&SK_thb_cm,"SK.thb_cm/D");
   t1.Branch("SK.thx_cm",&SK_thx_cm,"SK.thx_cm/D");
   t1.Branch("SK.tx",&SK_tx,"SK.tx/D");
   t1.Branch("SK.tx_cm",&SK_tx_cm,"SK.tx_cm/D");
   t1.Branch("SK.xangle",&SK_xangle,"SK.xangle/D");



   Int_t nentries = Int_t(fChain->GetEntriesFast());

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry % 1000 == 0) cout << "." << flush;
      if (Cut(ientry) < 0) continue;
      //  if (fEvtHdr.fEvtType==3)
	t1.Fill();
   }
   cout << "I got to here! 1 " << endl;
   out->cd();
   cout << "I got to here! 2 " << endl;
   t1.Write();
   cout << "I got to here! 3 " << endl;
   //   out->Write();
}
