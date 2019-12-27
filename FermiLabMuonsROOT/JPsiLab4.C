//This is the main loop used to interact with the ROOT 
//software at CERN and identify particles from the tracks they leave 
//in the detector.

//It consists of a nested for loops. If the particle we are looking for 
//decays into three known products, we would use three nested for loops for each
//detectable particle decay. 

#define JPsiLab_cxx
#include "JPsiLab.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include "LsqFit.hh"
#include "LsqFit.cc"
double MU_Mass_Gev = 0.105658369;
double K_Mass = 0.493677;
double Mu_Mass_Gev = 0.105658369;
double K_mass = 0.493677;
//globals, known masses 
void JPsiLab::Loop()
{
  //---
  //--- The next section contains details of creating the output tree
  //---
  TFile *f2 = new TFile("phimble.root", "RECREATE");  //make a new root file
  f2->cd();
  TTree *output = new TTree("output","calculated values");  //make a new tree inside the file

  struct outvar {         //here we are creating data structures of new variables.
    Float_t jpmass;       //having an existing structure makes it MUCH simpler to create new branches
    
  };                      //and leaves in the tree.  the variables will be filled and written to the
  outvar out;             //tree later, but for now we just create space in memory for them.

  struct P_struct {
    Float_t B_mass;
    Float_t phimass;
    double vx;
    double vy;
  };
  P_struct Part_Struct;
  struct vtxvar {
    Float_t xintersect;
    Float_t yintersect;
    Float_t primx;
    Float_t primy;
    Float_t jpmomz;
    Float_t jpmomt;
    Float_t jpcosth;
    Float_t pthdecay;
  };
  vtxvar outvtx;
 
  //make new branches and leaves in the "output" tree:
  TBranch *jpmass = output->Branch("jpmass",&out.jpmass,"jpmass/F"); 
  TBranch *Bmeson = output->Branch("BMass",&Part_Struct.B_mass,"B_mass/F:phimass/F");
  TBranch *vertex = output->Branch("vertex",&outvtx.xintersect,"xintersect/F:yintersect:primx:primy:jpmomz/F:jpmomt/F:jpcosth/F:pthdecay/F:vx/d:vy/d");
  if (jpmass) ;    // This just avoids a 'variable not used' message
  if (vertex) ;    // so does this
  if (Bmeson) ;

  //---
  //--- Create new histograms here, by copy-pasting the example
  //--- below and changing the jpm variable to something unique
  //---
  //TH1F *hMOMpt = new TH1F("hMOMpt", "Transverse momentum of muons",100,0.0,5.20);
  //TH1F *hMOMpz = new TH1F("hMOMpz", "Parallel momentum of muons",100,-5.20,5.20);
  //TH1F *hMOMth = new TH1F("hMOMth", "Angle from z of muon momentum",100,-1,1);
  //TH1F *hMOMthp = new TH1F("hMOMthp", "Angle from z of Positive muon momentum",100,-1,1);
  //TH1F *hMOMthn = new TH1F("hMOMthn", "Angle from z of Negative muon momentum",100,-1,1);
  //TH1F *hMOMno = new TH1F("hMOMno", "Number of Muons in event",100,0.0,8.0);
  //TH1F *jpm = new TH1F("jpm", "Dimuon Invariant Mass",100,3.00,3.20); //make a new histogram too
  //TH1F *jpmomz_hist = new TH1F("jpmomz","Dimuon z Momentum",100,-5.20,5.20);
  //TH1F *jpmomt_hist = new TH1F("jpmomt","Dimuon Transverse Momentum",100,0,10.40);
  //TH1F *jpcosth_hist = new TH1F("jpcosth", "Cosine angle of J/Psi z momentum",100,-1,1);
  //TH1F *pthdecay_hist = new TH1F("pthdecay","Cosine of angle between J/Psi and Muon",100,-1,1);
  
  TH1F *BMass_hist = new TH1F("BMass", "Meson Mass",100,5,5.5); //make a new histogram too

  //now retrieve the input data and start looping over entries:
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Int_t nrunnew = 1;
   Int_t nrunold = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (nrunnew != nrunold) {  //this output helps you be sure your code is running over long intervals!
	std::cout<<"run number ="<<GLB_nrun<<"  entry number ="<<jentry<<std::endl; 
	nrunold = nrunnew;
      }
      nrunnew = Int_t(jentry/100000);

     

      if (MU_nMuon<2.0) continue; //do not bother if fewer than two muons in the event

      //hMOMno->Fill(MU_nMuon);

      Float_t primarysvxx = VERTEX_vtx_svx_beamx+ZVTX_zvtx_zPos[0]*VERTEX_vtx_svx_slopex; //find primary vertex
      Float_t primarysvxy = VERTEX_vtx_svx_beamy+ZVTX_zvtx_zPos[0]*VERTEX_vtx_svx_slopey;  
      outvtx.primx = primarysvxx;  //store the primary vtx coordinates in leaves
      outvtx.primy = primarysvxy;

      for (Int_t nMup=0;nMup<MU_nMuon;nMup++){ //muon+ loop
        Int_t trkp = MU_mu_trkIndex[nMup];    //find track corresponding to the muon
	if (trkp>299 || trkp<0) continue;    //eliminating unphysical track numbers
        Int_t mtypep = 0;  //classify muons by detector
        if (MU_mu_fidBMU[nMup]) continue;
        if (MU_mu_fidCMX[nMup]) mtypep = 2;
        if (MU_mu_fidCMU[nMup]) mtypep = 3;
        if (MU_mu_fidCMP[nMup]) mtypep = 5;
        if (MU_mu_fidCMU[nMup]&& MU_mu_fidCMX[nMup]) mtypep = 4;
        if (MU_mu_fidCMU[nMup]&& MU_mu_fidCMP[nMup]) mtypep = 6;
        if (mtypep == 0) continue;
        if (mtypep < 5) continue;  //use only muons that made it to the most shielded detector
        if (TRACK_trk_nonbc_nsvxHits[trkp] < 1) continue;  //SVX tracks only
	
       
        if (TRACK_trk_nonbc_curv[trkp]<0) continue; //Positive muon, uncomment to count all muons

        double costh=-1.02;  //Invalid in case p=0
        double p = TMath::Sqrt(MOM_pt[trkp]*MOM_pt[trkp] + MOM_pz[trkp]*MOM_pz[trkp]);
        if (p > 0.) costh=MOM_pz[trkp]/p;
        
	
	//if (TRACK_trk_nonbc_curv[trkp]<0) hMOMthp->Fill(costh); //positive muon

        //if (TRACK_trk_nonbc_curv[trkp]>0) hMOMthn->Fill(costh); //negative muon
	
	
        
	//---
        //--- This is the place to insert histograms on single-muon quantities
        //---
        //hMOMpt->Fill(MOM_pt[trkp]);  //Plot the muon transverse momentum
        //hMOMpz->Fill(MOM_pz[trkp]);  //Plot the muon transverse momentum
	//hMOMth->Fill(costh);




        for (Int_t nMum=0;nMum<MU_nMuon;nMum++){ //muon- loop
          Int_t trkm = MU_mu_trkIndex[nMum];
	  if (trkm>299 || trkm<0) continue;
          if (trkp == trkm) continue;//probably not needed
          Int_t mtypem = 0;
          if (MU_mu_fidBMU[nMum]) continue;
          if (MU_mu_fidCMX[nMum]) mtypem = 2;
          if (MU_mu_fidCMU[nMum]) mtypem = 3;
          if (MU_mu_fidCMP[nMum]) mtypem = 5;
          if (MU_mu_fidCMU[nMum]&& MU_mu_fidCMX[nMum]) mtypem = 4;
          if (MU_mu_fidCMU[nMum]&& MU_mu_fidCMP[nMum]) mtypem = 6;
          if (mtypem == 0) continue;
          if (mtypem <5) continue;
          if (TRACK_trk_nonbc_nsvxHits[trkm] < 1) continue;         
	  if (TRACK_trk_nonbc_curv[trkm]>0) continue; //negative muon

          //---
          //--- This is the place to insert histograms on paired muon quantities
          //--- and calculate variables for the output tree
          //---



	  double p_m = TMath::Sqrt(MOM_pt[trkm]*MOM_pt[trkm] + MOM_pz[trkm]*MOM_pz[trkm]);
	  double p_p = TMath::Sqrt(MOM_pt[trkp]*MOM_pt[trkp] + MOM_pz[trkp]*MOM_pz[trkp]);
          //Float_t psimass=0;
	  Float_t Mu_term_m =  2*MU_Mass_Gev*MU_Mass_Gev;
	  Float_t Mu_term_p = MOM_px[trkp]*MOM_px[trkm]  + MOM_py[trkp]*MOM_py[trkm]+ MOM_pz[trkp]*MOM_pz[trkm];
          Float_t psimass2 = Mu_term_m + 2*(TMath::Sqrt((MU_Mass_Gev*MU_Mass_Gev + p_p*p_p)*(MU_Mass_Gev*MU_Mass_Gev + p_m*p_m)) - Mu_term_p);//calculate mass of jpsi candidate
	  Float_t psimass = TMath::Sqrt(psimass2);


	  if (psimass>3.2 || psimass<3.0) continue;  //loose jpsi mass cut
	  out.jpmass=psimass;  //store calculated values in leaves

	  
	  Float_t ppz = MOM_pz[trkm] + MOM_pz[trkp];
	  outvtx.jpmomz = ppz;
	  Float_t ppx = MOM_px[trkm] + MOM_px[trkp];
	  Float_t ppy = MOM_py[trkm] + MOM_py[trkp];
	  Float_t ppt = TMath::Sqrt(ppx*ppx + ppy*ppy);
	  Float_t ppp = TMath::Sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
	  outvtx.jpmomt = ppt;

	  double costh = ppz/TMath::Sqrt((ppz*ppz + ppt*ppt));
	  outvtx.jpcosth = costh;

	  double pangle = (MOM_pz[trkp]*ppz + MOM_px[trkp]*ppx + MOM_py[trkp]*ppy) / (p_p * ppp);

	  outvtx.xintersect=0;   //change these to calculate muon track intersection
	  outvtx.yintersect=0;
	  outvtx.pthdecay = pangle;

	  if(MOM_pt[trkm] <0.8 || MOM_pt[trkp] <0.8) continue; //limiting lower bound trans mom


	  if(trkp > 75||trkm > 75)
	    {
	    std::cout<<"Value Error"<<std::endl;
	    return;
	    }
	  LsqFit l;
	  l.push(TRACK_trk_nonbc_phi[trkp] , TRACK_trk_nonbc_d0[trkp] ,TMath::Sqrt(TRKDET_trkdet_nonbc_sigD02[trkp]),MOM_pt[trkp]);

	  l.push(TRACK_trk_nonbc_phi[trkm], TRACK_trk_nonbc_d0[trkm],TMath::Sqrt(TRKDET_trkdet_nonbc_sigD02[trkm]),MOM_pt[trkm]);

	  l.fit();

	  Part_Struct.vx = l.x()- outvtx.primx;
	  Part_Struct.vy = l.y() - outvtx.primy;
	 
	  

	  for (Int_t nKa=0;nKa<TRACK_ntrk;nKa++){ //kaon loop
	    if(nKa == trkp)
	      continue;
               if(nKa == trkm)
	      continue;
          if (TRACK_trk_nonbc_nsvxHits[nKa] < 1) continue;  //SVX tracks only
	   if (TRACK_trk_nonbc_curv[nKa]<0) continue; //Positive Kaons

	  Float_t momt_1 = TMath::Sqrt(MOM_px[nKa]*MOM_px[nKa] + MOM_py[nKa]*MOM_py[nKa]);

	  
	  Float_t mom_1 = TMath::Sqrt(momt_1*momt_1 + MOM_pz[nKa]*MOM_pz[nKa]);

	  
	  //Float_t Bmtr = TMath::Sqrt((ppx+MOM_px[nKa])*(ppx+MOM_px[nKa]) + (ppy+MOM_py[nKa])*(ppy+MOM_py[nKa]));
	  //if (Bmtr < 6) continue;

	  for (Int_t mKa=0;mKa<TRACK_ntrk;mKa++){ //kaon loop
	    if(mKa == trkp)
	      continue;
             if(mKa == trkm)
	      continue;
                if(mKa == nKa)
	      continue;
          if (TRACK_trk_nonbc_nsvxHits[mKa] < 1) continue;  //SVX tracks only
	   if (TRACK_trk_nonbc_curv[mKa]>0) continue; //Negative kaons
	  Float_t momt_2 = TMath::Sqrt(MOM_px[mKa]*MOM_px[mKa] + MOM_py[mKa]*MOM_py[mKa]);
	  Float_t mom_2 = TMath::Sqrt(momt_2*momt_2 + MOM_pz[mKa]*MOM_pz[mKa]);
	  Float_t dot_Ka = MOM_px[mKa]*MOM_px[nKa] +MOM_py[mKa]*MOM_py[nKa] + MOM_pz[mKa]*MOM_pz[nKa];
	  Float_t Phimass2 = 2*K_Mass*K_Mass + 2*TMath::Sqrt((K_Mass*K_Mass + mom_1*mom_1)*(K_Mass*K_Mass + mom_2*mom_2)) - 2*dot_Ka;
	  Float_t Phimass = TMath::Sqrt(Phimass);
          if(0>Phimass || 1.15<Phimass) continue;
	  Float_t PhiMomx = MOM_px[mKa]+MOM_px[nKa];
          Float_t PhiMomy = MOM_py[mKa]+MOM_py[nKa];
          Float_t PhiMomz = MOM_pz[mKa]+MOM_pz[nKa];
          Float_t PhiMOM = TMath::Sqrt(PhiMomx*PhiMomx +PhiMomy*PhiMomy +PhiMomz*PhiMomz);
          Float_t dot_Phi = PhiMomx*ppx +PhiMomy*ppy + PhiMomz*ppz;
          //Float_t BSmass = Phimass*Phimass + psimass*psimass + 2*TMath::Sqrt((Phimass*Phimass + PhiMOM*PhiMOM)*(psimass*psimass + ppp*ppp)) - 2*dot_Phi;
	  Float_t BS_Ent1 =  2*TMath::Sqrt((K_mass*K_mass + mom_1*mom_1)*(Mu_Mass_Gev*Mu_Mass_Gev + p_m*p_m)) + 2*TMath::Sqrt((K_mass*K_mass + mom_1*mom_1)*(Mu_Mass_Gev*Mu_Mass_Gev + p_p*p_p));
          Float_t BS_Ent2 =  2*TMath::Sqrt((K_mass*K_mass + mom_2*mom_2)*(Mu_Mass_Gev*Mu_Mass_Gev + p_m*p_m)) + 2*TMath::Sqrt((K_mass*K_mass + mom_2*mom_2)*(Mu_Mass_Gev*Mu_Mass_Gev + p_p*p_p));
	  Float_t BS_dot1 = (MOM_px[mKa]*MOM_px[trkm]  + MOM_py[mKa]*MOM_py[trkm]+ MOM_pz[mKa]*MOM_pz[trkm]) +  (MOM_px[mKa]*MOM_px[trkp]  + MOM_py[mKa]*MOM_py[trkp]+ MOM_pz[mKa]*MOM_pz[trkp]);
            Float_t BS_dot2 = (MOM_px[nKa]*MOM_px[trkm]  + MOM_py[nKa]*MOM_py[trkm]+ MOM_pz[nKa]*MOM_pz[trkm]) +  (MOM_px[nKa]*MOM_px[trkp]  + MOM_py[nKa]*MOM_py[trkp]+ MOM_pz[nKa]*MOM_pz[trkp]);
	    Float_t BSmass2 = Phimass2 + psimass2 + BS_Ent1 + BS_Ent2 - 2*BS_dot1 - 2*BS_dot1;
	  	  Float_t BSmass = TMath::Sqrt(BSmass2);
	  

	  if(BSmass>200)
	    continue;
	  
	  Part_Struct.phimass = Phimass;
	  Part_Struct.B_mass = BSmass;
	  BMass_hist->Fill(BSmass);

	  
	  //Jpcosth_hist->Fill(costh);
	  //jpmomz_hist->Fill(ppz);
	  //jpmomt_hist->Fill(ppt);
	  //pthdecay_hist->Fill(pangle);
	  //jpm->Fill(psimass); //fill the histogram with mass values
	  output->Fill();  //fill the tree with the calculated leaf values
         	} //K2 loop
         }//K1 loop
        } //muon- loop
      }// muon+ loop
   }
   f2->cd();
   output->Write();  //once the tree has been fully filled, write to file

   //---
   //--- Insert Write() statements for each of the histograms you add here
   //---
   //jpm->Write();     //write mass plot histogram
   //jpmomz_hist->Write();
   //jpmomt_hist->Write();
   //jpcosth_hist->Write();
   //pthdecay_hist->Write();
   //hMOMpt->Write(); //write pt histogram
   //hMOMpz->Write(); //write pz histogram
   //hMOMth->Write(); //write theta histogram
   //hMOMthp->Write(); //write theta histogram
   //hMOMthn->Write(); //write theta histogram
   //hMOMno->Write(); //write theta histogram
   BMass_hist->Write();
   f2->Close();      //close the file
   std::cout<<"loop finished"<<std::endl;
}

