#define ana_mc_clas12_cxx
// The class definition in ana_mc_clas12.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("ana_mc_clas12.C")
// root> T->Process("ana_mc_clas12.C","some options")
// root> T->Process("ana_mc_clas12.C+")
//

#include "ana_mc_clas12.h"
// #include <TH2.h>
// #include <TStyle.h>
// #include <TH3.h>
// #include <TTree.h>
// #include <TSelector.h>
// #include <TString.h>
// #include <fstream>    // for fstream
// #include <stdio.h>
// #include <string.h>
// #include <math.h>
// #include <iostream>   // for cout, endl
// #include <map>
// #include <string>
// #include <TTreeReader.h>
// #include <TTreeReaderValue.h>
// #include <TTreeReaderArray.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "Rtypes.h"
using namespace std;

void ana_mc_clas12::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString Option = GetOption();
   outfile = new TFile("ana_dvcs_MC.root","RECREATE");

   AddBranches();
   AddMCBranches();

   firstCall = kTRUE;

   Pmass = 0.938272;
   Nmass = 0.93957;
   Dmass = 1.8756;
   Elmass = 0.00051;
   LightSpeed = 29.9792458;
   Ebeam = 10.6;
   
   tstart = 124.25;

   ElectronBeam.SetXYZT(0,0,Ebeam,Ebeam);
   Target_Vec.SetXYZT(0,0,0,Dmass); 
   PTarget_Vec.SetXYZT(0,0,0,Pmass); 
   NTarget_Vec.SetXYZT(0,0,0,Nmass); 
   //   DefineHistos();
   firstCall = kFALSE;
}

void ana_mc_clas12::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  //   TString option = GetOption();
 }

Bool_t ana_mc_clas12::Process(Long64_t entry)
{

   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ana_mc_clas12::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  fReader.SetLocalEntry(entry);         
  nelec=0;
  i_el = -1;
  nprot=0;
  i_pr = -1;
  nneut=0;
  i_n=-1;
  nphotEC=0;
  i_phEC = -1;
  nphotFT=0;
  i_phFT = -1;
  nphot=0;
  i_ph=-1;
  npion=0;
  FD_proton=0;
  CD_proton=0;
  FD_neutron=0;
  CD_neutron=0;
  ncharged=0;
  Float_t pmom_photon;
  pmom_photon=1.;
  Float_t mean_p_neutron;
  mean_p_neutron=0.67;
  Float_t delta_p_neutron;
  delta_p_neutron=2.4;
  TVector3 cluster;
  Float_t delta_theta_cluster;
  delta_theta_cluster=10;
  Float_t delta_phi_cluster;
  delta_phi_cluster=10;
  Float_t theta_cluster;
  theta_cluster=0;  
  Float_t phi_cluster;
  phi_cluster=0;
  veto=-1;
  nneut_matched=0;
  Calc_kine_MC(); 
  MC_tree->Fill();
  
  int size = REC_Particle_pid.GetSize();
  //  int size_event = REC_Event_NRUN.GetSize();
  int size_cluster = __CND_Clusters_veto.GetSize();

  for(int j=0;j<size;j++)
    {
      Float_t pmom[size];
      pmom[j]=TMath::Sqrt((REC_Particle_px)[j]*(REC_Particle_px)[j]+(REC_Particle_py)[j]*(REC_Particle_py)[j]+(REC_Particle_pz)[j]*(REC_Particle_pz)[j]);

      //electron ID
      if((REC_Particle_charge)[j]<0 && (REC_Particle_pid)[j]==11 && pmom[j]>1.)
	{
	  nelec++;
	  i_el=j;
	  El_Vec.SetPxPyPzE((REC_Particle_px)[j],(REC_Particle_py)[j],(REC_Particle_pz)[j],TMath::Sqrt(pmom[j]*pmom[j]+Elmass*Elmass));
	  Q2 = 4*Ebeam*El_Vec.P()*TMath::Power(TMath::Sin(El_Vec.Theta()/2),2.);
	  Float_t nu = Ebeam-El_P;
	  W = TMath::Sqrt(Pmass*Pmass + 2*Pmass*nu - Q2);
	  El_P = pmom[j];
	  Xbj = Q2/(2*Pmass*(Ebeam-El_P));
	  El_Theta = El_Vec.Theta()*180./TMath::Pi();
	  El_Phi = El_Vec.Phi()*180./TMath::Pi();
	  El_vz = (REC_Particle_vz)[j];
	  El_vx = (REC_Particle_vx)[j];
	  El_vy = (REC_Particle_vy)[j];
	}
  
      if(nelec>=1)
	{
	  if(j!=i_el && (REC_Particle_charge)[j]!=0)ncharged++;
	  
	  //photon ID 
	  if(j!=i_el && (REC_Particle_charge)[j]==0 && (REC_Particle_pid)[j]==22)
	    {
	      nphot++;
	      if(pmom[j]>pmom_photon) {
	      i_ph=j;
	      pmom_photon=pmom[j];
	      Ph_Vec.SetPxPyPzE((REC_Particle_px)[j],(REC_Particle_py)[j],(REC_Particle_pz)[j],pmom[j]);
	      Ph_P = Ph_Vec.E();
	      Ph_Theta = Ph_Vec.Theta()*180./TMath::Pi();
	      Ph_Phi = Ph_Vec.Phi()*180./TMath::Pi();  
	      }}	  
	  if(j!=i_el && j!=i_ph && (REC_Particle_pid)[j]==2212 && (REC_Particle_charge)[j]>0)
	    {
	      //proton ID
	      nprot++;
	      i_pr=j;
	      Pr_Vec.SetPxPyPzE((REC_Particle_px)[j],(REC_Particle_py)[j],(REC_Particle_pz)[j],TMath::Sqrt(pmom[j]*pmom[j]+Pmass*Pmass));
	      Pr_P = pmom[j];
	      Pr_Theta = Pr_Vec.Theta()*180./TMath::Pi();
	      Pr_Phi = Pr_Vec.Phi()*180./TMath::Pi();
	      Pr_vz = (REC_Particle_vz)[j];
	      if((REC_Particle_status)[j]>=2000 && (REC_Particle_status)[j]<4000)FD_proton=1;
	      else if((REC_Particle_status)[j]>=4000)CD_proton=1;
	    }
	  if(j!=i_el && j!=i_ph && TMath::Abs((REC_Particle_pid)[j])==211)
	    {
	      npion++;
	    }

	  if(j!=i_el && j!=i_ph && j!=i_pr && (REC_Particle_charge)[j]==0 && (REC_Particle_pid)[j]==2112 && pmom[j]>0)
	    {
              nneut++;
	      if(pmom[j]>0.2 && pmom[j]<3 && TMath::Abs(pmom[j]-mean_p_neutron)<delta_p_neutron) {		
	      i_n=j;
	      delta_p_neutron=TMath::Abs(pmom[j]-mean_p_neutron);
	      N_Vec.SetPxPyPzE((REC_Particle_px)[j],(REC_Particle_py)[j],(REC_Particle_pz)[j],TMath::Sqrt(pmom[j]*pmom[j]+Nmass*Nmass));
	      N_P = pmom[j];
	      N_Theta = N_Vec.Theta()*180./TMath::Pi();
	      N_Phi = N_Vec.Phi()*180./TMath::Pi();
	      if((REC_Particle_status)[j]>=2000 && (REC_Particle_status)[j]<4000)FD_neutron=1;
	      else if((REC_Particle_status)[j]>=4000)
		{
		  CD_neutron=1;
		  if(size_cluster!=0)
		    {
		      for(int k=0;k<size_cluster;k++)
			{
			  cluster.SetXYZ((__CND_Clusters_x)[k],(__CND_Clusters_y)[k],(__CND_Clusters_z)[k]);
			  theta_cluster=cluster.Theta()*180./TMath::Pi();
			  phi_cluster=cluster.Phi()*180./TMath::Pi();
			  if(TMath::Abs(theta_cluster-N_Theta)<=delta_theta_cluster && TMath::Abs(phi_cluster-N_Phi)<=delta_phi_cluster)
			    {
			      delta_theta_cluster=TMath::Abs(N_Theta-theta_cluster);
			      delta_phi_cluster=TMath::Abs(N_Phi-phi_cluster);
			      veto=(__CND_Clusters_veto)[k];
			    }
			}
		    }
		}	      
	      }
	    }
	}
    }
  if(nelec>=1)el_tree->Fill();
  if(nelec>=1 && nprot>=1 && nphot>=1 && i_ph>=0)
    {
      Calc_kine_P();
      pDVCS_tree->Fill();
    }
  if(nelec>=1 && nneut>=1 && nphot>=1 && i_ph>=0 && i_n>=0)
    {
      Calc_kine_N();
      nDVCS_tree->Fill();
    }
  return kTRUE;
}

void ana_mc_clas12::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void ana_mc_clas12::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  outfile->cd();
  el_tree->Write();
  pDVCS_tree->Write();
  nDVCS_tree->Write();
  MC_tree->Write();
  outfile->Write();
  outfile->Close();

}

void ana_mc_clas12::Calc_kine_P()
{
    TVector3 VelectronIn,VelectronOut,VprotonOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
    TVector3 VX;

    VelectronIn    = ElectronBeam.Vect();
    VelectronOut   = El_Vec.Vect();
    VprotonOut     = Pr_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-El_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VprotonOut.Cross(Vvirtualphoton);
    VhadroPP       = VprotonOut.Cross(VphotonOut);

    Phi_Pr = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VprotonOut)>0.)  Phi_Pr = 360.-Phi_Pr;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_Pr  = (Pr_Vec-PTarget_Vec).M2();
    t_Ph  = (ElectronBeam-El_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+PTarget_Vec-Ph_Vec-El_Vec-Pr_Vec;
    mm2_epg = (ElectronBeam+Target_Vec-Pr_Vec-El_Vec-Ph_Vec).M2();
    mm2_epg_p = (ElectronBeam+PTarget_Vec-Pr_Vec-El_Vec-Ph_Vec).M2();
    
    Xbal  = BalV.X();
    Ybal  = BalV.Y();
    Zbal  = BalV.Z();
    Ebal  = BalV.E();
}
void ana_mc_clas12::Calc_kine_N()
{
  
    TVector3 VelectronIn,VelectronOut,VneutronOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
    TVector3 VX;

    VelectronIn    = ElectronBeam.Vect();
    VelectronOut   = El_Vec.Vect();
    VneutronOut     = N_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-El_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VneutronOut.Cross(Vvirtualphoton);
    VhadroPP       = VneutronOut.Cross(VphotonOut);

    Phi_N = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VneutronOut)>0.)  Phi_N = 360.-Phi_N;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_N  = (N_Vec-NTarget_Vec).M2();
    t_Ph  = (ElectronBeam-El_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+Target_Vec-Ph_Vec-El_Vec-N_Vec;

    TLorentzVector BalV_n = ElectronBeam+NTarget_Vec-Ph_Vec-El_Vec-N_Vec;
    
    TVector3 missing_gamma_n = (ElectronBeam+NTarget_Vec-El_Vec-N_Vec).Vect();
    
    mm2_eng = (ElectronBeam+Target_Vec-N_Vec-El_Vec-Ph_Vec).M2();
    mm2_eng_n = (ElectronBeam+NTarget_Vec-N_Vec-El_Vec-Ph_Vec).M2();

    Xbal  = BalV.X();
    Ybal  = BalV.Y();
    Zbal  = BalV.Z();
    Ebal  = BalV.E();

    Xbal_n  = BalV_n.X();
    Ybal_n  = BalV_n.Y();
    Zbal_n  = BalV_n.Z();
    Ebal_n  = BalV_n.E();
    
    theta_gamma_X=180./TMath::Pi()*TMath::ACos((missing_gamma_n.Dot(VphotonOut))/(missing_gamma_n.Mag()*VphotonOut.Mag()));
    theta_gamma_e=180./TMath::Pi()*TMath::ACos((VphotonOut.Dot(VelectronOut))/(VphotonOut.Mag()*VelectronOut.Mag()));
    theta_n_e=180./TMath::Pi()*TMath::ACos((VneutronOut.Dot(VelectronOut))/(VneutronOut.Mag()*VelectronOut.Mag()));
}

void ana_mc_clas12::Calc_kine_MC()
{
  
   TVector3 VelectronIn,VelectronOut,VprotonOut,VneutronOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
     TVector3 VX;

     p_spec=0;
     n_spec=0;
     // Phi_MC          = -900.;          // initialization - impossible value
     // t_Pr_MC = -900.;
     // t_N_MC = -900.;
     VelectronIn.SetXYZ(0.,0.,Ebeam);
     // VneutronOut.SetXYZ(-500.,-500.,-500.);
     // VprotonOut.SetXYZ(-500.,-500.,-500.);
     // Pr_Vec_MC.SetPxPyPzE(-500.,-500.,-500.,-500.);
     // N_Vec_MC.SetPxPyPzE(-500.,-500.,-500.,-500.);
     if (MC_Particle_px.IsEmpty() || MC_Particle_py.IsEmpty() || MC_Particle_pz.IsEmpty() || MC_Particle_pid.IsEmpty())
     {
       cout<<"Skipped because of an event with empty: "
 	  << (MC_Particle_px.IsEmpty() ? "MC_Particle_px " : " ")
 	  << (MC_Particle_py.IsEmpty() ? "MC_Particle_py " : " ")
 	  << (MC_Particle_pz.IsEmpty() ? "MC_Particle_pz " : " ")
 	  << (MC_Particle_pid.IsEmpty() ? "MC_Particle_pid " : " ")
 	  << std::endl;
       return;
     }
    
     VelectronOut.SetXYZ((MC_Particle_px)[0],(MC_Particle_py)[0],(MC_Particle_pz)[0]);

     if((MC_Particle_pid)[1]==2212.)
       {
 	p_spec=1;
	VneutronOut.SetXYZ((MC_Particle_px)[2],(MC_Particle_py)[2],(MC_Particle_pz)[2]);
	VprotonOut.SetXYZ((MC_Particle_px)[1],(MC_Particle_py)[1],(MC_Particle_pz)[1]);
      }
    if((MC_Particle_pid)[1]==2112.)
      {
	n_spec=1;
	VprotonOut.SetXYZ((MC_Particle_px)[2],(MC_Particle_py)[2],(MC_Particle_pz)[2]);
	VneutronOut.SetXYZ((MC_Particle_px)[1],(MC_Particle_py)[1],(MC_Particle_pz)[1]);
      }
    VphotonOut.SetXYZ((MC_Particle_px)[3],(MC_Particle_py)[3],(MC_Particle_pz)[3]);
    Vvirtualphoton = (ElectronBeam-El_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    if((MC_Particle_pid)[1]==2112){
      Vhadro         = VprotonOut.Cross(Vvirtualphoton);
      VhadroPP       = VprotonOut.Cross(VphotonOut);
    }

    if((MC_Particle_pid)[1]==2212){
      Vhadro         = VneutronOut.Cross(Vvirtualphoton);
      VhadroPP       = VneutronOut.Cross(VphotonOut);
    }
    Phi_MC = 180./TMath::Pi()*Vlepto.Angle(Vhadro);

    if((MC_Particle_pid)[1]==2112 && Vlepto.Dot(VprotonOut)>0.)  Phi_MC = 360.-Phi_MC;
    if((MC_Particle_pid)[1]==2212 && Vlepto.Dot(VneutronOut)>0.)  Phi_MC = 360.-Phi_MC;
    VX = (ElectronBeam + Target_Vec - Pr_Vec - El_Vec).Vect();

    Float_t pmom[4];
    for(Int_t i=0;i<4;i++)
      {
	pmom[i]=TMath::Sqrt((MC_Particle_px)[i]*(MC_Particle_px)[i]+(MC_Particle_py)[i]*(MC_Particle_py)[i]+(MC_Particle_pz)[i]*(MC_Particle_pz)[i]);
      }

    El_Vec_MC.SetPxPyPzE((MC_Particle_px)[0],(MC_Particle_py)[0],(MC_Particle_pz)[0],pmom[0]);
    El_P_MC=El_Vec_MC.P();
    El_Theta_MC=180./TMath::Pi()*El_Vec_MC.Theta();
    El_Phi_MC=180./TMath::Pi()*El_Vec_MC.Phi();

    if((MC_Particle_pid)[1]==2112.){
      Pr_Vec_MC.SetPxPyPzE((MC_Particle_px)[2],(MC_Particle_py)[2],(MC_Particle_pz)[2],TMath::Sqrt(pmom[2]*pmom[2]+Pmass*Pmass));
      Pr_P_MC=Pr_Vec_MC.P();
      Pr_Theta_MC=180./TMath::Pi()*Pr_Vec_MC.Theta();
      Pr_Phi_MC=180./TMath::Pi()*Pr_Vec_MC.Phi();
      N_Vec_MC.SetPxPyPzE((MC_Particle_px)[1],(MC_Particle_py)[1],(MC_Particle_pz)[1],TMath::Sqrt(pmom[1]*pmom[1]+Nmass*Nmass));
      N_P_MC=N_Vec_MC.P();
      N_Theta_MC=180./TMath::Pi()*N_Vec_MC.Theta();
      N_Phi_MC=180./TMath::Pi()*N_Vec_MC.Phi();
    }

    if((MC_Particle_pid)[1]==2212.){
      N_Vec_MC.SetPxPyPzE((MC_Particle_px)[2],(MC_Particle_py)[2],(MC_Particle_pz)[2],TMath::Sqrt(pmom[2]*pmom[2]+Nmass*Nmass));
      N_P_MC=N_Vec_MC.P();
      N_Theta_MC=180./TMath::Pi()*N_Vec_MC.Theta();
      N_Phi_MC=180./TMath::Pi()*N_Vec_MC.Phi();
      Pr_Vec_MC.SetPxPyPzE((MC_Particle_px)[1],(MC_Particle_py)[1],(MC_Particle_pz)[1],TMath::Sqrt(pmom[1]*pmom[1]+Pmass*Pmass));
      Pr_P_MC=Pr_Vec_MC.P();
      Pr_Theta_MC=180./TMath::Pi()*Pr_Vec_MC.Theta();
      Pr_Phi_MC=180./TMath::Pi()*Pr_Vec_MC.Phi();
    }

    // Ph_Vec_MC.SetPxPyPzE((MC_Particle_px)[2],(MC_Particle_py)[2],(MC_Particle_pz)[2],pmom[2]);
    // Ph_P_MC=Ph_Vec_MC.P();
    // Ph_Theta_MC=180./TMath::Pi()*Ph_Vec_MC.Theta();
    // Ph_Phi_MC=180./TMath::Pi()*Ph_Vec_MC.Phi();
    Ph_Vec_MC.SetPxPyPzE((MC_Particle_px)[3],(MC_Particle_py)[3],(MC_Particle_pz)[3],pmom[3]);
    Ph_P_MC=Ph_Vec_MC.P();
    Ph_Theta_MC=180./TMath::Pi()*Ph_Vec_MC.Theta();
    Ph_Phi_MC=180./TMath::Pi()*Ph_Vec_MC.Phi();

    Float_t nu = Ebeam-El_P_MC;
    Q2_MC = 4*Ebeam*El_Vec_MC.P()*TMath::Power(TMath::Sin(El_Vec_MC.Theta()/2),2.);
    Xbj_MC = Q2_MC/(2*Pmass*(Ebeam-El_P_MC));
    if((MC_Particle_pid)[1]==2112.)W_MC = TMath::Sqrt(Pmass*Pmass + 2*Pmass*nu - Q2_MC);
    if((MC_Particle_pid)[1]==2212.)W_MC = TMath::Sqrt(Nmass*Nmass + 2*Nmass*nu - Q2_MC);
    if((MC_Particle_pid)[1]==2112.)t_Pr_MC  = (Pr_Vec_MC-PTarget_Vec).M2();
    if((MC_Particle_pid)[1]==2212.)t_N_MC = (N_Vec_MC-NTarget_Vec).M2();
    t_Ph_MC  = (ElectronBeam-El_Vec_MC-Ph_Vec_MC).M2();

    theta_n_e_MC=180./TMath::Pi()*TMath::ACos(((N_Vec_MC.Vect()).Dot((El_Vec_MC.Vect()))/((N_Vec_MC.Vect()).Mag()*(El_Vec_MC.Vect()).Mag())));
    theta_gamma_e_MC=180./TMath::Pi()*TMath::ACos(((Ph_Vec_MC.Vect()).Dot((El_Vec_MC.Vect()))/((Ph_Vec_MC.Vect()).Mag()*(El_Vec_MC.Vect()).Mag())));
}
void ana_mc_clas12::AddBranches()
{
   el_tree = (TTree*) new TTree("el","el");
   el_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   el_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   //   el_tree->Branch("bpola",&bpola,"bpola/F");
   el_tree->Branch("helicity",&helicity,"helicity/I");
   el_tree->Branch("TarPol",&TarPol,"TarPol/F");
   el_tree->Branch("Q2",&Q2,"Q2/F");
   el_tree->Branch("Xbj", &Xbj, "Xbj/F");
   el_tree->Branch("W", &W, "W/F");
   el_tree->Branch("El_P",&El_P,"El_P/F");
   el_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   el_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   el_tree->Branch("El_vx",&El_vx,"El_vx/F");
   el_tree->Branch("El_vy",&El_vy,"El_vy/F");
   el_tree->Branch("El_vz",&El_vz,"El_vz/F");
   el_tree->Branch("El_vx_MC",&El_vx_MC,"El_vx_MC/F");
   el_tree->Branch("El_vy_MC",&El_vy_MC,"El_vy_MC/F");
   el_tree->Branch("El_vz_MC",&El_vz_MC,"El_vz_MC/F");
   el_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   el_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   el_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   el_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   el_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   el_tree->Branch("W_MC", &W_MC, "W_MC/F");
   
   pDVCS_tree = (TTree*) new TTree("pDVCS","pDVCS");
   pDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   pDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   pDVCS_tree->Branch("helicity",&helicity,"helicity/I");
   pDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
   pDVCS_tree->Branch("t_Pr",&t_Pr,"t_Pr/F");
   pDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
   pDVCS_tree->Branch("Phi_Pr",&Phi_Pr,"Phi_Pr/F");
   pDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
   pDVCS_tree->Branch("Q2",&Q2,"Q2/F");
   pDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
   pDVCS_tree->Branch("W", &W, "W/F");
   pDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
   pDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
   pDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
   pDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
   pDVCS_tree->Branch("nelec",&nelec,"nelec/I");
   pDVCS_tree->Branch("El_P",&El_P,"El_P/F");
   pDVCS_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   pDVCS_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   pDVCS_tree->Branch("El_vz",&El_vz,"El_vz/F");
   pDVCS_tree->Branch("El_vx",&El_vx,"El_vx/F");
   pDVCS_tree->Branch("El_vy",&El_vy,"El_vy/F");
   pDVCS_tree->Branch("nprot",&nprot,"nprot/I");
   pDVCS_tree->Branch("Pr_P",&Pr_P,"Pr_P/F");
   pDVCS_tree->Branch("Pr_Theta",&Pr_Theta,"Pr_Theta/F");
   pDVCS_tree->Branch("Pr_Phi",&Pr_Phi,"Pr_Phi/F");
   pDVCS_tree->Branch("Pr_vz",&Pr_vz,"Pr_vz/F");
   pDVCS_tree->Branch("FD_proton",&FD_proton,"FD_proton/I");
   pDVCS_tree->Branch("CD_proton",&CD_proton,"CD_proton/I");
   pDVCS_tree->Branch("nphotEC",&nphotEC,"nphotEC/I");
   pDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
   pDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
   pDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
   pDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
   pDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
   pDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
   pDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
   pDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
   pDVCS_tree->Branch("nphot",&nphot,"nphot/I");
   pDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
   pDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
   pDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
   pDVCS_tree->Branch("mm2_epg",&mm2_epg,"mm2_epg/F");
   pDVCS_tree->Branch("mm2_epg_p",&mm2_epg_p,"mm2_epg_p/F");
   pDVCS_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   pDVCS_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   pDVCS_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   pDVCS_tree->Branch("Pr_P_MC",&Pr_P_MC,"Pr_P_MC/F");
   pDVCS_tree->Branch("Pr_Theta_MC",&Pr_Theta_MC,"Pr_Theta_MC/F");
   pDVCS_tree->Branch("Pr_Phi_MC",&Pr_Phi_MC,"Pr_Phi_MC/F");
   pDVCS_tree->Branch("Ph_P_MC",&Ph_P_MC,"Ph_P_MC/F");
   pDVCS_tree->Branch("Ph_Theta_MC",&Ph_Theta_MC,"Ph_Theta_MC/F");
   pDVCS_tree->Branch("Ph_Phi_MC",&Ph_Phi_MC,"Ph_Phi_MC/F");
   pDVCS_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   pDVCS_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   pDVCS_tree->Branch("W_MC", &W_MC, "W_MC/F");
   pDVCS_tree->Branch("t_Pr_MC",&t_Pr_MC,"t_Pr_MC/F");
   pDVCS_tree->Branch("t_N_MC",&t_N_MC,"t_N_MC/F");
   pDVCS_tree->Branch("t_Ph_MC",&t_Ph_MC,"t_Ph_MC/F");
   pDVCS_tree->Branch("Phi_MC",&Phi_MC,"Phi_MC/F");
   pDVCS_tree->Branch("n_spec",&n_spec,"n_spec/I");
   pDVCS_tree->Branch("p_spec",&p_spec,"p_spec/I");

   nDVCS_tree = (TTree*) new TTree("nDVCS","nDVCS");
   nDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   nDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   nDVCS_tree->Branch("helicity",&helicity,"helicity/I");
   nDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
   nDVCS_tree->Branch("Q2",&Q2,"Q2/F");
   nDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
   nDVCS_tree->Branch("W", &W, "W/F");
   nDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
   nDVCS_tree->Branch("t_N",&t_N,"t_N/F");
   nDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
   nDVCS_tree->Branch("Phi_N",&Phi_N,"Phi_N/F");
   nDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
   nDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
   nDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
   nDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
   nDVCS_tree->Branch("Xbal_n",&Xbal_n,"Xbal_n/F");
   nDVCS_tree->Branch("Ybal_n",&Ybal_n,"Ybal_n/F");
   nDVCS_tree->Branch("Zbal_n",&Zbal_n,"Zbal_n/F");
   nDVCS_tree->Branch("Ebal_n",&Ebal_n,"Ebal_n/F");
   nDVCS_tree->Branch("theta_gamma_X",&theta_gamma_X,"theta_gamma_X/F");
   nDVCS_tree->Branch("theta_gamma_e",&theta_gamma_e,"theta_gamma_e/F");
   nDVCS_tree->Branch("theta_n_e",&theta_n_e,"theta_n_e/F");
   nDVCS_tree->Branch("nelec",&nelec,"nelec/I");
   nDVCS_tree->Branch("El_P",&El_P,"El_P/F");
   nDVCS_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   nDVCS_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   nDVCS_tree->Branch("El_vz",&El_vz,"El_vz/F");
   nDVCS_tree->Branch("El_vx",&El_vx,"El_vx/F");
   nDVCS_tree->Branch("El_vy",&El_vy,"El_vy/F");
   nDVCS_tree->Branch("nneut",&nneut,"nneut/I");
   nDVCS_tree->Branch("nneut_matched",&nneut_matched,"nneut_matched/I");
   nDVCS_tree->Branch("veto",&veto,"veto/I");
   nDVCS_tree->Branch("FD_neutron",&FD_neutron,"FD_neutron/I");
   nDVCS_tree->Branch("CD_neutron",&CD_neutron,"CD_neutron/I");
   nDVCS_tree->Branch("N_P",&N_P,"N_P/F");
   nDVCS_tree->Branch("N_Theta",&N_Theta,"N_Theta/F");
   nDVCS_tree->Branch("N_Phi",&N_Phi,"N_Phi/F");
   nDVCS_tree->Branch("nphotEC",&nphotEC,"nphotEC/I");
   nDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
   nDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
   nDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
   nDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
   nDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
   nDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
   nDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
   nDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
   nDVCS_tree->Branch("nphot",&nphot,"nphot/I");
   nDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
   nDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
   nDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
   nDVCS_tree->Branch("mm2_eng",&mm2_eng,"mm2_eng/F");
   nDVCS_tree->Branch("mm2_eng_n",&mm2_eng_n,"mm2_eng_n/F");
   nDVCS_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   nDVCS_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   nDVCS_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   nDVCS_tree->Branch("Pr_P_MC",&Pr_P_MC,"Pr_P_MC/F");
   nDVCS_tree->Branch("Pr_Theta_MC",&Pr_Theta_MC,"Pr_Theta_MC/F");
   nDVCS_tree->Branch("Pr_Phi_MC",&Pr_Phi_MC,"Pr_Phi_MC/F");
   nDVCS_tree->Branch("N_P_MC",&N_P_MC,"N_P_MC/F");
   nDVCS_tree->Branch("N_Theta_MC",&N_Theta_MC,"N_Theta_MC/F");
   nDVCS_tree->Branch("N_Phi_MC",&N_Phi_MC,"N_Phi_MC/F");
   nDVCS_tree->Branch("Ph_P_MC",&Ph_P_MC,"Ph_P_MC/F");
   nDVCS_tree->Branch("Ph_Theta_MC",&Ph_Theta_MC,"Ph_Theta_MC/F");
   nDVCS_tree->Branch("Ph_Phi_MC",&Ph_Phi_MC,"Ph_Phi_MC/F");
   nDVCS_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   nDVCS_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   nDVCS_tree->Branch("W_MC", &W_MC, "W_MC/F");
   nDVCS_tree->Branch("t_Pr_MC",&t_Pr_MC,"t_Pr_MC/F");
   nDVCS_tree->Branch("t_N_MC",&t_N_MC,"t_N_MC/F");
   nDVCS_tree->Branch("t_Ph_MC",&t_Ph_MC,"t_Ph_MC/F");
   nDVCS_tree->Branch("theta_n_e_MC",&theta_n_e_MC,"theta_n_e_MC/F");
   nDVCS_tree->Branch("theta_gamma_e_MC",&theta_gamma_e_MC,"theta_gamma_e_MC/F");
   nDVCS_tree->Branch("Phi_MC",&Phi_MC,"Phi_MC/F");
   nDVCS_tree->Branch("n_spec",&n_spec,"n_spec/I");
   nDVCS_tree->Branch("p_spec",&p_spec,"p_spec/I");
   nDVCS_tree->Branch("nprot",&nprot,"nprot/I");
   nDVCS_tree->Branch("npion",&npion,"npion/I");
   nDVCS_tree->Branch("ncharged",&ncharged,"ncharged/I");
}
void ana_mc_clas12::AddMCBranches()
{
   MC_tree = (TTree*) new TTree("MC","MC");
   MC_tree->Branch("El_vx_MC",&El_vx_MC,"El_vx_MC/F");
   MC_tree->Branch("El_vy_MC",&El_vy_MC,"El_vy_MC/F");
   MC_tree->Branch("El_vz_MC",&El_vz_MC,"El_vz_MC/F");
   MC_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   MC_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   MC_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   MC_tree->Branch("Pr_P_MC",&Pr_P_MC,"Pr_P_MC/F");
   MC_tree->Branch("Pr_Theta_MC",&Pr_Theta_MC,"Pr_Theta_MC/F");
   MC_tree->Branch("Pr_Phi_MC",&Pr_Phi_MC,"Pr_Phi_MC/F");
   MC_tree->Branch("N_P_MC",&N_P_MC,"N_P_MC/F");
   MC_tree->Branch("N_Theta_MC",&N_Theta_MC,"N_Theta_MC/F");
   MC_tree->Branch("N_Phi_MC",&N_Phi_MC,"N_Phi_MC/F");
   MC_tree->Branch("Ph_P_MC",&Ph_P_MC,"Ph_P_MC/F");
   MC_tree->Branch("Ph_Theta_MC",&Ph_Theta_MC,"Ph_Theta_MC/F");
   MC_tree->Branch("Ph_Phi_MC",&Ph_Phi_MC,"Ph_Phi_MC/F");
   MC_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   MC_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   MC_tree->Branch("W_MC", &W_MC, "W_MC/F");
   MC_tree->Branch("t_Pr_MC",&t_Pr_MC,"t_Pr_MC/F");
   MC_tree->Branch("t_N_MC",&t_N_MC,"t_N_MC/F");
   MC_tree->Branch("t_Ph_MC",&t_Ph_MC,"t_Ph_MC/F");
   MC_tree->Branch("theta_n_e_MC",&theta_n_e_MC,"theta_n_e_MC/F");
   MC_tree->Branch("Phi_MC",&Phi_MC,"Phi_MC/F");
   MC_tree->Branch("n_spec",&n_spec,"n_spec/I");
   MC_tree->Branch("p_spec",&p_spec,"p_spec/I");
}
