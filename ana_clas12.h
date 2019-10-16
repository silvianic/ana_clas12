//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep  6 13:41:00 2019 by ROOT version 6.14/06
// from TTree clas12/clas12
// found on file: calib_clas_006595.evio.00773.hipo.root
//////////////////////////////////////////////////////////

#ifndef ana_clas12_h
#define ana_clas12_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVirtualIndex.h>
#include <iostream>             // std::cout, std::endl
#include <fstream>              // std::ifstream
#include <sstream>   
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <cstdlib>
#include <map>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;

class ana_clas12 : public TSelector {
public :
   TTree *el_tree,*pDVCS_tree,*nDVCS_tree; 
   TTree *MC_tree; 
   TFile *outfile; 
   Bool_t   firstCall; 
   void    AddBranches();
   void    AddMCBranches();
   void    Calc_kine_P();
   void    Calc_kine_N();   
   void    Calc_kine_MC();


// Fixed size dimensions of array or collection stored in the TTree if any.
   Int_t nelec,nprot,nneut,nneut_matched,nphotFT,nphotEC,nphot,npion,CD_proton,FD_proton,CD_neutron,FD_neutron,ncharged;  
  Int_t i_el,i_pr,i_n,i_n_cnd,i_phEC,i_phFT,i_ph;
  Int_t isortEC[40],isortIC[40];
  Int_t GoodECPhotIndex,GoodICPhotIndex;
  Int_t RunNumber,FileNumber;
  Int_t pol_sign,helicity;
  Int_t Ph_det;
  Float_t TarPol;
  Float_t Ebeam, Pmass, Nmass, Dmass,Elmass,LightSpeed;
  Float_t Q2, Xbj, W,t_Pr,t_Ph,t_N,El_P,Phi_Pr,Phi_Ph, Phi_N,Angl_X_g,Angl_hg_hp,Xbal,Ybal,Zbal,Ebal,Xbal_n,Ybal_n,Zbal_n,Ebal_n,El_Theta, El_Phi,El_vx,El_vy,El_vz,Pr_P,Pr_Theta,Pr_Phi,Pr_vz,N_P,N_Theta,N_Phi,Ph_EC_P[40],Ph_EC_Theta[40],Ph_EC_Phi[40],Ph_FT_P[40],Ph_FT_Theta[40],Ph_FT_Phi[40],Ph_P,Ph_Theta,Ph_Phi;
  Float_t mm2_epg,mm2_eng,mm2_ep,mm2_eg,mm2_eng_n,mm2_epg_p;
  Float_t theta_gamma_X,theta_gamma_e,theta_n_e;
  Float_t El_vx_MC,El_vy_MC,El_vz_MC;
  Float_t El_P_MC,El_Theta_MC,El_Phi_MC;
  Float_t Pr_P_MC,Pr_Theta_MC,Pr_Phi_MC;
  Float_t N_P_MC,N_Theta_MC,N_Phi_MC;
  Float_t Ph_P_MC,Ph_Theta_MC,Ph_Phi_MC; 
  Float_t Q2_MC, Xbj_MC, W_MC,t_Pr_MC,t_Ph_MC,t_N_MC,Phi_MC, theta_n_e_MC,theta_gamma_e_MC;
  Float_t tstart;
  Int_t n_spec,p_spec;
  Int_t veto;
  
  TLorentzVector ElectronBeam, Target_Vec,PTarget_Vec,NTarget_Vec;
  TLorentzVector El_Vec,Pr_Vec,N_Vec,Ph_IC_Vec[40],Ph_EC_Vec[40],Ph_Vec;
  TLorentzVector El_Vec_MC,Pr_Vec_MC,N_Vec_MC,Ph_Vec_MC;

   TVector3 pp;
   TVector3 PP[4];



   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> RUN_Config_run = {fReader, "RUN_Config_run"};
   TTreeReaderValue<Int_t> RUN_Config_event = {fReader, "RUN_Config_event"};
   TTreeReaderValue<Int_t> RUN_Config_unixtime = {fReader, "RUN_Config_unixtime"};
   TTreeReaderValue<Float_t> RUN_Config_trigger = {fReader, "RUN_Config_trigger"};
   TTreeReaderValue<Float_t> RUN_Config_timestamp = {fReader, "RUN_Config_timestamp"};
   TTreeReaderValue<Int_t> RUN_Config_type = {fReader, "RUN_Config_type"};
   TTreeReaderValue<Int_t> RUN_Config_mode = {fReader, "RUN_Config_mode"};
   TTreeReaderValue<Float_t> RUN_Config_torus = {fReader, "RUN_Config_torus"};
   TTreeReaderValue<Float_t> RUN_Config_solenoid = {fReader, "RUN_Config_solenoid"};
   TTreeReaderValue<Float_t> REC_Event_category = {fReader, "REC_Event_category"};
   TTreeReaderValue<Float_t> REC_Event_topology = {fReader, "REC_Event_topology"};
   TTreeReaderValue<Float_t> REC_Event_beamCharge = {fReader, "REC_Event_beamCharge"};
   TTreeReaderValue<Double_t> REC_Event_liveTime = {fReader, "REC_Event_liveTime"};
   TTreeReaderValue<Float_t> REC_Event_startTime = {fReader, "REC_Event_startTime"};
   TTreeReaderValue<Float_t> REC_Event_RFTime = {fReader, "REC_Event_RFTime"};
   TTreeReaderValue<Int_t> REC_Event_helicity = {fReader, "REC_Event_helicity"};
   TTreeReaderValue<Int_t> REC_Event_helicityRaw = {fReader, "REC_Event_helicityRaw"};
   TTreeReaderValue<Float_t> REC_Event_procTime = {fReader, "REC_Event_procTime"};
   TTreeReaderValue<Int_t> HEL_Online_helicity = {fReader, "HEL_Online_helicity"};
   TTreeReaderValue<Int_t> HEL_Online_helicityRaw = {fReader, "HEL_Online_helicityRaw"};
   TTreeReaderArray<int> REC_Particle_pid = {fReader, "REC_Particle_pid"};
   TTreeReaderArray<float> REC_Particle_px = {fReader, "REC_Particle_px"};
   TTreeReaderArray<float> REC_Particle_py = {fReader, "REC_Particle_py"};
   TTreeReaderArray<float> REC_Particle_pz = {fReader, "REC_Particle_pz"};
   TTreeReaderArray<float> REC_Particle_vx = {fReader, "REC_Particle_vx"};
   TTreeReaderArray<float> REC_Particle_vy = {fReader, "REC_Particle_vy"};
   TTreeReaderArray<float> REC_Particle_vz = {fReader, "REC_Particle_vz"};
   TTreeReaderArray<float> REC_Particle_vt = {fReader, "REC_Particle_vt"};
   TTreeReaderArray<int> REC_Particle_charge = {fReader, "REC_Particle_charge"};
   TTreeReaderArray<float> REC_Particle_beta = {fReader, "REC_Particle_beta"};
   TTreeReaderArray<float> REC_Particle_chi2pid = {fReader, "REC_Particle_chi2pid"};
   TTreeReaderArray<int> REC_Particle_status = {fReader, "REC_Particle_status"};
   TTreeReaderArray<float> __CND_Clusters_energy = {fReader, "__CND_Clusters_energy"};
   TTreeReaderArray<int> __CND_Clusters_id = {fReader, "__CND_Clusters_id"};
   TTreeReaderArray<int> __CND_Clusters_sector = {fReader, "__CND_Clusters_sector"};
   TTreeReaderArray<int> __CND_Clusters_layer = {fReader, "__CND_Clusters_layer"};
   TTreeReaderArray<int> __CND_Clusters_component = {fReader, "__CND_Clusters_component"};
   TTreeReaderArray<int> __CND_Clusters_nhits = {fReader, "__CND_Clusters_nhits"};
   TTreeReaderArray<int> __CND_Clusters_status = {fReader, "__CND_Clusters_status"};
   TTreeReaderArray<int> __CND_Clusters_veto = {fReader, "__CND_Clusters_veto"};
   TTreeReaderArray<int> __CND_Clusters_layermultip = {fReader, "__CND_Clusters_layermultip"};
   TTreeReaderArray<int> __CND_Clusters_layer1 = {fReader, "__CND_Clusters_layer1"};
   TTreeReaderArray<int> __CND_Clusters_layer2 = {fReader, "__CND_Clusters_layer2"};
   TTreeReaderArray<int> __CND_Clusters_layer3 = {fReader, "__CND_Clusters_layer3"};
   TTreeReaderArray<float> __CND_Clusters_x = {fReader, "__CND_Clusters_x"};
   TTreeReaderArray<float> __CND_Clusters_y = {fReader, "__CND_Clusters_y"};
   TTreeReaderArray<float> __CND_Clusters_z = {fReader, "__CND_Clusters_z"};
   TTreeReaderArray<float> __CND_Clusters_time = {fReader, "__CND_Clusters_time"};

   ana_clas12(TTree * /*tree*/ =0) { }
   virtual ~ana_clas12() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(ana_clas12,0);

};

#endif

#ifdef ana_clas12_cxx
void ana_clas12::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t ana_clas12::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef ana_clas12_cxx
