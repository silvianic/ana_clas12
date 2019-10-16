//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 26 11:09:57 2018 by ROOT version 6.02/02
// from TTree hipo2root/CLAS12 banks in ROOT
// found on file: out__171218_095035_events_deut_run31127632.dat.root
//////////////////////////////////////////////////////////

#ifndef ana_mc_clas12_h
#define ana_mc_clas12_h

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

class ana_mc_clas12 : public TSelector {
public :
  //   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree *el_tree,*pDVCS_tree,*nDVCS_tree; 
   TTree *MC_tree; 
   TFile *outfile; 
   Bool_t   firstCall; 
   void    AddBranches();
   void    AddMCBranches();
   void    Calc_kine_P();
   void    Calc_kine_N();   
   void    Calc_kine_MC();


// Fixed size dimensions of array or collections stored in the TTree if any.
   Int_t nelec,nprot,nneut,nneut_matched,nphotFT,nphotEC,nphot,npion,CD_proton,FD_proton,CD_neutron,FD_neutron,ncharged;  
  Int_t i_el,i_pr,i_n,i_n_cnd,i_phEC,i_phFT,i_ph;
  Int_t isortEC[40],isortIC[40];
  Int_t GoodECPhotIndex,GoodICPhotIndex;
  Int_t RunNumber,FileNumber;
  Int_t pol_sign,helicity;
  Int_t Ph_det;
  Float_t TarPol;
  Float_t Ebeam, Pmass, Nmass, Dmass, Elmass, LightSpeed;
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
   /* TTreeReaderArray<int> MC_Header_run = {fReader, "MC_Header_run"}; */
   /* TTreeReaderArray<int> MC_Header_event = {fReader, "MC_Header_event"}; */
   /* TTreeReaderArray<int> MC_Header_type = {fReader, "MC_Header_type"}; */
   /* TTreeReaderArray<float> MC_Header_helicity = {fReader, "MC_Header_helicity"}; */
   /* TTreeReaderArray<int> MC_Event_npart = {fReader, "MC_Event_npart"}; */
   /* TTreeReaderArray<float> MC_Event_ebeam = {fReader, "MC_Event_ebeam"}; */
   /* TTreeReaderArray<float> MC_Event_weight = {fReader, "MC_Event_weight"}; */
   TTreeReaderArray<int> MC_Particle_pid = {fReader, "MC_Particle_pid"};
   TTreeReaderArray<float> MC_Particle_px = {fReader, "MC_Particle_px"};
   TTreeReaderArray<float> MC_Particle_py = {fReader, "MC_Particle_py"};
   TTreeReaderArray<float> MC_Particle_pz = {fReader, "MC_Particle_pz"};
   TTreeReaderArray<float> MC_Particle_vx = {fReader, "MC_Particle_vx"};
   TTreeReaderArray<float> MC_Particle_vy = {fReader, "MC_Particle_vy"};
   TTreeReaderArray<float> MC_Particle_vz = {fReader, "MC_Particle_vz"};
   TTreeReaderArray<float> MC_Particle_vt = {fReader, "MC_Particle_vt"};
   /* TTreeReaderArray<int> MC_Lund_pid = {fReader, "MC_Lund_pid"}; */
   /* TTreeReaderArray<float> MC_Lund_mass = {fReader, "MC_Lund_mass"}; */
   /* TTreeReaderArray<float> MC_Lund_E = {fReader, "MC_Lund_E"}; */
   /* TTreeReaderArray<float> MC_Lund_px = {fReader, "MC_Lund_px"}; */
   /* TTreeReaderArray<float> MC_Lund_py = {fReader, "MC_Lund_py"}; */
   /* TTreeReaderArray<float> MC_Lund_pz = {fReader, "MC_Lund_pz"}; */
   /* TTreeReaderArray<float> MC_Lund_vx = {fReader, "MC_Lund_vx"}; */
   /* TTreeReaderArray<float> MC_Lund_vy = {fReader, "MC_Lund_vy"}; */
   /* TTreeReaderArray<float> MC_Lund_vz = {fReader, "MC_Lund_vz"}; */
   TTreeReaderArray<int> REC_Calorimeter_index = {fReader, "REC_Calorimeter_index"};
   TTreeReaderArray<int> REC_Calorimeter_pindex = {fReader, "REC_Calorimeter_pindex"};
   TTreeReaderArray<int> REC_Calorimeter_detector = {fReader, "REC_Calorimeter_detector"};
   TTreeReaderArray<int> REC_Calorimeter_sector = {fReader, "REC_Calorimeter_sector"};
   TTreeReaderArray<int> REC_Calorimeter_layer = {fReader, "REC_Calorimeter_layer"};
   TTreeReaderArray<float> REC_Calorimeter_energy = {fReader, "REC_Calorimeter_energy"};
   TTreeReaderArray<float> REC_Calorimeter_time = {fReader, "REC_Calorimeter_time"};
   TTreeReaderArray<float> REC_Calorimeter_path = {fReader, "REC_Calorimeter_path"};
   TTreeReaderArray<float> REC_Calorimeter_chi2 = {fReader, "REC_Calorimeter_chi2"};
   TTreeReaderArray<float> REC_Calorimeter_x = {fReader, "REC_Calorimeter_x"};
   TTreeReaderArray<float> REC_Calorimeter_y = {fReader, "REC_Calorimeter_y"};
   TTreeReaderArray<float> REC_Calorimeter_z = {fReader, "REC_Calorimeter_z"};
   TTreeReaderArray<float> REC_Calorimeter_hx = {fReader, "REC_Calorimeter_hx"};
   TTreeReaderArray<float> REC_Calorimeter_hy = {fReader, "REC_Calorimeter_hy"};
   TTreeReaderArray<float> REC_Calorimeter_hz = {fReader, "REC_Calorimeter_hz"};
   TTreeReaderArray<float> REC_Calorimeter_lu = {fReader, "REC_Calorimeter_lu"};
   TTreeReaderArray<float> REC_Calorimeter_lv = {fReader, "REC_Calorimeter_lv"};
   TTreeReaderArray<float> REC_Calorimeter_lw = {fReader, "REC_Calorimeter_lw"};
   TTreeReaderArray<float> REC_Calorimeter_du = {fReader, "REC_Calorimeter_du"};
   TTreeReaderArray<float> REC_Calorimeter_dv = {fReader, "REC_Calorimeter_dv"};
   TTreeReaderArray<float> REC_Calorimeter_dw = {fReader, "REC_Calorimeter_dw"};
   TTreeReaderArray<float> REC_Calorimeter_m2u = {fReader, "REC_Calorimeter_m2u"};
   TTreeReaderArray<float> REC_Calorimeter_m2v = {fReader, "REC_Calorimeter_m2v"};
   TTreeReaderArray<float> REC_Calorimeter_m2w = {fReader, "REC_Calorimeter_m2w"};
   TTreeReaderArray<float> REC_Calorimeter_m3u = {fReader, "REC_Calorimeter_m3u"};
   TTreeReaderArray<float> REC_Calorimeter_m3v = {fReader, "REC_Calorimeter_m3v"};
   TTreeReaderArray<float> REC_Calorimeter_m3w = {fReader, "REC_Calorimeter_m3w"};
   TTreeReaderArray<int> REC_Calorimeter_status = {fReader, "REC_Calorimeter_status"};
   TTreeReaderArray<int> REC_Cherenkov_index = {fReader, "REC_Cherenkov_index"};
   TTreeReaderArray<int> REC_Cherenkov_pindex = {fReader, "REC_Cherenkov_pindex"};
   TTreeReaderArray<int> REC_Cherenkov_detector = {fReader, "REC_Cherenkov_detector"};
   TTreeReaderArray<int> REC_Cherenkov_sector = {fReader, "REC_Cherenkov_sector"};
   TTreeReaderArray<float> REC_Cherenkov_nphe = {fReader, "REC_Cherenkov_nphe"};
   TTreeReaderArray<float> REC_Cherenkov_time = {fReader, "REC_Cherenkov_time"};
   TTreeReaderArray<float> REC_Cherenkov_path = {fReader, "REC_Cherenkov_path"};
   TTreeReaderArray<float> REC_Cherenkov_chi2 = {fReader, "REC_Cherenkov_chi2"};
   TTreeReaderArray<float> REC_Cherenkov_x = {fReader, "REC_Cherenkov_x"};
   TTreeReaderArray<float> REC_Cherenkov_y = {fReader, "REC_Cherenkov_y"};
   TTreeReaderArray<float> REC_Cherenkov_z = {fReader, "REC_Cherenkov_z"};
   TTreeReaderArray<float> REC_Cherenkov_dtheta = {fReader, "REC_Cherenkov_dtheta"};
   TTreeReaderArray<float> REC_Cherenkov_dphi = {fReader, "REC_Cherenkov_dphi"};
   TTreeReaderArray<int> REC_Cherenkov_status = {fReader, "REC_Cherenkov_status"};
   TTreeReaderArray<int> REC_ForwardTagger_index = {fReader, "REC_ForwardTagger_index"};
   TTreeReaderArray<int> REC_ForwardTagger_pindex = {fReader, "REC_ForwardTagger_pindex"};
   TTreeReaderArray<int> REC_ForwardTagger_detector = {fReader, "REC_ForwardTagger_detector"};
   TTreeReaderArray<int> REC_ForwardTagger_layer = {fReader, "REC_ForwardTagger_layer"};
   TTreeReaderArray<float> REC_ForwardTagger_energy = {fReader, "REC_ForwardTagger_energy"};
   TTreeReaderArray<float> REC_ForwardTagger_time = {fReader, "REC_ForwardTagger_time"};
   TTreeReaderArray<float> REC_ForwardTagger_path = {fReader, "REC_ForwardTagger_path"};
   TTreeReaderArray<float> REC_ForwardTagger_chi2 = {fReader, "REC_ForwardTagger_chi2"};
   TTreeReaderArray<float> REC_ForwardTagger_x = {fReader, "REC_ForwardTagger_x"};
   TTreeReaderArray<float> REC_ForwardTagger_y = {fReader, "REC_ForwardTagger_y"};
   TTreeReaderArray<float> REC_ForwardTagger_z = {fReader, "REC_ForwardTagger_z"};
   TTreeReaderArray<float> REC_ForwardTagger_dx = {fReader, "REC_ForwardTagger_dx"};
   TTreeReaderArray<float> REC_ForwardTagger_dy = {fReader, "REC_ForwardTagger_dy"};
   TTreeReaderArray<float> REC_ForwardTagger_radius = {fReader, "REC_ForwardTagger_radius"};
   TTreeReaderArray<int> REC_ForwardTagger_size = {fReader, "REC_ForwardTagger_size"};
   TTreeReaderArray<int> REC_ForwardTagger_status = {fReader, "REC_ForwardTagger_status"};
   TTreeReaderArray<int> REC_Scintillator_index = {fReader, "REC_Scintillator_index"};
   TTreeReaderArray<int> REC_Scintillator_pindex = {fReader, "REC_Scintillator_pindex"};
   TTreeReaderArray<int> REC_Scintillator_detector = {fReader, "REC_Scintillator_detector"};
   TTreeReaderArray<int> REC_Scintillator_sector = {fReader, "REC_Scintillator_sector"};
   TTreeReaderArray<int> REC_Scintillator_layer = {fReader, "REC_Scintillator_layer"};
   TTreeReaderArray<int> REC_Scintillator_component = {fReader, "REC_Scintillator_component"};
   TTreeReaderArray<float> REC_Scintillator_energy = {fReader, "REC_Scintillator_energy"};
   TTreeReaderArray<float> REC_Scintillator_time = {fReader, "REC_Scintillator_time"};
   TTreeReaderArray<float> REC_Scintillator_path = {fReader, "REC_Scintillator_path"};
   TTreeReaderArray<float> REC_Scintillator_chi2 = {fReader, "REC_Scintillator_chi2"};
   TTreeReaderArray<float> REC_Scintillator_x = {fReader, "REC_Scintillator_x"};
   TTreeReaderArray<float> REC_Scintillator_y = {fReader, "REC_Scintillator_y"};
   TTreeReaderArray<float> REC_Scintillator_z = {fReader, "REC_Scintillator_z"};
   TTreeReaderArray<float> REC_Scintillator_hx = {fReader, "REC_Scintillator_hx"};
   TTreeReaderArray<float> REC_Scintillator_hy = {fReader, "REC_Scintillator_hy"};
   TTreeReaderArray<float> REC_Scintillator_hz = {fReader, "REC_Scintillator_hz"};
   TTreeReaderArray<int> REC_Scintillator_status = {fReader, "REC_Scintillator_status"};
   TTreeReaderArray<int> REC_Track_index = {fReader, "REC_Track_index"};
   TTreeReaderArray<int> REC_Track_pindex = {fReader, "REC_Track_pindex"};
   TTreeReaderArray<int> REC_Track_detector = {fReader, "REC_Track_detector"};
   TTreeReaderArray<int> REC_Track_sector = {fReader, "REC_Track_sector"};
   TTreeReaderArray<int> REC_Track_status = {fReader, "REC_Track_status"};
   TTreeReaderArray<int> REC_Track_q = {fReader, "REC_Track_q"};
   TTreeReaderArray<float> REC_Track_chi2 = {fReader, "REC_Track_chi2"};
   TTreeReaderArray<int> REC_Track_NDF = {fReader, "REC_Track_NDF"};
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

   ana_mc_clas12(TTree * /*tree*/ =0) { }
   virtual ~ana_mc_clas12() { }
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

   ClassDef(ana_mc_clas12,0);
};

#endif

#ifdef ana_mc_clas12_cxx
void ana_mc_clas12::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

 
   fReader.SetTree(tree);

}


Bool_t ana_mc_clas12::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ana_clas12_cxx
