#include "TROOT.h"
#include "TRint.h"
#include "TTree.h"
#include <TH2.h>
#include "TSelector.h"
#include <TStyle.h>
#include <TString.h>
#include "TLorentzVector.h"
#include "ana_mc_clas12.h"
#include <fstream>    // for fstream
#include <iostream>   // for cout, endl
#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <stdio.h>
#include <string>

using namespace std;

int main(Int_t argc, Char_t *argv[])
{

ana_mc_clas12 *basic = new ana_mc_clas12();

TString OutputListName = "LAST_OUTPUT";
ofstream OutputList(OutputListName.Data());

TString FileNameIn = "runList";
Int_t Nlist;
Nlist =  atoi(argv[1]);

FileNameIn.Append(Form("%d",Nlist));
ifstream  ListFile(FileNameIn.Data());

Int_t Nruns;
Nruns =  atoi(argv[2]);

/*********************************************************/
/* d'abord on cherche h10 ou h25 dans le premier fichier */
/*********************************************************/

TString FileName;
ListFile>>FileName;
//FileName.Prepend("");
//FileName.Prepend("/work/clas/disk1/hyonsuk/dvcs/production/pass1/v1/root/ep/");
OutputList<<"Checking whether "<<FileName.Data()<<" is a good file."<<endl;
TFile* FirstFile = new TFile(FileName,"READ");
TList *FileContent;
FileContent = FirstFile->GetListOfKeys();
TChain *chain;
if(FileContent->Contains("clas12"))
   {
      cout << "making an hipo2root chain to analyse the " << Nruns << " runs"<< endl;
      chain = new TChain("clas12","Analysis chain");
      if(FirstFile){OutputList<<"Adding "<<FileName.Data()<<" to tree."<<endl;
      chain->Add(FileName.Data());}
   }
else 
   {
      cout << "Problem indetermining the the chain ID...fix problem " << endl;
      return 1;
}

/***************************************************************/
/* ensuite on ajoute tous les autres arbres a la chaine        */
/* en supposant qu'on ne travaille que avec des h?? identiques */
/* en particulier ici on commence a 1 et pas a 0 dans la liste */
/***************************************************************/

for(Int_t i=1;i<Nruns;i++)
  {
    TString FileName;
    ListFile>>FileName;

    /******************/
    /* chemin au CCPN */
    /******************/
//    FileName.Prepend("/sps/clas/fx/");
//    FileName.Prepend("rootFiles/");
    OutputList<<"Checking whether "<<FileName.Data()<<" is a good file."<<endl;
    TFile* file = new TFile(FileName,"READ");
    if(file){OutputList<<"Adding "<<FileName.Data()<<" to tree."<<endl;
      delete file;
      chain->Add(FileName.Data());}
  }



OutputList << "Number of events in chain : " << chain->GetEntries() << endl;

//chain->Process(basic,"",200000);
chain->Process(basic);

 return 0;
}
