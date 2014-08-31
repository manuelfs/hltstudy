// smallntuple.C: Makes small ntuples from the HLT babies

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include "TCollection.h"
#include "TError.h" // Turns off "no dictionary for class" warnings

using namespace std;
using std::cout;
using std::endl;

vector<TString> dirlist(TString folder);
float cross_section(TString file);

void smallntuple(TString folder="/hadoop/cms/store/user/manuelf/HLT/"){
  gErrorIgnoreLevel=kError; // Turns off "no dictionary for class" warnings
  float onmet, onht;
  TChain chain("Events");
  TTree tree("tree", "tree");
  tree.Branch("onmet", &onmet);
  tree.Branch("onht", &onht);

  TString filename, rootfiles;
  vector<TString> dirs = dirlist(folder);
  //for(unsigned idir(0); idir < dirs.size(); idir++){
  for(unsigned idir(0); idir < 2; idir++){
    filename = "ntuples/"+dirs[idir]+".root";
    TFile rootfile(filename, "recreate");
    rootfiles = folder+"/"+dirs[idir]+"/*_80_*.root";
    chain.Add(rootfiles);
    chain.SetMakeClass(1);
    TBranch *b_ht = chain.GetBranch(chain.GetAlias("pf_ht"));
    TBranch *b_met = chain.GetBranch(chain.GetAlias("met_pt"));
    b_ht->SetAddress(&onht);
    b_met->SetAddress(&onmet);
    chain.SetMakeClass(0);
    long entries(chain.GetEntries());
    cout<<endl<<"Added "<<rootfiles<<" with "<<entries<<" entries"<<endl;
    //entries = 100;
    for(int entry(0); entry<entries; entry++){
      if(entry%10000==0) cout<<"Doing entry "<<entry<<" of "<<entries<<endl;
      b_ht->GetEntry(entry);b_met->GetEntry(entry);
      if(onht<=0) continue;
      tree.Fill();
    }
    tree.Write();
    rootfile.Close();
    cout<<"Written tree "<<filename<<endl;
    chain.Reset(); tree.Reset();
  }
}

// Returns list of directorites in folder
vector<TString> dirlist(TString folder){
  vector<TString> v_dirs;
  TSystemDirectory dir(folder, folder);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (file->IsDirectory() && !fname.Contains(".")) v_dirs.push_back(fname);
    }
  } // if(files)

  return v_dirs;
}

// Returns cross section of sample in pb
float cross_section(TString file){
  float xsec(0.);

  // From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVgluglu
  if(file.Contains("T1tttt") && file.Contains("825_"))   xsec = 1.2167;
  if(file.Contains("T1tttt") && file.Contains("1025_"))  xsec = 0.272778;
  if(file.Contains("T1tttt") && file.Contains("1150_"))  xsec = 0.117687;
  if(file.Contains("T1tttt") && file.Contains("1200_"))  xsec = 0.0856418;
  if(file.Contains("T1tttt") && file.Contains("1500_"))  xsec = 0.0141903;

  if(file.Contains("TT"))  xsec = 818.8;

  // From https://cms-pdmv.cern.ch/mcm
  if(file.Contains("WJetsToLNu_HT-100to200"))  xsec = 1817.0;
  if(file.Contains("WJetsToLNu_HT-200to400"))  xsec = 471.6;
  if(file.Contains("WJetsToLNu_HT-400to600"))  xsec = 55.61;
  if(file.Contains("WJetsToLNu_HT-600toInf"))  xsec = 18.81;

  if(file.Contains("QCD_Pt-5to10"))	   xsec = 80710000000;
  if(file.Contains("QCD_Pt-10to15"))	   xsec = 7528000000;
  if(file.Contains("QCD_Pt-15to30"))	   xsec = 2237000000;
  if(file.Contains("QCD_Pt-30to50"))	   xsec = 161500000;
  if(file.Contains("QCD_Pt-50to80"))	   xsec = 22110000;
  if(file.Contains("QCD_Pt-80to120"))	   xsec = 3000114;
  if(file.Contains("QCD_Pt-120to170"))   xsec = 493200;
  if(file.Contains("QCD_Pt-170to300"))   xsec = 120300;
  if(file.Contains("QCD_Pt-300to470"))   xsec = 7475;
  if(file.Contains("QCD_Pt-470to600"))   xsec = 587.1;
  if(file.Contains("QCD_Pt-600to800"))   xsec = 167;
  if(file.Contains("QCD_Pt-800to1000"))  xsec = 28.25;
  if(file.Contains("QCD_Pt-1000to1400")) xsec = 8.195;
  if(file.Contains("QCD_Pt-1400to1800")) xsec = 0.7346;
  if(file.Contains("QCD_Pt-1800to2400")) xsec = 0.102;
  if(file.Contains("QCD_Pt-2400to3200")) xsec = 0.00644;
  if(file.Contains("QCD_Pt-3200"))       xsec = 0.000163;

  return xsec;
}

