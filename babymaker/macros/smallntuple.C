// smallntuple.C: Makes small ntuples from the HLT babies

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
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

vector<TString> dirlist(TString folder, TString inname="dir");
float cross_section(TString file);

void smallntuple(TString folder="/hadoop/cms/store/user/manuelf/HLT/"){
  gErrorIgnoreLevel=kError; // Turns off "no dictionary for class" warnings
  float onmet, onht, weight, wl1ht200, genht;
  vector<float> genjets_pt, genjets_eta;
  const float luminosity = 19600;
  TChain chain("Events");
  TTree tree("tree", "tree");
  tree.Branch("onmet", &onmet);
  tree.Branch("onht", &onht);
  tree.Branch("weight", &weight);
  tree.Branch("wl1ht200", &wl1ht200);
  tree.Branch("genht", &genht);

  TString filename, rootname, sampledir;
  vector<TString> dirs = dirlist(folder);
  for(unsigned idir(0); idir < dirs.size(); idir++){
    //for(unsigned idir(0); idir < 3; idir++){
    sampledir = folder+"/"+dirs[idir];
    chain.Add(sampledir+"/*.root");
    long totentries(chain.GetEntries()), totentry(0);
    chain.Reset();
    float xsection = cross_section(dirs[idir]);

    vector<TString> rootfiles = dirlist(sampledir, ".root");
    cout<<endl<<"Adding the "<<rootfiles.size()<<" files in "<<sampledir
	<<" with "<<totentries<<" entries"<<endl;
    filename = "ntuples/"+dirs[idir]+".root";
    TFile rootfile(filename, "recreate");
    // Using SetMakeClass == 1 can only be done in indivicual TTrees
    for(unsigned ifile(0); ifile < rootfiles.size(); ifile++){
      rootname = sampledir+"/"+rootfiles[ifile];
      chain.Add(rootname);
      chain.SetMakeClass(1);
      TBranch *b_ht = chain.GetBranch(chain.GetAlias("pf_ht"));
      TBranch *b_met = chain.GetBranch(chain.GetAlias("met_pt"));
      TBranch *b_genjets_pt = chain.GetBranch(chain.GetAlias("genjets_pt"));
      TBranch *b_genjets_eta = chain.GetBranch(chain.GetAlias("genjets_eta"));
      b_ht->SetAddress(&onht);
      b_met->SetAddress(&onmet);
      b_genjets_pt->SetAddress(&genjets_pt);
      b_genjets_eta->SetAddress(&genjets_eta);
      chain.SetMakeClass(0);
      long entries(chain.GetEntries());
      for(int entry(0); entry<entries; entry++, totentry++){
	if((totentry+1)%500000==0) cout<<"Done "<<totentry+1<<" entries"<<endl;
	b_ht->GetEntry(entry);b_met->GetEntry(entry);
	b_genjets_pt->GetEntry(entry);b_genjets_eta->GetEntry(entry);
	if(onht<=0) continue;
	genht = 0;
	for(unsigned ijet(0); ijet < genjets_pt.size(); ijet++) 
	  if(genjets_pt.at(ijet)>40 && genjets_eta.at(ijet)<3) genht += genjets_pt.at(ijet);
	wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(genht-(3.02695e+02)))+0.5);
	weight = wl1ht200*xsection*luminosity / static_cast<double>(totentries);
	tree.Fill();
      }
      chain.Reset(); 
    } // Loop over files in the sample folder
    tree.Write();
    rootfile.Close();
    cout<<"Written tree "<<filename<<endl;
    tree.Reset();
  } // Loop over sample folders
}

// Returns list of directorites or files in folder
vector<TString> dirlist(TString folder, TString inname){
  vector<TString> v_dirs;
  TSystemDirectory dir(folder, folder);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (inname=="dir") {
	if ((file->IsDirectory() && !fname.Contains("."))) v_dirs.push_back(fname);
      } else  if(fname.Contains(inname)) v_dirs.push_back(fname);
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

  if(file.Contains("QCD_Pt-5to10"))	 xsec = 80710000000;
  if(file.Contains("QCD_Pt-10to15"))	 xsec = 7528000000;
  if(file.Contains("QCD_Pt-15to30"))	 xsec = 2237000000;
  if(file.Contains("QCD_Pt-30to50"))	 xsec = 161500000;
  if(file.Contains("QCD_Pt-50to80"))	 xsec = 22110000;
  if(file.Contains("QCD_Pt-80to120"))	 xsec = 3000114;
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

