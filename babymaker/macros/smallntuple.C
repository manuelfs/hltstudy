// smallntuple.C: Makes small ntuples from the HLT babies

#include "hlt_class.cc"
#include <utility>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>    // std::sort
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
#include "TSystem.h"

using namespace baby;
using namespace std;
using std::cout;
using std::endl;

vector<TString> dirlist(TString folder, TString inname="dir", TString tag="");
float cross_section(TString file);
void sortlists(int &nlist, vector<double> *pt, vector<double> *eta, vector<double> *phi,
	       const vector<float> treept, const vector<float> treeeta, const vector<float> treephi);
double deltaphi(double phi1, double phi2);
float dR(float eta1, float eta2, float phi1, float phi2);
void ngenleptons(TString filename, vector<int> &nori_genels, vector<int> &nori_genmus, vector<int> genlep_thresh);

void smallntuple(TString folder="/hadoop/cms/store/user/jaehyeok/HLT/", TString subfolder="test",
		 TString tagFolders="", int maxfiles=-1){
  gErrorIgnoreLevel=kError; // Turns off "no dictionary for class" warnings
  TString ntupledir = "ntuples/"+subfolder+"/";
  gSystem->mkdir(ntupledir);

  int totentries;
  vector<int> genlep_thresh, nori_genels, nori_genmus;
  float onmet, onmet_phi, onht, weight, wl1ht200, genht, genmet, mindr_mu, mindr2_mu, mindr_genmu;
  //vector<int_double> sorted; 
  int nels, ngenels, nmus, ngenmus, njets, ngenjets;
  vector<double> elspt, elseta, elsphi, genelspt, genelseta, genelsphi;
  vector<double> muspt, museta, musphi, genmuspt, genmuseta, genmusphi;
  vector<double> jetspt, jetseta, jetsphi, genjetspt, genjetseta, genjetsphi;
  const float luminosity = 19600;
  TChain chain("Events");
  TTree tree("tree", "tree");
  tree.Branch("mindr_mu", &mindr_mu);
  tree.Branch("mindr2_mu", &mindr2_mu);
  tree.Branch("mindr_genmu", &mindr_genmu);
  tree.Branch("onmet", &onmet);
  tree.Branch("onmet_phi", &onmet_phi);
  tree.Branch("onht", &onht);
  tree.Branch("weight", &weight);
  tree.Branch("wl1ht200", &wl1ht200);
  tree.Branch("genht", &genht);
  tree.Branch("genmet", &genmet);
  tree.Branch("els_pt", &elspt);
  tree.Branch("els_eta", &elseta);
  tree.Branch("els_phi", &elsphi);
  tree.Branch("genels_pt", &genelspt);
  tree.Branch("genels_eta", &genelseta);
  tree.Branch("genels_phi", &genelsphi);
  tree.Branch("mus_pt", &muspt);
  tree.Branch("mus_eta", &museta);
  tree.Branch("mus_phi", &musphi);
  tree.Branch("genmus_pt", &genmuspt);
  tree.Branch("genmus_eta", &genmuseta);
  tree.Branch("genmus_phi", &genmusphi);
  tree.Branch("jets_pt", &jetspt);
  tree.Branch("jets_eta", &jetseta);
  tree.Branch("jets_phi", &jetsphi);
  tree.Branch("genjets_pt", &genjetspt);
  tree.Branch("genjets_eta", &genjetseta);
  tree.Branch("genjets_phi", &genjetsphi);
  tree.Branch("nels", &nels);
  tree.Branch("ngenels", &ngenels);
  tree.Branch("nmus", &nmus);
  tree.Branch("ngenmus", &ngenmus);
  tree.Branch("njets", &njets);
  tree.Branch("ngenjets", &ngenjets);

  genlep_thresh.push_back(0); genlep_thresh.push_back(10); 
  genlep_thresh.push_back(15); genlep_thresh.push_back(17); 
  genlep_thresh.push_back(20); genlep_thresh.push_back(25); 
  for(unsigned ind(0); ind<genlep_thresh.size(); ind++){
    nori_genmus.push_back(0); nori_genels.push_back(0);
  }
  TTree treeglobal("treeglobal", "treeglobal");
  treeglobal.Branch("noriginal", &totentries);
  treeglobal.Branch("nori_genmus", &nori_genmus);
  treeglobal.Branch("nori_genels", &nori_genels);
  treeglobal.Branch("genlep_thresh", &genlep_thresh);

  TString filename, rootname, sampledir;
  vector<TString> dirs = dirlist(folder, "dir", tagFolders);
  for(unsigned idir(0); idir < dirs.size(); idir++){
    sampledir = folder+"/"+dirs[idir];
    chain.Add(sampledir+"/*.root");
    totentries = chain.GetEntries();
    int totentry(0);
    chain.Reset();
    float xsection = cross_section(dirs[idir]);

    vector<TString> rootfiles = dirlist(sampledir, ".root");
    cout<<endl<<"Adding the "<<rootfiles.size()<<" files in "<<sampledir
	<<" with "<<totentries<<" entries"<<endl;
    filename = ntupledir+dirs[idir]+".root";
    TFile rootfile(filename, "recreate");
    int nfiles = static_cast<int>(rootfiles.size());
    if(maxfiles>0 && maxfiles<nfiles) nfiles = maxfiles;
    for(int ifile(0); ifile < nfiles; ifile++){
      rootname = sampledir+"/"+rootfiles[ifile];
      chain.Add(rootname);
      hlt.Init(&chain);
      long entries(chain.GetEntries());
      for(int entry(0); entry<entries; entry++, totentry++){
	if((totentry+1)%500000==0) cout<<"Done "<<totentry+1<<" entries"<<endl;
	hlt.GetEntry(entry);
	// Saving only events with at least one lepton
	if(mus_pt().size()==0 && els_pt().size()==0) continue;
	onht = pf_ht();
	onmet = met_pt();
	onmet_phi = met_phi();
	genht = 0;
	for(unsigned ijet(0); ijet < genjets_pt().size(); ijet++)
	  if(genjets_pt().at(ijet)>40 && genjets_eta().at(ijet)<3) genht += genjets_pt().at(ijet);
	//genht = gen_ht();
	//genmet = gen_met();
	// Sort object lists in terms of pt, and save them
	sortlists(nmus, &muspt, &museta, &musphi, mus_pt(), mus_eta(), mus_phi());
	mindr_mu = 999;
	mindr2_mu = 999;
	if(museta.size() > 0){
	  float mueta(museta[0]), muphi(musphi[0]);
	  for(unsigned ijet(0); ijet < pfjets_pt().size(); ijet++){
	    if(pfjets_pt()[ijet] < 50) continue;
	    float dRjet(dR(mueta, pfjets_eta()[ijet], muphi, pfjets_phi()[ijet]));
	    if(dRjet < mindr2_mu) {
	      if(dRjet < mindr_mu){
		mindr2_mu = mindr_mu;
		mindr_mu = dRjet;
	      } else mindr2_mu = dRjet;
	    }
	  }
	}
	sortlists(ngenmus, &genmuspt, &genmuseta, &genmusphi, genmus_pt(), genmus_eta(), genmus_phi());
	mindr_genmu = 999;
	if(genmuseta.size() > 0){
	  float genmueta(genmuseta[0]), genmuphi(genmusphi[0]);
	  for(unsigned ijet(0); ijet < genjets_pt().size(); ijet++){
	    if(genjets_pt()[ijet] < 50) continue;
	    float dRjet(dR(genmueta, genjets_eta()[ijet], genmuphi, genjets_phi()[ijet]));
	    if(dRjet < mindr_genmu) mindr_genmu = dRjet;
	  }
	}
	sortlists(nels, &elspt, &elseta, &elsphi, els_pt(), els_eta(), els_phi());
	sortlists(ngenels, &genelspt, &genelseta, &genelsphi, genels_pt(), genels_eta(), genels_phi());
	sortlists(njets, &jetspt, &jetseta, &jetsphi, pfjets_pt(), pfjets_eta(), pfjets_phi());
	sortlists(ngenjets, &genjetspt, &genjetseta, &genjetsphi, genjets_pt(), genjets_eta(), genjets_phi());
	wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(genht-(3.02695e+02)))+0.5);
	weight = wl1ht200*xsection*luminosity / static_cast<double>(totentries);
	tree.Fill();
      }
      chain.Reset(); 
    } // Loop over files in the sample folder
    tree.Write();
    ngenleptons(filename, nori_genels, nori_genmus, genlep_thresh);
    treeglobal.Fill();
    treeglobal.Write();
    rootfile.Close();
    cout<<"Written tree "<<filename<<endl;
    tree.Reset(); treeglobal.Reset();
  } // Loop over sample folders
}

void ngenleptons(TString filename, vector<int> &nori_genels, vector<int> &nori_genmus, vector<int> genlep_thresh){
  if(filename.Contains("T1tttt")) {
    for(unsigned ind(0); ind<nori_genels.size(); ind++){
      nori_genels[ind] = 0; nori_genmus[ind] = 0;
    }
    TChain chain("Events");
    if(filename.Contains("1025_")) chain.Add("ntuples/babies/gen/ntuple_gen_T1tttt_1025.root");
    else chain.Add("ntuples/babies/gen/ntuple_gen_T1tttt_825.root");
    hlt_class hltgen;
    hltgen.Init(&chain);
    long entries(chain.GetEntries());
    for(int entry(0); entry<entries; entry++){
      hltgen.GetEntry(entry);
      for(unsigned ithresh(0); ithresh<genlep_thresh.size(); ithresh++){
	for(unsigned lep(0); lep<hltgen.genels_pt().size(); lep++)
	  if(hltgen.genels_pt()[lep] >= genlep_thresh[ithresh]) {
	    nori_genels[ithresh]++;
	    break;
	  }
	for(unsigned lep(0); lep<hltgen.genmus_pt().size(); lep++)
	  if(hltgen.genmus_pt()[lep] >= genlep_thresh[ithresh]) {
	    nori_genmus[ithresh]++;
	    break;
	  }
      } // Loop over lepton pT thresholds
    } // Loop over chain entries
  } else return;
}

// Definitions to sort vectors
typedef std::pair<int,double> int_double;
bool id_big2small(const int_double& left, const int_double& right){ 
  return left.second > right.second; 
}  

void sortlists(int &nlist, vector<double> *pt, vector<double> *eta, vector<double> *phi,
	       const vector<float> treept, const vector<float> treeeta, const vector<float> treephi){

  nlist = static_cast<int>(treept.size());
  vector<int_double> sorted; 
  pt->resize(0); eta->resize(0); phi->resize(0); sorted.resize(0);
  for(int ind(0); ind<nlist; ind++)
    sorted.push_back(make_pair(ind,treept[ind]));
  sort(sorted.begin(), sorted.end(), id_big2small);
  for(int ind(0); ind<nlist; ind++){
    pt->push_back(treept[sorted[ind].first]);
    eta->push_back(treeeta[sorted[ind].first]);
    phi->push_back(treephi[sorted[ind].first]);
  }
}

// Returns list of directorites or files in folder
vector<TString> dirlist(TString folder, TString inname, TString tag){
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
	if ((file->IsDirectory() && !fname.Contains(".") && fname.EndsWith(tag))) v_dirs.push_back(fname);
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

double deltaphi(double phi1, double phi2){
  double result = fabs(phi1-phi2);
  while (result>TMath::Pi()) result -= 2*TMath::Pi();
  while (result<=-TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return sqrt(pow(eta1-eta2, 2) + pow(deltaphi(phi1,phi2), 2)) ;
}

