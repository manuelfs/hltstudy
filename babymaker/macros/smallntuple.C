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

typedef std::pair<int,double> int_double;
vector<TString> dirlist(TString folder, TString inname="dir", TString tag="");
float cross_section(TString file);
vector<int_double> sortlists(int &nlist, vector<double> *pt, vector<double> *eta, vector<double> *phi,
			     const vector<float> treept, const vector<float> treeeta, const vector<float> treephi);
double deltaphi(double phi1, double phi2);
float dR(float eta1, float eta2, float phi1, float phi2);
void ngenleptons(TString filename, vector<int> &nori_genels, vector<int> &nori_genmus, vector<int> genlep_thresh);

void smallntuple(TString folder="/hadoop/cms/store/user/manuelf/HLT/", TString subfolder="test",
		 TString tagFolders="", int maxfiles=-1){
  gErrorIgnoreLevel=kError; // Turns off "no dictionary for class" warnings
  TString ntupledir = "ntuples/"+subfolder+"/";
  gSystem->mkdir(ntupledir);

  int totentries, nori_genmu0, nori_genel0;
  vector<int> genlep_thresh, v_nori_genels, v_nori_genmus;
  vector<int> *nori_genels = &v_nori_genels, *nori_genmus = &v_nori_genmus;
  float onmet, onmet_phi, onht, caloht, weight, wl1ht200, genht, genmet, wlumi;
  //vector<int_double> sorted; 
  int nels, ngenels, nmus, ngenmus, njets, nbjets, ngenjets;
  vector<double> elspt, elseta, elsphi, genelspt, genelseta, genelsphi;
  vector<double> elsclustershape, elshe, elseminusp, elsdeta, elsdphi, elsmindr;
  vector<double> muspt, museta, musphi, genmuspt, genmuseta, genmusphi, musmindr;
  vector<double> musreliso, musgenpt, elsreliso, elstrackiso, elsecaliso, elshcaliso, elsgenpt;
  vector<double> jetspt, jetseta, jetsphi, genjetspt, genjetseta, genjetsphi;
  vector<double> bjetspt, bjetseta, bjetsphi, bjetscsv;
  bool not_pu;
  const float luminosity = 19600, mindrcut = 0.4;
  TChain chain("Events");
  TTree tree("tree", "tree");
  tree.Branch("onmet", &onmet);
  tree.Branch("onmet_phi", &onmet_phi);
  tree.Branch("onht", &onht);
  tree.Branch("caloht", &caloht);
  tree.Branch("weight", &weight);
  tree.Branch("wl1ht200", &wl1ht200);
  tree.Branch("wlumi", &wlumi);
  tree.Branch("genht", &genht);
  tree.Branch("genmet", &genmet);
  tree.Branch("els_pt", &elspt);
  tree.Branch("els_eta", &elseta);
  tree.Branch("els_phi", &elsphi);
  tree.Branch("els_reliso", &elsreliso);
  tree.Branch("els_trackiso", &elstrackiso);
  tree.Branch("els_ecaliso", &elsecaliso);
  tree.Branch("els_hcaliso", &elshcaliso);
  tree.Branch("els_genpt", &elsgenpt);
  tree.Branch("els_clustershape", &elsclustershape);
  tree.Branch("els_he", &elshe);
  tree.Branch("els_eminusp", &elseminusp);
  tree.Branch("els_deta", &elsdeta);
  tree.Branch("els_dphi", &elsdphi);
  tree.Branch("els_mindr", &elsmindr);
  tree.Branch("genels_pt", &genelspt);
  tree.Branch("genels_eta", &genelseta);
  tree.Branch("genels_phi", &genelsphi);
  tree.Branch("mus_pt", &muspt);
  tree.Branch("mus_eta", &museta);
  tree.Branch("mus_phi", &musphi);
  tree.Branch("mus_reliso", &musreliso);
  tree.Branch("mus_genpt", &musgenpt);
  tree.Branch("mus_mindr", &musmindr);
  tree.Branch("genmus_pt", &genmuspt);
  tree.Branch("genmus_eta", &genmuseta);
  tree.Branch("genmus_phi", &genmusphi);
  tree.Branch("jets_pt", &jetspt);
  tree.Branch("jets_eta", &jetseta);
  tree.Branch("jets_phi", &jetsphi);
  tree.Branch("bjets_pt", &bjetspt);
  tree.Branch("bjets_eta", &bjetseta);
  tree.Branch("bjets_phi", &bjetsphi);
  tree.Branch("bjets_csv", &bjetscsv);
  tree.Branch("genjets_pt", &genjetspt);
  tree.Branch("genjets_eta", &genjetseta);
  tree.Branch("genjets_phi", &genjetsphi);
  tree.Branch("nels", &nels);
  tree.Branch("ngenels", &ngenels);
  tree.Branch("nmus", &nmus);
  tree.Branch("ngenmus", &ngenmus);
  tree.Branch("njets", &njets);
  tree.Branch("nbjets", &nbjets);
  tree.Branch("ngenjets", &ngenjets);
  tree.Branch("not_pu", &not_pu);

  genlep_thresh.push_back(0); genlep_thresh.push_back(10); 
  genlep_thresh.push_back(15); genlep_thresh.push_back(17); 
  genlep_thresh.push_back(20); genlep_thresh.push_back(25); 
  for(unsigned ind(0); ind<genlep_thresh.size(); ind++){
    v_nori_genmus.push_back(0); v_nori_genels.push_back(0);
  }
  TTree treeglobal("treeglobal", "treeglobal");
  treeglobal.Branch("noriginal", &totentries);
  treeglobal.Branch("nori_genmus", &nori_genmus);
  treeglobal.Branch("nori_genels", &nori_genels);
  treeglobal.Branch("nori_genmu0", &nori_genmu0);
  treeglobal.Branch("nori_genel0", &nori_genel0);
  treeglobal.Branch("genlep_thresh", &genlep_thresh);

  vector<int_double> mus_sorted, els_sorted, bjets_sorted;
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

	caloht = calo_ht();
	onht = pf_ht();
	onmet = met_pt();
	onmet_phi = met_phi();
	genht = 0;
	// If the leading jet is not matched to the leading genjet the event is likely to be PU
	// https://indico.cern.ch/event/351559/contribution/4/material/slides/0.pdf
	if(genjets_eta().size()>0 && pfjets_eta().size()>0)
	  not_pu = (sqrt(pow(genjets_eta().at(0)-pfjets_eta().at(0),2)+
			 pow(genjets_phi().at(0)-pfjets_phi().at(0),2)) < 0.1);
	else not_pu = false;
	for(unsigned ijet(0); ijet < genjets_pt().size(); ijet++)
	  if(genjets_pt().at(ijet)>40 && genjets_eta().at(ijet)<3) genht += genjets_pt().at(ijet);
	genmet = gen_met();

	// Sort object lists in terms of pt, and save them
	mus_sorted = sortlists(nmus, &muspt, &museta, &musphi, mus_pt(), mus_eta(), mus_phi());
	musreliso.resize(0); musgenpt.resize(0); musmindr.resize(0);
	for(int ilep(0); ilep < nmus; ilep++){
	  musreliso.push_back(mus_iso()[mus_sorted[ilep].first]);
	  float mindr(999.);
	  int imin(-1);
	  for(unsigned igen(0); igen < genmus_pt().size(); igen++){
	    if(abs(genmus_mom_id()[igen]) == 24 || abs(genmus_gmom_id()[igen]) == 24 || 
	       abs(genmus_ggmom_id()[igen]) == 24){
	      float dr = abs(dR(museta[ilep], genmus_eta()[igen], musphi[ilep], genmus_phi()[igen]));
	      float mingenpt = genmus_pt()[igen];
	      if(dr < mindr && dr < mindrcut && abs((mingenpt-muspt[ilep])/muspt[ilep])<0.3) 
		{mindr = dr; imin = igen;}
	    }
	  } // Loop over genmus
	  if(imin>=0) musgenpt.push_back(genmus_pt()[imin]);
	  else musgenpt.push_back(-99);
	  mindr = 999.;
	  for(unsigned igen(0); igen < gentop_pt().size(); igen++){
	    float dr = dR(museta[ilep], gentop_eta()[igen], musphi[ilep], gentop_phi()[igen]);
	    if(dr < mindr) mindr = dr;
	  }
	  musmindr.push_back(mindr);
	}
	sortlists(ngenmus, &genmuspt, &genmuseta, &genmusphi, genmus_pt(), genmus_eta(), genmus_phi());
	els_sorted = sortlists(nels, &elspt, &elseta, &elsphi, els_pt(), els_eta(), els_phi());
	elsreliso.resize(0); elstrackiso.resize(0); elsecaliso.resize(0); elshcaliso.resize(0);
	elsclustershape.resize(0); elshe.resize(0); elseminusp.resize(0); elsdeta.resize(0); elsdphi.resize(0); 
	elsgenpt.resize(0); elsmindr.resize(0);
	for(int ilep(0); ilep < nels; ilep++){
	  elstrackiso.push_back(els_track_iso()[els_sorted[ilep].first]);
	  elsecaliso.push_back(els_ecal_iso()[els_sorted[ilep].first]);
	  elshcaliso.push_back(els_hcal_iso()[els_sorted[ilep].first]);
	  elsreliso.push_back(elstrackiso[ilep]+elsecaliso[ilep]+elshcaliso[ilep]);
	  elsclustershape.push_back(els_clustershape()[els_sorted[ilep].first]);
	  elshe.push_back(els_he()[els_sorted[ilep].first]);
	  elseminusp.push_back(els_eminusp()[els_sorted[ilep].first]);
	  elsdeta.push_back(els_deta()[els_sorted[ilep].first]);
	  elsdphi.push_back(els_dphi()[els_sorted[ilep].first]);
	  float mindr(999.);
	  int imin(-1);
	  for(unsigned igen(0); igen < genels_pt().size(); igen++){
	    if(abs(genels_mom_id()[igen]) == 24 || abs(genels_gmom_id()[igen]) == 24 || 
	       abs(genels_ggmom_id()[igen]) == 24){
	      float dr = abs(dR(elseta[ilep], genels_eta()[igen], elsphi[ilep], genels_phi()[igen]));
	      float mingenpt = genels_pt()[igen];
	      if(dr < mindr && dr < mindrcut && abs((mingenpt-elspt[ilep])/elspt[ilep])<0.3) 
		{mindr = dr; imin = igen;}
	    }
	  } // Loop over genels
	  if(imin>=0) elsgenpt.push_back(genels_pt()[imin]);
	  else elsgenpt.push_back(-99);
	  mindr = 999.;
	  for(unsigned igen(0); igen < gentop_pt().size(); igen++){
	    float dr = dR(elseta[ilep], gentop_eta()[igen], elsphi[ilep], gentop_phi()[igen]);
	    if(dr < mindr) mindr = dr;
	  }
	  elsmindr.push_back(mindr);
	}

	sortlists(ngenels, &genelspt, &genelseta, &genelsphi, genels_pt(), genels_eta(), genels_phi());
	sortlists(njets, &jetspt, &jetseta, &jetsphi, pfjets_pt(), pfjets_eta(), pfjets_phi());
	bjets_sorted = sortlists(nbjets, &bjetspt, &bjetseta, &bjetsphi, bjets_pt(), bjets_eta(), bjets_phi());
	bjetscsv.resize(0);
	for(int ilep(0); ilep < nbjets; ilep++)
	  bjetscsv.push_back(bjets_csv()[bjets_sorted[ilep].first]);
	sortlists(ngenjets, &genjetspt, &genjetseta, &genjetsphi, genjets_pt(), genjets_eta(), genjets_phi());
	wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(genht-(3.02695e+02)))+0.5);
	wlumi = xsection*luminosity / static_cast<double>(totentries);
	weight = wl1ht200*wlumi;
	tree.Fill();
      }
      chain.Reset(); 
    } // Loop over files in the sample folder
    tree.Write();
    ngenleptons(filename, v_nori_genels, v_nori_genmus, genlep_thresh);
    nori_genel0 = v_nori_genels[0]; nori_genmu0 = v_nori_genmus[0]; 
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
    TString gluinoMass(""), folder("ntuples/babies/gen/");
    if(filename.Contains("1025_")) gluinoMass = "1025";
    else if(filename.Contains("825_")) gluinoMass = "825";
    else if(filename.Contains("1500_")) gluinoMass = "1500_";
    else if(filename.Contains("1200_")) gluinoMass = "1200_";
    else {cout<<"No gen file for "<<filename<<". Exiting"<<endl; return;}

    vector<TString> rootfiles = dirlist(folder, gluinoMass);
    int nfiles = static_cast<int>(rootfiles.size());
    for(int ifile(0); ifile < nfiles; ifile++){
      TChain chain("Events");
      TString rootname = folder+"/"+rootfiles[ifile];
      chain.Add(rootname);
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
    } // Loop over signal files
  } else return;
}

// Definitions to sort vectors
bool id_big2small(const int_double& left, const int_double& right){ 
  return left.second > right.second; 
}  

vector<int_double> sortlists(int &nlist, vector<double> *pt, vector<double> *eta, vector<double> *phi,
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
  return sorted;
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

  if(file.Contains("WToENu"))   xsec = 16000.0;
  if(file.Contains("WToMuNu"))  xsec = 16100.0;

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

