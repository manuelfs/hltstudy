// rates_hlt.C: Compares HLT rates of lepton pt, HT, and MET cuts

#define INT_ROOT
#include "styles.hpp"
#include "styles.cpp"
#include "ucsb_utils.hpp"
#include "ucsb_utils.cpp"


#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include "TChain.h"
#include "TFile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TBox.h"

#define NCuts 1

using namespace std;
using std::cout;
using std::endl;

void ReadChains(vector<TChain*> &chain, vector<int> &entries, TString folder, vector<TString> filenames);

void rates_hlt(TString isocut = "0.8", TString folder="root/15-03-14/el15/", 
	       bool applyCSV=false, bool onlyRates=false){
  styles style("2Dtitle"); style.setDefaultStyle(); gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetHatchesLineWidth(2);
  TCanvas can;
  bool isEl = folder.Contains("el15")?true:false;

  //Files
  vector<TString> filenames, Htag;
  filenames.push_back("/*QCD*");	        Htag.push_back("QCD");
  filenames.push_back("/*TT*");		        Htag.push_back("tt");
  filenames.push_back("/W*");                   Htag.push_back("Wjets");
  filenames.push_back("/*T1tttt*1500_*PU20*");  Htag.push_back("sig1500");
  filenames.push_back("/*T1tttt*1025_*");       Htag.push_back("sig1025");
  //filenames.push_back("/*T1tttt*1200_*PU20*");  Htag.push_back("sig1200");
  Htag.push_back("Total");
  vector<TChain*> chain;
  vector<int> noriginal;
  ReadChains(chain, noriginal, folder, filenames);
  if(chain.size()==0) return;
  vector<int> indchain;
  indchain.push_back(0);
  indchain.push_back(1);
  indchain.push_back(2);
  indchain.push_back(3);
  indchain.push_back(4);
  indchain.push_back(Htag.size()-1);
  const int NSam = indchain.size();

  // Histograms
  TString Hname, totCut, Pnamebase, Pname, Title;
  TString xTitle = "Minimum HLT PF H_{T} (GeV)";
  TString yTitle = "Minimum HLT PF E_{T,miss} (GeV)";
  TString zTitle = "HLT rate (Hz)";
  TH2D *hRate[2][NCuts][20];
  TBox box; box.SetLineColor(9); box.SetLineWidth(4); box.SetFillStyle(0); 

  // Cuts
  TString Cuts[] = {"Max$(mus_pt*((mus_reliso)<isocut))>15&&not_pu"};
  //TString Cuts[] = {"Max$(mus_pt*((mus_reliso)<isocut))>15"};
  TString TitleCuts[] = {", RelIso_{R=0.2} < isocut"};
  if(isEl)
    //Cuts[0] = "Max$(els_pt*((els_trackiso+els_hcaliso+els_ecaliso)<isocut))>15&&not_pu";
    Cuts[0] = "Max$(els_pt*(els_trackiso<0.4*isocut/0.8&&els_hcaliso<0.6*isocut/0.8&&els_ecaliso<0.5*isocut/0.8))>15&&not_pu";
  //Cuts[0] = "Max$(els_pt*((els_trackiso+els_hcaliso+els_ecaliso)<isocut))>15";
  Cuts[0].ReplaceAll("isocut",isocut);
  // Cuts for precise rates
  vector<TString> cutTrigger, nameTrigger;
  cutTrigger.push_back("onht>=600"); nameTrigger.push_back("HT600");
  cutTrigger.push_back("onht>=400&&abs(onmet)>=70"); nameTrigger.push_back("HT400_MET70");
  cutTrigger.push_back("onht>=400&&Max$(bjets_csv*(bjets_pt>35))>0.7"); nameTrigger.push_back("HT400_CSV");
  cutTrigger.push_back("(0"); nameTrigger.push_back("\tAll");
  for(unsigned icut(0); icut<cutTrigger.size()-1; icut++)
    cutTrigger[cutTrigger.size()-1] += "||"+cutTrigger[icut];
  cutTrigger[cutTrigger.size()-1] += ")";
  for(unsigned icut(0); icut<cutTrigger.size(); icut++)
    cutTrigger[icut] = "1.4e-2/19600*weight*(" + Cuts[0] + "&&" + cutTrigger[icut] + ")";

  if(applyCSV){
    Cuts[0] += "&&Max$(bjets_csv*(bjets_pt>35))>0.7";
    TitleCuts[0] += ", 1 b-tag";
  }
  TitleCuts[0].ReplaceAll("isocut",isocut);
  isocut.ReplaceAll(".","p");


  TString Tags[] = {(isEl?"el_":"mu_")+isocut+(applyCSV?"_btag":"")};

  TString tagfolder = folder;
  if(tagfolder[tagfolder.Sizeof()-2] == '/') tagfolder.Remove(tagfolder.Sizeof()-2);
  tagfolder.Remove(0,tagfolder.Last('/')+1);

  int nBinsHt = 14, nBinsMet = 14;
  float minHt=200, maxHt=900, minMet=0, maxMet=140;

  /////////////////// Calculating specific trigger rates  ////////////////
  TH1D histo1d("histo1d","",1,0,1);
  vector<double> trigger_rate(cutTrigger.size(),0.);
  cout<<endl;
  for(unsigned icut(0); icut<cutTrigger.size(); icut++){
    for(int sam(0); sam < NSam-1; sam++){
      if(Htag[indchain[sam]].Contains("sig")) continue;
      chain[indchain[sam]]->Project("histo1d", "onht", cutTrigger[icut]);
      trigger_rate[icut] += histo1d.Integral(0,2);
    } // Loop over samples
    cout<<nameTrigger[icut]<<": "<<RoundNumber(trigger_rate[icut],1)<<"\t";
  }
  cout<<endl;
  if(onlyRates) {cout<<endl; return;}
  //////////////////////////////////////////////////////////////////////

  for(int icut(0); icut < NCuts; icut++){
    totCut = "1.4e-2/19600*weight*(" + Cuts[icut] + ")";

    Title = "p_{T}^{#mu} > 15 GeV"+TitleCuts[icut];
    if(isEl) {
      Title.ReplaceAll("#mu","e");
    }
    cout<<endl<<"Doing cuts "<<totCut<<endl;
    for(int sam(0); sam < NSam; sam++){
      for(int his(1); his >= 0; his--) {
	Hname = Htag[indchain[sam]]; Hname += icut; Hname += his;
	hRate[his][icut][sam] = new TH2D(Hname,Title+"       L = 1.4 #times 10^{34} cm^{-2}s^{-1}",
					 nBinsHt, minHt, maxHt, nBinsMet, minMet, maxMet);
      }
      hRate[0][icut][sam]->SetXTitle(xTitle);
      hRate[0][icut][sam]->SetYTitle(yTitle);
      hRate[0][icut][sam]->SetZTitle(Htag[indchain[sam]]+" "+zTitle);
      hRate[1][icut][sam]->SetMarkerStyle(20);

      TString actualCuts(totCut);
      if(Htag[indchain[sam]].Contains("sig")) {
	hRate[0][icut][sam]->SetTitle(Title);
	actualCuts.ReplaceAll("1.4e-2/19600*weight*(","wl1ht200*(Max$(genmus_pt)>0&&");
	if(isEl) actualCuts.ReplaceAll("genmus", "genels");
	hRate[1][icut][sam]->SetMarkerSize(2.2);
      } else hRate[1][icut][sam]->SetMarkerSize(2.5);
      if(Htag[indchain[sam]]=="sig1025") hRate[0][icut][sam]->SetZTitle("T1tttt(1025,625) efficiency (%)");	
      if(Htag[indchain[sam]]=="sig1200") hRate[0][icut][sam]->SetZTitle("T1tttt(1200,800) efficiency (%)");	
      if(Htag[indchain[sam]]=="sig1500") hRate[0][icut][sam]->SetZTitle("T1tttt(1500,100) efficiency (%)");

      
      // Projecting chain
      if(Htag[indchain[sam]]!="Total") {
	chain[indchain[sam]]->Project(Hname, "abs(onmet):onht", actualCuts);
      } else hRate[0][icut][sam]->SetZTitle("HLT rate [QCD + tt + W#rightarrowl#nu] (Hz)");
    }

    // Finding yields for different cuts
    unsigned int iht(0); int itag(3);
    for(int htbin(1); htbin <= nBinsHt; htbin++){
      for(int metbin(1); metbin <= nBinsMet; metbin++){
    	float bkgrate(0), rate;
    	for(int sam(0); sam < NSam-1; sam++){
    	  rate = hRate[0][icut][sam]->Integral(htbin, nBinsHt+1, metbin, nBinsMet+1);
    	  if(!Htag[indchain[sam]].Contains("sig")) bkgrate += rate;
    	  else rate *= (100./(double)noriginal[indchain[sam]]);
    	  for(int his(0); his < 2; his++) 
    	    hRate[his][icut][sam]->SetBinContent(htbin, metbin, rate);
    	}
    	for(int his(0); his < 2; his++) 
    	  hRate[his][icut][NSam-1]->SetBinContent(htbin, metbin, bkgrate);
      } // Loop over MET bins
    } // Loop over HT bins

    // Plotting rates and efficiencies
    Pnamebase = "plots//rate2d_"+Tags[icut]+ ".eps";
    for(int sam(0); sam < NSam; sam++){
      Pname = Pnamebase; 
      if(Htag[indchain[sam]].Contains("sig")) {
    	Pname.ReplaceAll("rate2d", Htag[indchain[sam]]+"effi2d");
    	can.SetLogz(0);
      } else {
    	//hRate[0][icut][sam]->SetMinimum(hRate[0][icut][2]->GetMinimum());
    	//hRate[0][icut][sam]->SetMaximum(2*hRate[0][icut][NSam]->GetMaximum());
    	Pname.ReplaceAll("rate2d", Htag[indchain[sam]]+"rate2d");
    	can.SetLogz(1);
      }
      hRate[0][icut][sam]->Draw("colz");
      hRate[1][icut][sam]->Draw("text same");
      if(!applyCSV){
	box.DrawBox(600, 0, 652, 10);
	box.DrawBox(400, 70, 452, 80);
      } else box.DrawBox(400, 0, 452, 10);
      can.SaveAs(Pname);
    }
  } // Loop over cuts

  for(int icut(0); icut < NCuts; icut++){
    for(int tag(0); tag < NSam; tag++){
      for(int his(0); his < 2; his++)
  	if(hRate[his][icut][tag]) delete hRate[his][icut][tag];
    }
  }
}


void ReadChains(vector<TChain*> &chain, vector<int> &entries, TString folder, vector<TString> FileNames){

  bool isEl = folder.Contains("el15")?true:false;
  for(unsigned sam(0); sam < FileNames.size(); sam++){
    chain.push_back(new TChain("tree"));
    if(chain[sam]->Add(folder+FileNames[sam])==0){
      cout<<folder+FileNames[sam]<<" does not exist"<<endl;
      chain.clear();
      return;
    }
    if(FileNames[sam].Contains("T1tttt")) {
      int noriginal, nori_genmu0, nori_genel0;
      TChain tglobal("treeglobal");
      tglobal.Add(folder+FileNames[sam]);
      tglobal.SetBranchAddress("noriginal",&noriginal);
      tglobal.SetBranchAddress("nori_genel0",&nori_genel0);
      tglobal.SetBranchAddress("nori_genmu0",&nori_genmu0);
      tglobal.GetEntry(0);

      if(isEl) entries.push_back(nori_genel0);
      else entries.push_back(nori_genmu0);
    } else entries.push_back(0);
  }// Loop over samples
}


void singlerate_hlt(TString cuts="onht>600||onht>400&&onmet>70||onht>400&&Max$(bjets_csv*(bjets_pt>40))>0.7", 
		    TString muisocut = "1.2", TString eisocut = "0.8", 
		    TString folder="root/15-03-14/"){
  //Files
  vector<TString> filenames, Htag;
  vector<TChain*> chain;
  vector<int> noriginal;
  vector<int> indchain;
  vector<double> trigger_rate(2,0.);
  TString leptag[] = {"mu15", "el15"};
  TString Cuts[] = {"Max$(mus_pt*((mus_reliso)<"+muisocut+"))>15&&not_pu", 
		    "Max$(els_pt*(els_trackiso<0.4*"+eisocut+
		    "/0.8&&els_hcaliso<0.6*"+eisocut+"/0.8&&els_ecaliso<0.5*"+eisocut+"/0.8))>15&&not_pu"};
  //"not_pu"};
  //"Max$(els_pt*((els_trackiso+els_hcaliso+els_ecaliso)<"+eisocut+"))>15&&not_pu"};
  cout<<endl;
  for(int lep(0); lep<2; lep++){
    filenames.clear(); Htag.clear(); chain.clear(); noriginal.clear(); indchain.clear();

    filenames.push_back("/*QCD*");	        Htag.push_back("QCD");
    filenames.push_back("/*TT*");	        Htag.push_back("tt");
    filenames.push_back("/W*");                 Htag.push_back("Wjets");
    Htag.push_back("Total");
    ReadChains(chain, noriginal, folder+"/"+leptag[lep], filenames);
    if(chain.size()==0) return;
    indchain.push_back(0);
    indchain.push_back(1);
    indchain.push_back(2);
    indchain.push_back(Htag.size()-1);
    const int NSam = indchain.size();

    /////////////////// Calculating specific trigger rates  ////////////////
    TH1D histo1d("histo1d","",1,0,1);
    TString totCut("1.4e-2/19600*weight*(" + Cuts[lep] + "&&("+cuts+"))");
    for(int sam(0); sam < NSam-1; sam++){
      if(Htag[indchain[sam]].Contains("sig")) continue;
      chain[indchain[sam]]->Project("histo1d", "onht", totCut);
      trigger_rate[lep] += histo1d.Integral(0,2);
    } // Loop over samples
    cout<<leptag[lep]<<": "<<RoundNumber(trigger_rate[lep],1)<<" Hz\t";
  }
  cout<<"Total: "<< RoundNumber(trigger_rate[0]+trigger_rate[1],1)<<" Hz"<<endl<<endl;
}
