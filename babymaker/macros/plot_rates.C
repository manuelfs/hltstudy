// plot_rates: Plots HLT rates as a function of HT, MET

#define INT_ROOT
#include "styles.hpp"
#include "styles.cpp"
#include "plot_rates.hpp"

#include <stdexcept>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TH1D.h"


using namespace std;

void plot_rates(bool isEl=true){
  styles style("Standard"); style.setDefaultStyle();
  gStyle->SetPadTickX(1); // Tickmarks on the top
  gStyle->SetPadTickY(1); // Tickmarks on the right
  TCanvas can;

  //Files
  vector<TString> filenames, Htag;
  filenames.push_back("/*QCD*");	        Htag.push_back("QCD");
  filenames.push_back("/*TT*");		        Htag.push_back("tt");
  filenames.push_back("/W*");                   Htag.push_back("Wjets");
  vector<sample_class> folders;
  folders.push_back(sample_class("~/code/hltstudy/babymaker/macros/root/old_720/", "720pre8", 1)); 
  folders.push_back(sample_class("~/code/hltstudy/babymaker/macros/root/15-03-14/","740pre8 with multifit",2)); 
  folders.push_back(sample_class("~/code/hltstudy/babymaker/macros/root/15-04-01/", "740pre9 with multifit + method 3", 4)); 

  TString leptag(isEl?"el15":"mu15");
  vector<TChain*> chain;
  for(unsigned dir(0); dir < folders.size(); dir++){
    chain.push_back(new TChain("tree"));
    for(unsigned sam(0); sam < filenames.size(); sam++){
      TString files(folders[dir].files); files += leptag;
      files += filenames[sam];
      chain[dir]->Add(files);
    }
  } // Loop over folders

  // Legend
  double legX = 0.45, legY = 0.92, legSingle = 0.069;
  double legW = 0.12, legH = legSingle*folders.size();
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.062); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  // Histograms
  int nbins(500);
  TH1D *hisdenom(0);
  TString cuts="1", hname, totcut, Pname;
  vector<vector<TH1D> > histos;
  vector<var_class> vars;
  vars.push_back(var_class("onmet",0,300,"HLT PF MET threshold (GeV)"));
  vars.push_back(var_class("onht",200,1000,"HLT PF H_{T} threshold (GeV)"));

  for(unsigned var(0); var < vars.size(); var++){
    float maxhisto(-1);
    histos.push_back(vector<TH1D>());
    leg.Clear();
    for(unsigned dir(0); dir < folders.size(); dir++){
      hname = "histo"; hname += dir; hname += var;
      totcut = "1.4e-2/19600*weight*("+cuts+")";
      histos[var].push_back(TH1D(hname,"",nbins,vars[var].minx,vars[var].maxx));
      chain[dir]->Project(hname, vars[var].varname, totcut);
    } // Needed this intermediate loop to avoid crash
    for(unsigned dir(0); dir < folders.size(); dir++){
      histos[var][dir].SetLineColor(folders[dir].color);
      histos[var][dir].SetLineWidth(3);
      histos[var][dir].SetXTitle(vars[var].title);
      histos[var][dir].SetYTitle("HLT rate (Hz)");
      leg.AddEntry(&(histos[var][dir]), folders[dir].label, "l");

      for(int bin(nbins); bin>=1; bin--)
      	histos[var][dir].SetBinContent(bin, histos[var][dir].GetBinContent(bin)+
      				       histos[var][dir].GetBinContent(bin+1));
      if(maxhisto<histos[var][dir].GetMaximum()) maxhisto = histos[var][dir].GetMaximum();
      if(dir==0) histos[var][dir].Draw("c");
      else histos[var][dir].Draw("c same");
    } // Loop over folders
    histos[var][0].SetMaximum(1.2*maxhisto);
    histos[var][0].SetMinimum(0);
    leg.Draw();
    Pname = "plots/rate_"+vars[var].varname+"_"+leptag+".eps";
    can.SetLogy(0);
    can.SaveAs(Pname);
    can.SetLogy(1);
    histos[var][0].SetMinimum(0.001);
    Pname.ReplaceAll("rate","lograte");
    can.SaveAs(Pname);

    // Plotting ratios
    for(unsigned dir(0); dir < folders.size(); dir++){
      histos[var][dir].SetYTitle("Ratio with respect to 720pre8");
      if(dir==0) hisdenom = static_cast<TH1D*>(histos[var][dir].Clone("hisdenom"));
      histos[var][dir].Divide(hisdenom);
      if(dir==0) histos[var][dir].Draw("c");
      else histos[var][dir].Draw("c same");
    } // Loop over folders
    histos[var][0].SetMaximum(1.5);
    histos[var][0].SetMinimum(0);
    leg.Draw();
    Pname = "plots/ratio_"+vars[var].varname+"_"+leptag+".eps";
    can.SetLogy(0);
    can.SaveAs(Pname);
    can.SetLogy(1);
    histos[var][0].SetMinimum(0.001);
    Pname.ReplaceAll("ratio","logratio");
    can.SaveAs(Pname);
  } // Loop over variables


  for(unsigned dir(0); dir < folders.size(); dir++)
    chain[dir]->Delete();

}

var_class::var_class(TString ivarname, float iminx, float imaxx, TString ititle){
  varname = ivarname; minx = iminx; maxx = imaxx; title = ititle;
}

sample_class::sample_class(TString ifiles, TString ilabel, int icolor){
  files = ifiles; label = ilabel; color = icolor;
}


