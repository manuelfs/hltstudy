// -*- C++ -*-
#ifndef hlt_class_H
#define hlt_class_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector>
#include <unistd.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#define PARANOIA

using namespace std; 
class hlt_class {
private: 
protected: 
	unsigned int index;
	float calo_ht_;
	TBranch *calo_ht_branch;
	bool calo_ht_isLoaded;
	float calo_mht_eta_;
	TBranch *calo_mht_eta_branch;
	bool calo_mht_eta_isLoaded;
	float calo_mht_phi_;
	TBranch *calo_mht_phi_branch;
	bool calo_mht_phi_isLoaded;
	float calo_mht_pt_;
	TBranch *calo_mht_pt_branch;
	bool calo_mht_pt_isLoaded;
	float gen_ht_;
	TBranch *gen_ht_branch;
	bool gen_ht_isLoaded;
	float gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float gen_metcalo_;
	TBranch *gen_metcalo_branch;
	bool gen_metcalo_isLoaded;
	float gen_metcalononprompt_;
	TBranch *gen_metcalononprompt_branch;
	bool gen_metcalononprompt_isLoaded;
	float met_eta_;
	TBranch *met_eta_branch;
	bool met_eta_isLoaded;
	float met_phi_;
	TBranch *met_phi_branch;
	bool met_phi_isLoaded;
	float met_pt_;
	TBranch *met_pt_branch;
	bool met_pt_isLoaded;
	float pf_ht_;
	TBranch *pf_ht_branch;
	bool pf_ht_isLoaded;
	float pf_mht_eta_;
	TBranch *pf_mht_eta_branch;
	bool pf_mht_eta_isLoaded;
	float pf_mht_phi_;
	TBranch *pf_mht_phi_branch;
	bool pf_mht_phi_isLoaded;
	float pf_mht_pt_;
	TBranch *pf_mht_pt_branch;
	bool pf_mht_pt_isLoaded;
	float wl1ht200_;
	TBranch *wl1ht200_branch;
	bool wl1ht200_isLoaded;
	vector<float> bjets_csv_;
	TBranch *bjets_csv_branch;
	bool bjets_csv_isLoaded;
	vector<float> bjets_eta_;
	TBranch *bjets_eta_branch;
	bool bjets_eta_isLoaded;
	vector<float> bjets_phi_;
	TBranch *bjets_phi_branch;
	bool bjets_phi_isLoaded;
	vector<float> bjets_pt_;
	TBranch *bjets_pt_branch;
	bool bjets_pt_isLoaded;
	vector<float> ele_eta_;
	TBranch *ele_eta_branch;
	bool ele_eta_isLoaded;
	vector<float> ele_phi_;
	TBranch *ele_phi_branch;
	bool ele_phi_isLoaded;
	vector<float> ele_pt_;
	TBranch *ele_pt_branch;
	bool ele_pt_isLoaded;
	vector<float> els_clustershape_;
	TBranch *els_clustershape_branch;
	bool els_clustershape_isLoaded;
	vector<float> els_deta_;
	TBranch *els_deta_branch;
	bool els_deta_isLoaded;
	vector<float> els_dphi_;
	TBranch *els_dphi_branch;
	bool els_dphi_isLoaded;
	vector<float> els_ecal_iso_;
	TBranch *els_ecal_iso_branch;
	bool els_ecal_iso_isLoaded;
	vector<float> els_eminusp_;
	TBranch *els_eminusp_branch;
	bool els_eminusp_isLoaded;
	vector<float> els_eta_;
	TBranch *els_eta_branch;
	bool els_eta_isLoaded;
	vector<float> els_hcal_iso_;
	TBranch *els_hcal_iso_branch;
	bool els_hcal_iso_isLoaded;
	vector<float> els_he_;
	TBranch *els_he_branch;
	bool els_he_isLoaded;
	vector<float> els_phi_;
	TBranch *els_phi_branch;
	bool els_phi_isLoaded;
	vector<float> els_pt_;
	TBranch *els_pt_branch;
	bool els_pt_isLoaded;
	vector<float> els_track_iso_;
	TBranch *els_track_iso_branch;
	bool els_track_iso_isLoaded;
	vector<float> genels_eta_;
	TBranch *genels_eta_branch;
	bool genels_eta_isLoaded;
	vector<float> genels_phi_;
	TBranch *genels_phi_branch;
	bool genels_phi_isLoaded;
	vector<float> genels_pt_;
	TBranch *genels_pt_branch;
	bool genels_pt_isLoaded;
	vector<float> genjets_eta_;
	TBranch *genjets_eta_branch;
	bool genjets_eta_isLoaded;
	vector<float> genjets_phi_;
	TBranch *genjets_phi_branch;
	bool genjets_phi_isLoaded;
	vector<float> genjets_pt_;
	TBranch *genjets_pt_branch;
	bool genjets_pt_isLoaded;
	vector<float> genmus_eta_;
	TBranch *genmus_eta_branch;
	bool genmus_eta_isLoaded;
	vector<float> genmus_phi_;
	TBranch *genmus_phi_branch;
	bool genmus_phi_isLoaded;
	vector<float> genmus_pt_;
	TBranch *genmus_pt_branch;
	bool genmus_pt_isLoaded;
	vector<float> gentop_eta_;
	TBranch *gentop_eta_branch;
	bool gentop_eta_isLoaded;
	vector<float> gentop_phi_;
	TBranch *gentop_phi_branch;
	bool gentop_phi_isLoaded;
	vector<float> gentop_pt_;
	TBranch *gentop_pt_branch;
	bool gentop_pt_isLoaded;
	vector<float> mus_eta_;
	TBranch *mus_eta_branch;
	bool mus_eta_isLoaded;
	vector<float> mus_iso_;
	TBranch *mus_iso_branch;
	bool mus_iso_isLoaded;
	vector<float> mus_phi_;
	TBranch *mus_phi_branch;
	bool mus_phi_isLoaded;
	vector<float> mus_pt_;
	TBranch *mus_pt_branch;
	bool mus_pt_isLoaded;
	vector<float> pfjets_eta_;
	TBranch *pfjets_eta_branch;
	bool pfjets_eta_isLoaded;
	vector<float> pfjets_phi_;
	TBranch *pfjets_phi_branch;
	bool pfjets_phi_isLoaded;
	vector<float> pfjets_pt_;
	TBranch *pfjets_pt_branch;
	bool pfjets_pt_isLoaded;
	vector<float> reco_els_eta_;
	TBranch *reco_els_eta_branch;
	bool reco_els_eta_isLoaded;
	vector<float> reco_els_phi_;
	TBranch *reco_els_phi_branch;
	bool reco_els_phi_isLoaded;
	vector<float> reco_els_pt_;
	TBranch *reco_els_pt_branch;
	bool reco_els_pt_isLoaded;
	vector<float> reco_genjet_eta_;
	TBranch *reco_genjet_eta_branch;
	bool reco_genjet_eta_isLoaded;
	vector<float> reco_genjet_phi_;
	TBranch *reco_genjet_phi_branch;
	bool reco_genjet_phi_isLoaded;
	vector<float> reco_genjet_pt_;
	TBranch *reco_genjet_pt_branch;
	bool reco_genjet_pt_isLoaded;
	vector<float> reco_jet_eta_;
	TBranch *reco_jet_eta_branch;
	bool reco_jet_eta_isLoaded;
	vector<float> reco_jet_phi_;
	TBranch *reco_jet_phi_branch;
	bool reco_jet_phi_isLoaded;
	vector<float> reco_jet_pt_;
	TBranch *reco_jet_pt_branch;
	bool reco_jet_pt_isLoaded;
	vector<float> reco_mus_eta_;
	TBranch *reco_mus_eta_branch;
	bool reco_mus_eta_isLoaded;
	vector<float> reco_mus_phi_;
	TBranch *reco_mus_phi_branch;
	bool reco_mus_phi_isLoaded;
	vector<float> reco_mus_pt_;
	TBranch *reco_mus_pt_branch;
	bool reco_mus_pt_isLoaded;
	vector<int> genels_ggmom_id_;
	TBranch *genels_ggmom_id_branch;
	bool genels_ggmom_id_isLoaded;
	vector<int> genels_gmom_id_;
	TBranch *genels_gmom_id_branch;
	bool genels_gmom_id_isLoaded;
	vector<int> genels_mom_id_;
	TBranch *genels_mom_id_branch;
	bool genels_mom_id_isLoaded;
	vector<int> genmus_ggmom_id_;
	TBranch *genmus_ggmom_id_branch;
	bool genmus_ggmom_id_isLoaded;
	vector<int> genmus_gmom_id_;
	TBranch *genmus_gmom_id_branch;
	bool genmus_gmom_id_isLoaded;
	vector<int> genmus_mom_id_;
	TBranch *genmus_mom_id_branch;
	bool genmus_mom_id_isLoaded;
	vector<int> gentop_id_;
	TBranch *gentop_id_branch;
	bool gentop_id_isLoaded;
	vector<int> gentop_mom_id_;
	TBranch *gentop_mom_id_branch;
	bool gentop_mom_id_isLoaded;
	unsigned int ngenlep_;
	TBranch *ngenlep_branch;
	bool ngenlep_isLoaded;
public: 
void Init(TTree *tree) {
  tree->SetMakeClass(1);
	calo_ht_branch = 0;
	if (tree->GetAlias("calo_ht") != 0) {
		calo_ht_branch = tree->GetBranch(tree->GetAlias("calo_ht"));
		if (calo_ht_branch) {calo_ht_branch->SetAddress(&calo_ht_);}
	}
	calo_mht_eta_branch = 0;
	if (tree->GetAlias("calo_mht_eta") != 0) {
		calo_mht_eta_branch = tree->GetBranch(tree->GetAlias("calo_mht_eta"));
		if (calo_mht_eta_branch) {calo_mht_eta_branch->SetAddress(&calo_mht_eta_);}
	}
	calo_mht_phi_branch = 0;
	if (tree->GetAlias("calo_mht_phi") != 0) {
		calo_mht_phi_branch = tree->GetBranch(tree->GetAlias("calo_mht_phi"));
		if (calo_mht_phi_branch) {calo_mht_phi_branch->SetAddress(&calo_mht_phi_);}
	}
	calo_mht_pt_branch = 0;
	if (tree->GetAlias("calo_mht_pt") != 0) {
		calo_mht_pt_branch = tree->GetBranch(tree->GetAlias("calo_mht_pt"));
		if (calo_mht_pt_branch) {calo_mht_pt_branch->SetAddress(&calo_mht_pt_);}
	}
	gen_ht_branch = 0;
	if (tree->GetAlias("gen_ht") != 0) {
		gen_ht_branch = tree->GetBranch(tree->GetAlias("gen_ht"));
		if (gen_ht_branch) {gen_ht_branch->SetAddress(&gen_ht_);}
	}
	gen_met_branch = 0;
	if (tree->GetAlias("gen_met") != 0) {
		gen_met_branch = tree->GetBranch(tree->GetAlias("gen_met"));
		if (gen_met_branch) {gen_met_branch->SetAddress(&gen_met_);}
	}
	gen_metcalo_branch = 0;
	if (tree->GetAlias("gen_metcalo") != 0) {
		gen_metcalo_branch = tree->GetBranch(tree->GetAlias("gen_metcalo"));
		if (gen_metcalo_branch) {gen_metcalo_branch->SetAddress(&gen_metcalo_);}
	}
	gen_metcalononprompt_branch = 0;
	if (tree->GetAlias("gen_metcalononprompt") != 0) {
		gen_metcalononprompt_branch = tree->GetBranch(tree->GetAlias("gen_metcalononprompt"));
		if (gen_metcalononprompt_branch) {gen_metcalononprompt_branch->SetAddress(&gen_metcalononprompt_);}
	}
	met_eta_branch = 0;
	if (tree->GetAlias("met_eta") != 0) {
		met_eta_branch = tree->GetBranch(tree->GetAlias("met_eta"));
		if (met_eta_branch) {met_eta_branch->SetAddress(&met_eta_);}
	}
	met_phi_branch = 0;
	if (tree->GetAlias("met_phi") != 0) {
		met_phi_branch = tree->GetBranch(tree->GetAlias("met_phi"));
		if (met_phi_branch) {met_phi_branch->SetAddress(&met_phi_);}
	}
	met_pt_branch = 0;
	if (tree->GetAlias("met_pt") != 0) {
		met_pt_branch = tree->GetBranch(tree->GetAlias("met_pt"));
		if (met_pt_branch) {met_pt_branch->SetAddress(&met_pt_);}
	}
	pf_ht_branch = 0;
	if (tree->GetAlias("pf_ht") != 0) {
		pf_ht_branch = tree->GetBranch(tree->GetAlias("pf_ht"));
		if (pf_ht_branch) {pf_ht_branch->SetAddress(&pf_ht_);}
	}
	pf_mht_eta_branch = 0;
	if (tree->GetAlias("pf_mht_eta") != 0) {
		pf_mht_eta_branch = tree->GetBranch(tree->GetAlias("pf_mht_eta"));
		if (pf_mht_eta_branch) {pf_mht_eta_branch->SetAddress(&pf_mht_eta_);}
	}
	pf_mht_phi_branch = 0;
	if (tree->GetAlias("pf_mht_phi") != 0) {
		pf_mht_phi_branch = tree->GetBranch(tree->GetAlias("pf_mht_phi"));
		if (pf_mht_phi_branch) {pf_mht_phi_branch->SetAddress(&pf_mht_phi_);}
	}
	pf_mht_pt_branch = 0;
	if (tree->GetAlias("pf_mht_pt") != 0) {
		pf_mht_pt_branch = tree->GetBranch(tree->GetAlias("pf_mht_pt"));
		if (pf_mht_pt_branch) {pf_mht_pt_branch->SetAddress(&pf_mht_pt_);}
	}
	wl1ht200_branch = 0;
	if (tree->GetAlias("wl1ht200") != 0) {
		wl1ht200_branch = tree->GetBranch(tree->GetAlias("wl1ht200"));
		if (wl1ht200_branch) {wl1ht200_branch->SetAddress(&wl1ht200_);}
	}
	bjets_csv_branch = 0;
	if (tree->GetAlias("bjets_csv") != 0) {
		bjets_csv_branch = tree->GetBranch(tree->GetAlias("bjets_csv"));
		if (bjets_csv_branch) {bjets_csv_branch->SetAddress(&bjets_csv_);}
	}
	bjets_eta_branch = 0;
	if (tree->GetAlias("bjets_eta") != 0) {
		bjets_eta_branch = tree->GetBranch(tree->GetAlias("bjets_eta"));
		if (bjets_eta_branch) {bjets_eta_branch->SetAddress(&bjets_eta_);}
	}
	bjets_phi_branch = 0;
	if (tree->GetAlias("bjets_phi") != 0) {
		bjets_phi_branch = tree->GetBranch(tree->GetAlias("bjets_phi"));
		if (bjets_phi_branch) {bjets_phi_branch->SetAddress(&bjets_phi_);}
	}
	bjets_pt_branch = 0;
	if (tree->GetAlias("bjets_pt") != 0) {
		bjets_pt_branch = tree->GetBranch(tree->GetAlias("bjets_pt"));
		if (bjets_pt_branch) {bjets_pt_branch->SetAddress(&bjets_pt_);}
	}
	ele_eta_branch = 0;
	if (tree->GetAlias("ele_eta") != 0) {
		ele_eta_branch = tree->GetBranch(tree->GetAlias("ele_eta"));
		if (ele_eta_branch) {ele_eta_branch->SetAddress(&ele_eta_);}
	}
	ele_phi_branch = 0;
	if (tree->GetAlias("ele_phi") != 0) {
		ele_phi_branch = tree->GetBranch(tree->GetAlias("ele_phi"));
		if (ele_phi_branch) {ele_phi_branch->SetAddress(&ele_phi_);}
	}
	ele_pt_branch = 0;
	if (tree->GetAlias("ele_pt") != 0) {
		ele_pt_branch = tree->GetBranch(tree->GetAlias("ele_pt"));
		if (ele_pt_branch) {ele_pt_branch->SetAddress(&ele_pt_);}
	}
	els_clustershape_branch = 0;
	if (tree->GetAlias("els_clustershape") != 0) {
		els_clustershape_branch = tree->GetBranch(tree->GetAlias("els_clustershape"));
		if (els_clustershape_branch) {els_clustershape_branch->SetAddress(&els_clustershape_);}
	}
	els_deta_branch = 0;
	if (tree->GetAlias("els_deta") != 0) {
		els_deta_branch = tree->GetBranch(tree->GetAlias("els_deta"));
		if (els_deta_branch) {els_deta_branch->SetAddress(&els_deta_);}
	}
	els_dphi_branch = 0;
	if (tree->GetAlias("els_dphi") != 0) {
		els_dphi_branch = tree->GetBranch(tree->GetAlias("els_dphi"));
		if (els_dphi_branch) {els_dphi_branch->SetAddress(&els_dphi_);}
	}
	els_ecal_iso_branch = 0;
	if (tree->GetAlias("els_ecal_iso") != 0) {
		els_ecal_iso_branch = tree->GetBranch(tree->GetAlias("els_ecal_iso"));
		if (els_ecal_iso_branch) {els_ecal_iso_branch->SetAddress(&els_ecal_iso_);}
	}
	els_eminusp_branch = 0;
	if (tree->GetAlias("els_eminusp") != 0) {
		els_eminusp_branch = tree->GetBranch(tree->GetAlias("els_eminusp"));
		if (els_eminusp_branch) {els_eminusp_branch->SetAddress(&els_eminusp_);}
	}
	els_eta_branch = 0;
	if (tree->GetAlias("els_eta") != 0) {
		els_eta_branch = tree->GetBranch(tree->GetAlias("els_eta"));
		if (els_eta_branch) {els_eta_branch->SetAddress(&els_eta_);}
	}
	els_hcal_iso_branch = 0;
	if (tree->GetAlias("els_hcal_iso") != 0) {
		els_hcal_iso_branch = tree->GetBranch(tree->GetAlias("els_hcal_iso"));
		if (els_hcal_iso_branch) {els_hcal_iso_branch->SetAddress(&els_hcal_iso_);}
	}
	els_he_branch = 0;
	if (tree->GetAlias("els_he") != 0) {
		els_he_branch = tree->GetBranch(tree->GetAlias("els_he"));
		if (els_he_branch) {els_he_branch->SetAddress(&els_he_);}
	}
	els_phi_branch = 0;
	if (tree->GetAlias("els_phi") != 0) {
		els_phi_branch = tree->GetBranch(tree->GetAlias("els_phi"));
		if (els_phi_branch) {els_phi_branch->SetAddress(&els_phi_);}
	}
	els_pt_branch = 0;
	if (tree->GetAlias("els_pt") != 0) {
		els_pt_branch = tree->GetBranch(tree->GetAlias("els_pt"));
		if (els_pt_branch) {els_pt_branch->SetAddress(&els_pt_);}
	}
	els_track_iso_branch = 0;
	if (tree->GetAlias("els_track_iso") != 0) {
		els_track_iso_branch = tree->GetBranch(tree->GetAlias("els_track_iso"));
		if (els_track_iso_branch) {els_track_iso_branch->SetAddress(&els_track_iso_);}
	}
	genels_eta_branch = 0;
	if (tree->GetAlias("genels_eta") != 0) {
		genels_eta_branch = tree->GetBranch(tree->GetAlias("genels_eta"));
		if (genels_eta_branch) {genels_eta_branch->SetAddress(&genels_eta_);}
	}
	genels_phi_branch = 0;
	if (tree->GetAlias("genels_phi") != 0) {
		genels_phi_branch = tree->GetBranch(tree->GetAlias("genels_phi"));
		if (genels_phi_branch) {genels_phi_branch->SetAddress(&genels_phi_);}
	}
	genels_pt_branch = 0;
	if (tree->GetAlias("genels_pt") != 0) {
		genels_pt_branch = tree->GetBranch(tree->GetAlias("genels_pt"));
		if (genels_pt_branch) {genels_pt_branch->SetAddress(&genels_pt_);}
	}
	genjets_eta_branch = 0;
	if (tree->GetAlias("genjets_eta") != 0) {
		genjets_eta_branch = tree->GetBranch(tree->GetAlias("genjets_eta"));
		if (genjets_eta_branch) {genjets_eta_branch->SetAddress(&genjets_eta_);}
	}
	genjets_phi_branch = 0;
	if (tree->GetAlias("genjets_phi") != 0) {
		genjets_phi_branch = tree->GetBranch(tree->GetAlias("genjets_phi"));
		if (genjets_phi_branch) {genjets_phi_branch->SetAddress(&genjets_phi_);}
	}
	genjets_pt_branch = 0;
	if (tree->GetAlias("genjets_pt") != 0) {
		genjets_pt_branch = tree->GetBranch(tree->GetAlias("genjets_pt"));
		if (genjets_pt_branch) {genjets_pt_branch->SetAddress(&genjets_pt_);}
	}
	genmus_eta_branch = 0;
	if (tree->GetAlias("genmus_eta") != 0) {
		genmus_eta_branch = tree->GetBranch(tree->GetAlias("genmus_eta"));
		if (genmus_eta_branch) {genmus_eta_branch->SetAddress(&genmus_eta_);}
	}
	genmus_phi_branch = 0;
	if (tree->GetAlias("genmus_phi") != 0) {
		genmus_phi_branch = tree->GetBranch(tree->GetAlias("genmus_phi"));
		if (genmus_phi_branch) {genmus_phi_branch->SetAddress(&genmus_phi_);}
	}
	genmus_pt_branch = 0;
	if (tree->GetAlias("genmus_pt") != 0) {
		genmus_pt_branch = tree->GetBranch(tree->GetAlias("genmus_pt"));
		if (genmus_pt_branch) {genmus_pt_branch->SetAddress(&genmus_pt_);}
	}
	gentop_eta_branch = 0;
	if (tree->GetAlias("gentop_eta") != 0) {
		gentop_eta_branch = tree->GetBranch(tree->GetAlias("gentop_eta"));
		if (gentop_eta_branch) {gentop_eta_branch->SetAddress(&gentop_eta_);}
	}
	gentop_phi_branch = 0;
	if (tree->GetAlias("gentop_phi") != 0) {
		gentop_phi_branch = tree->GetBranch(tree->GetAlias("gentop_phi"));
		if (gentop_phi_branch) {gentop_phi_branch->SetAddress(&gentop_phi_);}
	}
	gentop_pt_branch = 0;
	if (tree->GetAlias("gentop_pt") != 0) {
		gentop_pt_branch = tree->GetBranch(tree->GetAlias("gentop_pt"));
		if (gentop_pt_branch) {gentop_pt_branch->SetAddress(&gentop_pt_);}
	}
	mus_eta_branch = 0;
	if (tree->GetAlias("mus_eta") != 0) {
		mus_eta_branch = tree->GetBranch(tree->GetAlias("mus_eta"));
		if (mus_eta_branch) {mus_eta_branch->SetAddress(&mus_eta_);}
	}
	mus_iso_branch = 0;
	if (tree->GetAlias("mus_iso") != 0) {
		mus_iso_branch = tree->GetBranch(tree->GetAlias("mus_iso"));
		if (mus_iso_branch) {mus_iso_branch->SetAddress(&mus_iso_);}
	}
	mus_phi_branch = 0;
	if (tree->GetAlias("mus_phi") != 0) {
		mus_phi_branch = tree->GetBranch(tree->GetAlias("mus_phi"));
		if (mus_phi_branch) {mus_phi_branch->SetAddress(&mus_phi_);}
	}
	mus_pt_branch = 0;
	if (tree->GetAlias("mus_pt") != 0) {
		mus_pt_branch = tree->GetBranch(tree->GetAlias("mus_pt"));
		if (mus_pt_branch) {mus_pt_branch->SetAddress(&mus_pt_);}
	}
	pfjets_eta_branch = 0;
	if (tree->GetAlias("pfjets_eta") != 0) {
		pfjets_eta_branch = tree->GetBranch(tree->GetAlias("pfjets_eta"));
		if (pfjets_eta_branch) {pfjets_eta_branch->SetAddress(&pfjets_eta_);}
	}
	pfjets_phi_branch = 0;
	if (tree->GetAlias("pfjets_phi") != 0) {
		pfjets_phi_branch = tree->GetBranch(tree->GetAlias("pfjets_phi"));
		if (pfjets_phi_branch) {pfjets_phi_branch->SetAddress(&pfjets_phi_);}
	}
	pfjets_pt_branch = 0;
	if (tree->GetAlias("pfjets_pt") != 0) {
		pfjets_pt_branch = tree->GetBranch(tree->GetAlias("pfjets_pt"));
		if (pfjets_pt_branch) {pfjets_pt_branch->SetAddress(&pfjets_pt_);}
	}
	reco_els_eta_branch = 0;
	if (tree->GetAlias("reco_els_eta") != 0) {
		reco_els_eta_branch = tree->GetBranch(tree->GetAlias("reco_els_eta"));
		if (reco_els_eta_branch) {reco_els_eta_branch->SetAddress(&reco_els_eta_);}
	}
	reco_els_phi_branch = 0;
	if (tree->GetAlias("reco_els_phi") != 0) {
		reco_els_phi_branch = tree->GetBranch(tree->GetAlias("reco_els_phi"));
		if (reco_els_phi_branch) {reco_els_phi_branch->SetAddress(&reco_els_phi_);}
	}
	reco_els_pt_branch = 0;
	if (tree->GetAlias("reco_els_pt") != 0) {
		reco_els_pt_branch = tree->GetBranch(tree->GetAlias("reco_els_pt"));
		if (reco_els_pt_branch) {reco_els_pt_branch->SetAddress(&reco_els_pt_);}
	}
	reco_genjet_eta_branch = 0;
	if (tree->GetAlias("reco_genjet_eta") != 0) {
		reco_genjet_eta_branch = tree->GetBranch(tree->GetAlias("reco_genjet_eta"));
		if (reco_genjet_eta_branch) {reco_genjet_eta_branch->SetAddress(&reco_genjet_eta_);}
	}
	reco_genjet_phi_branch = 0;
	if (tree->GetAlias("reco_genjet_phi") != 0) {
		reco_genjet_phi_branch = tree->GetBranch(tree->GetAlias("reco_genjet_phi"));
		if (reco_genjet_phi_branch) {reco_genjet_phi_branch->SetAddress(&reco_genjet_phi_);}
	}
	reco_genjet_pt_branch = 0;
	if (tree->GetAlias("reco_genjet_pt") != 0) {
		reco_genjet_pt_branch = tree->GetBranch(tree->GetAlias("reco_genjet_pt"));
		if (reco_genjet_pt_branch) {reco_genjet_pt_branch->SetAddress(&reco_genjet_pt_);}
	}
	reco_jet_eta_branch = 0;
	if (tree->GetAlias("reco_jet_eta") != 0) {
		reco_jet_eta_branch = tree->GetBranch(tree->GetAlias("reco_jet_eta"));
		if (reco_jet_eta_branch) {reco_jet_eta_branch->SetAddress(&reco_jet_eta_);}
	}
	reco_jet_phi_branch = 0;
	if (tree->GetAlias("reco_jet_phi") != 0) {
		reco_jet_phi_branch = tree->GetBranch(tree->GetAlias("reco_jet_phi"));
		if (reco_jet_phi_branch) {reco_jet_phi_branch->SetAddress(&reco_jet_phi_);}
	}
	reco_jet_pt_branch = 0;
	if (tree->GetAlias("reco_jet_pt") != 0) {
		reco_jet_pt_branch = tree->GetBranch(tree->GetAlias("reco_jet_pt"));
		if (reco_jet_pt_branch) {reco_jet_pt_branch->SetAddress(&reco_jet_pt_);}
	}
	reco_mus_eta_branch = 0;
	if (tree->GetAlias("reco_mus_eta") != 0) {
		reco_mus_eta_branch = tree->GetBranch(tree->GetAlias("reco_mus_eta"));
		if (reco_mus_eta_branch) {reco_mus_eta_branch->SetAddress(&reco_mus_eta_);}
	}
	reco_mus_phi_branch = 0;
	if (tree->GetAlias("reco_mus_phi") != 0) {
		reco_mus_phi_branch = tree->GetBranch(tree->GetAlias("reco_mus_phi"));
		if (reco_mus_phi_branch) {reco_mus_phi_branch->SetAddress(&reco_mus_phi_);}
	}
	reco_mus_pt_branch = 0;
	if (tree->GetAlias("reco_mus_pt") != 0) {
		reco_mus_pt_branch = tree->GetBranch(tree->GetAlias("reco_mus_pt"));
		if (reco_mus_pt_branch) {reco_mus_pt_branch->SetAddress(&reco_mus_pt_);}
	}
	genels_ggmom_id_branch = 0;
	if (tree->GetAlias("genels_ggmom_id") != 0) {
		genels_ggmom_id_branch = tree->GetBranch(tree->GetAlias("genels_ggmom_id"));
		if (genels_ggmom_id_branch) {genels_ggmom_id_branch->SetAddress(&genels_ggmom_id_);}
	}
	genels_gmom_id_branch = 0;
	if (tree->GetAlias("genels_gmom_id") != 0) {
		genels_gmom_id_branch = tree->GetBranch(tree->GetAlias("genels_gmom_id"));
		if (genels_gmom_id_branch) {genels_gmom_id_branch->SetAddress(&genels_gmom_id_);}
	}
	genels_mom_id_branch = 0;
	if (tree->GetAlias("genels_mom_id") != 0) {
		genels_mom_id_branch = tree->GetBranch(tree->GetAlias("genels_mom_id"));
		if (genels_mom_id_branch) {genels_mom_id_branch->SetAddress(&genels_mom_id_);}
	}
	genmus_ggmom_id_branch = 0;
	if (tree->GetAlias("genmus_ggmom_id") != 0) {
		genmus_ggmom_id_branch = tree->GetBranch(tree->GetAlias("genmus_ggmom_id"));
		if (genmus_ggmom_id_branch) {genmus_ggmom_id_branch->SetAddress(&genmus_ggmom_id_);}
	}
	genmus_gmom_id_branch = 0;
	if (tree->GetAlias("genmus_gmom_id") != 0) {
		genmus_gmom_id_branch = tree->GetBranch(tree->GetAlias("genmus_gmom_id"));
		if (genmus_gmom_id_branch) {genmus_gmom_id_branch->SetAddress(&genmus_gmom_id_);}
	}
	genmus_mom_id_branch = 0;
	if (tree->GetAlias("genmus_mom_id") != 0) {
		genmus_mom_id_branch = tree->GetBranch(tree->GetAlias("genmus_mom_id"));
		if (genmus_mom_id_branch) {genmus_mom_id_branch->SetAddress(&genmus_mom_id_);}
	}
	gentop_id_branch = 0;
	if (tree->GetAlias("gentop_id") != 0) {
		gentop_id_branch = tree->GetBranch(tree->GetAlias("gentop_id"));
		if (gentop_id_branch) {gentop_id_branch->SetAddress(&gentop_id_);}
	}
	gentop_mom_id_branch = 0;
	if (tree->GetAlias("gentop_mom_id") != 0) {
		gentop_mom_id_branch = tree->GetBranch(tree->GetAlias("gentop_mom_id"));
		if (gentop_mom_id_branch) {gentop_mom_id_branch->SetAddress(&gentop_mom_id_);}
	}
	ngenlep_branch = 0;
	if (tree->GetAlias("ngenlep") != 0) {
		ngenlep_branch = tree->GetBranch(tree->GetAlias("ngenlep"));
		if (ngenlep_branch) {ngenlep_branch->SetAddress(&ngenlep_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		calo_ht_isLoaded = false;
		calo_mht_eta_isLoaded = false;
		calo_mht_phi_isLoaded = false;
		calo_mht_pt_isLoaded = false;
		gen_ht_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metcalo_isLoaded = false;
		gen_metcalononprompt_isLoaded = false;
		met_eta_isLoaded = false;
		met_phi_isLoaded = false;
		met_pt_isLoaded = false;
		pf_ht_isLoaded = false;
		pf_mht_eta_isLoaded = false;
		pf_mht_phi_isLoaded = false;
		pf_mht_pt_isLoaded = false;
		wl1ht200_isLoaded = false;
		bjets_csv_isLoaded = false;
		bjets_eta_isLoaded = false;
		bjets_phi_isLoaded = false;
		bjets_pt_isLoaded = false;
		ele_eta_isLoaded = false;
		ele_phi_isLoaded = false;
		ele_pt_isLoaded = false;
		els_clustershape_isLoaded = false;
		els_deta_isLoaded = false;
		els_dphi_isLoaded = false;
		els_ecal_iso_isLoaded = false;
		els_eminusp_isLoaded = false;
		els_eta_isLoaded = false;
		els_hcal_iso_isLoaded = false;
		els_he_isLoaded = false;
		els_phi_isLoaded = false;
		els_pt_isLoaded = false;
		els_track_iso_isLoaded = false;
		genels_eta_isLoaded = false;
		genels_phi_isLoaded = false;
		genels_pt_isLoaded = false;
		genjets_eta_isLoaded = false;
		genjets_phi_isLoaded = false;
		genjets_pt_isLoaded = false;
		genmus_eta_isLoaded = false;
		genmus_phi_isLoaded = false;
		genmus_pt_isLoaded = false;
		gentop_eta_isLoaded = false;
		gentop_phi_isLoaded = false;
		gentop_pt_isLoaded = false;
		mus_eta_isLoaded = false;
		mus_iso_isLoaded = false;
		mus_phi_isLoaded = false;
		mus_pt_isLoaded = false;
		pfjets_eta_isLoaded = false;
		pfjets_phi_isLoaded = false;
		pfjets_pt_isLoaded = false;
		reco_els_eta_isLoaded = false;
		reco_els_phi_isLoaded = false;
		reco_els_pt_isLoaded = false;
		reco_genjet_eta_isLoaded = false;
		reco_genjet_phi_isLoaded = false;
		reco_genjet_pt_isLoaded = false;
		reco_jet_eta_isLoaded = false;
		reco_jet_phi_isLoaded = false;
		reco_jet_pt_isLoaded = false;
		reco_mus_eta_isLoaded = false;
		reco_mus_phi_isLoaded = false;
		reco_mus_pt_isLoaded = false;
		genels_ggmom_id_isLoaded = false;
		genels_gmom_id_isLoaded = false;
		genels_mom_id_isLoaded = false;
		genmus_ggmom_id_isLoaded = false;
		genmus_gmom_id_isLoaded = false;
		genmus_mom_id_isLoaded = false;
		gentop_id_isLoaded = false;
		gentop_mom_id_isLoaded = false;
		ngenlep_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (calo_ht_branch != 0) calo_ht();
	if (calo_mht_eta_branch != 0) calo_mht_eta();
	if (calo_mht_phi_branch != 0) calo_mht_phi();
	if (calo_mht_pt_branch != 0) calo_mht_pt();
	if (gen_ht_branch != 0) gen_ht();
	if (gen_met_branch != 0) gen_met();
	if (gen_metcalo_branch != 0) gen_metcalo();
	if (gen_metcalononprompt_branch != 0) gen_metcalononprompt();
	if (met_eta_branch != 0) met_eta();
	if (met_phi_branch != 0) met_phi();
	if (met_pt_branch != 0) met_pt();
	if (pf_ht_branch != 0) pf_ht();
	if (pf_mht_eta_branch != 0) pf_mht_eta();
	if (pf_mht_phi_branch != 0) pf_mht_phi();
	if (pf_mht_pt_branch != 0) pf_mht_pt();
	if (wl1ht200_branch != 0) wl1ht200();
	if (bjets_csv_branch != 0) bjets_csv();
	if (bjets_eta_branch != 0) bjets_eta();
	if (bjets_phi_branch != 0) bjets_phi();
	if (bjets_pt_branch != 0) bjets_pt();
	if (ele_eta_branch != 0) ele_eta();
	if (ele_phi_branch != 0) ele_phi();
	if (ele_pt_branch != 0) ele_pt();
	if (els_clustershape_branch != 0) els_clustershape();
	if (els_deta_branch != 0) els_deta();
	if (els_dphi_branch != 0) els_dphi();
	if (els_ecal_iso_branch != 0) els_ecal_iso();
	if (els_eminusp_branch != 0) els_eminusp();
	if (els_eta_branch != 0) els_eta();
	if (els_hcal_iso_branch != 0) els_hcal_iso();
	if (els_he_branch != 0) els_he();
	if (els_phi_branch != 0) els_phi();
	if (els_pt_branch != 0) els_pt();
	if (els_track_iso_branch != 0) els_track_iso();
	if (genels_eta_branch != 0) genels_eta();
	if (genels_phi_branch != 0) genels_phi();
	if (genels_pt_branch != 0) genels_pt();
	if (genjets_eta_branch != 0) genjets_eta();
	if (genjets_phi_branch != 0) genjets_phi();
	if (genjets_pt_branch != 0) genjets_pt();
	if (genmus_eta_branch != 0) genmus_eta();
	if (genmus_phi_branch != 0) genmus_phi();
	if (genmus_pt_branch != 0) genmus_pt();
	if (gentop_eta_branch != 0) gentop_eta();
	if (gentop_phi_branch != 0) gentop_phi();
	if (gentop_pt_branch != 0) gentop_pt();
	if (mus_eta_branch != 0) mus_eta();
	if (mus_iso_branch != 0) mus_iso();
	if (mus_phi_branch != 0) mus_phi();
	if (mus_pt_branch != 0) mus_pt();
	if (pfjets_eta_branch != 0) pfjets_eta();
	if (pfjets_phi_branch != 0) pfjets_phi();
	if (pfjets_pt_branch != 0) pfjets_pt();
	if (reco_els_eta_branch != 0) reco_els_eta();
	if (reco_els_phi_branch != 0) reco_els_phi();
	if (reco_els_pt_branch != 0) reco_els_pt();
	if (reco_genjet_eta_branch != 0) reco_genjet_eta();
	if (reco_genjet_phi_branch != 0) reco_genjet_phi();
	if (reco_genjet_pt_branch != 0) reco_genjet_pt();
	if (reco_jet_eta_branch != 0) reco_jet_eta();
	if (reco_jet_phi_branch != 0) reco_jet_phi();
	if (reco_jet_pt_branch != 0) reco_jet_pt();
	if (reco_mus_eta_branch != 0) reco_mus_eta();
	if (reco_mus_phi_branch != 0) reco_mus_phi();
	if (reco_mus_pt_branch != 0) reco_mus_pt();
	if (genels_ggmom_id_branch != 0) genels_ggmom_id();
	if (genels_gmom_id_branch != 0) genels_gmom_id();
	if (genels_mom_id_branch != 0) genels_mom_id();
	if (genmus_ggmom_id_branch != 0) genmus_ggmom_id();
	if (genmus_gmom_id_branch != 0) genmus_gmom_id();
	if (genmus_mom_id_branch != 0) genmus_mom_id();
	if (gentop_id_branch != 0) gentop_id();
	if (gentop_mom_id_branch != 0) gentop_mom_id();
	if (ngenlep_branch != 0) ngenlep();
}

	float &calo_ht()
	{
		if (not calo_ht_isLoaded) {
			if (calo_ht_branch != 0) {
				calo_ht_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_ht_)) {
					printf("branch calo_ht_branch contains a bad float: %f\n", calo_ht_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_ht_branch does not exist!\n");
				exit(1);
			}
			calo_ht_isLoaded = true;
		}
		return calo_ht_;
	}
	float &calo_mht_eta()
	{
		if (not calo_mht_eta_isLoaded) {
			if (calo_mht_eta_branch != 0) {
				calo_mht_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_mht_eta_)) {
					printf("branch calo_mht_eta_branch contains a bad float: %f\n", calo_mht_eta_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_mht_eta_branch does not exist!\n");
				exit(1);
			}
			calo_mht_eta_isLoaded = true;
		}
		return calo_mht_eta_;
	}
	float &calo_mht_phi()
	{
		if (not calo_mht_phi_isLoaded) {
			if (calo_mht_phi_branch != 0) {
				calo_mht_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_mht_phi_)) {
					printf("branch calo_mht_phi_branch contains a bad float: %f\n", calo_mht_phi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_mht_phi_branch does not exist!\n");
				exit(1);
			}
			calo_mht_phi_isLoaded = true;
		}
		return calo_mht_phi_;
	}
	float &calo_mht_pt()
	{
		if (not calo_mht_pt_isLoaded) {
			if (calo_mht_pt_branch != 0) {
				calo_mht_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_mht_pt_)) {
					printf("branch calo_mht_pt_branch contains a bad float: %f\n", calo_mht_pt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_mht_pt_branch does not exist!\n");
				exit(1);
			}
			calo_mht_pt_isLoaded = true;
		}
		return calo_mht_pt_;
	}
	float &gen_ht()
	{
		if (not gen_ht_isLoaded) {
			if (gen_ht_branch != 0) {
				gen_ht_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(gen_ht_)) {
					printf("branch gen_ht_branch contains a bad float: %f\n", gen_ht_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gen_ht_branch does not exist!\n");
				exit(1);
			}
			gen_ht_isLoaded = true;
		}
		return gen_ht_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) {
				gen_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(gen_met_)) {
					printf("branch gen_met_branch contains a bad float: %f\n", gen_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gen_met_branch does not exist!\n");
				exit(1);
			}
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &gen_metcalo()
	{
		if (not gen_metcalo_isLoaded) {
			if (gen_metcalo_branch != 0) {
				gen_metcalo_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(gen_metcalo_)) {
					printf("branch gen_metcalo_branch contains a bad float: %f\n", gen_metcalo_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gen_metcalo_branch does not exist!\n");
				exit(1);
			}
			gen_metcalo_isLoaded = true;
		}
		return gen_metcalo_;
	}
	float &gen_metcalononprompt()
	{
		if (not gen_metcalononprompt_isLoaded) {
			if (gen_metcalononprompt_branch != 0) {
				gen_metcalononprompt_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(gen_metcalononprompt_)) {
					printf("branch gen_metcalononprompt_branch contains a bad float: %f\n", gen_metcalononprompt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gen_metcalononprompt_branch does not exist!\n");
				exit(1);
			}
			gen_metcalononprompt_isLoaded = true;
		}
		return gen_metcalononprompt_;
	}
	float &met_eta()
	{
		if (not met_eta_isLoaded) {
			if (met_eta_branch != 0) {
				met_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(met_eta_)) {
					printf("branch met_eta_branch contains a bad float: %f\n", met_eta_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch met_eta_branch does not exist!\n");
				exit(1);
			}
			met_eta_isLoaded = true;
		}
		return met_eta_;
	}
	float &met_phi()
	{
		if (not met_phi_isLoaded) {
			if (met_phi_branch != 0) {
				met_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(met_phi_)) {
					printf("branch met_phi_branch contains a bad float: %f\n", met_phi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch met_phi_branch does not exist!\n");
				exit(1);
			}
			met_phi_isLoaded = true;
		}
		return met_phi_;
	}
	float &met_pt()
	{
		if (not met_pt_isLoaded) {
			if (met_pt_branch != 0) {
				met_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(met_pt_)) {
					printf("branch met_pt_branch contains a bad float: %f\n", met_pt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch met_pt_branch does not exist!\n");
				exit(1);
			}
			met_pt_isLoaded = true;
		}
		return met_pt_;
	}
	float &pf_ht()
	{
		if (not pf_ht_isLoaded) {
			if (pf_ht_branch != 0) {
				pf_ht_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_ht_)) {
					printf("branch pf_ht_branch contains a bad float: %f\n", pf_ht_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ht_branch does not exist!\n");
				exit(1);
			}
			pf_ht_isLoaded = true;
		}
		return pf_ht_;
	}
	float &pf_mht_eta()
	{
		if (not pf_mht_eta_isLoaded) {
			if (pf_mht_eta_branch != 0) {
				pf_mht_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_mht_eta_)) {
					printf("branch pf_mht_eta_branch contains a bad float: %f\n", pf_mht_eta_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_mht_eta_branch does not exist!\n");
				exit(1);
			}
			pf_mht_eta_isLoaded = true;
		}
		return pf_mht_eta_;
	}
	float &pf_mht_phi()
	{
		if (not pf_mht_phi_isLoaded) {
			if (pf_mht_phi_branch != 0) {
				pf_mht_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_mht_phi_)) {
					printf("branch pf_mht_phi_branch contains a bad float: %f\n", pf_mht_phi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_mht_phi_branch does not exist!\n");
				exit(1);
			}
			pf_mht_phi_isLoaded = true;
		}
		return pf_mht_phi_;
	}
	float &pf_mht_pt()
	{
		if (not pf_mht_pt_isLoaded) {
			if (pf_mht_pt_branch != 0) {
				pf_mht_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_mht_pt_)) {
					printf("branch pf_mht_pt_branch contains a bad float: %f\n", pf_mht_pt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_mht_pt_branch does not exist!\n");
				exit(1);
			}
			pf_mht_pt_isLoaded = true;
		}
		return pf_mht_pt_;
	}
	float &wl1ht200()
	{
		if (not wl1ht200_isLoaded) {
			if (wl1ht200_branch != 0) {
				wl1ht200_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(wl1ht200_)) {
					printf("branch wl1ht200_branch contains a bad float: %f\n", wl1ht200_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch wl1ht200_branch does not exist!\n");
				exit(1);
			}
			wl1ht200_isLoaded = true;
		}
		return wl1ht200_;
	}
	const vector<float> &bjets_csv()
	{
		if (not bjets_csv_isLoaded) {
			if (bjets_csv_branch != 0) {
				bjets_csv_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = bjets_csv_.begin(); i != bjets_csv_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch bjets_csv_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch bjets_csv_branch does not exist!\n");
				exit(1);
			}
			bjets_csv_isLoaded = true;
		}
		return bjets_csv_;
	}
	const vector<float> &bjets_eta()
	{
		if (not bjets_eta_isLoaded) {
			if (bjets_eta_branch != 0) {
				bjets_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = bjets_eta_.begin(); i != bjets_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch bjets_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch bjets_eta_branch does not exist!\n");
				exit(1);
			}
			bjets_eta_isLoaded = true;
		}
		return bjets_eta_;
	}
	const vector<float> &bjets_phi()
	{
		if (not bjets_phi_isLoaded) {
			if (bjets_phi_branch != 0) {
				bjets_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = bjets_phi_.begin(); i != bjets_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch bjets_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch bjets_phi_branch does not exist!\n");
				exit(1);
			}
			bjets_phi_isLoaded = true;
		}
		return bjets_phi_;
	}
	const vector<float> &bjets_pt()
	{
		if (not bjets_pt_isLoaded) {
			if (bjets_pt_branch != 0) {
				bjets_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = bjets_pt_.begin(); i != bjets_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch bjets_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch bjets_pt_branch does not exist!\n");
				exit(1);
			}
			bjets_pt_isLoaded = true;
		}
		return bjets_pt_;
	}
	const vector<float> &ele_eta()
	{
		if (not ele_eta_isLoaded) {
			if (ele_eta_branch != 0) {
				ele_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = ele_eta_.begin(); i != ele_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch ele_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ele_eta_branch does not exist!\n");
				exit(1);
			}
			ele_eta_isLoaded = true;
		}
		return ele_eta_;
	}
	const vector<float> &ele_phi()
	{
		if (not ele_phi_isLoaded) {
			if (ele_phi_branch != 0) {
				ele_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = ele_phi_.begin(); i != ele_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch ele_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ele_phi_branch does not exist!\n");
				exit(1);
			}
			ele_phi_isLoaded = true;
		}
		return ele_phi_;
	}
	const vector<float> &ele_pt()
	{
		if (not ele_pt_isLoaded) {
			if (ele_pt_branch != 0) {
				ele_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = ele_pt_.begin(); i != ele_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch ele_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ele_pt_branch does not exist!\n");
				exit(1);
			}
			ele_pt_isLoaded = true;
		}
		return ele_pt_;
	}
	const vector<float> &els_clustershape()
	{
		if (not els_clustershape_isLoaded) {
			if (els_clustershape_branch != 0) {
				els_clustershape_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_clustershape_.begin(); i != els_clustershape_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_clustershape_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_clustershape_branch does not exist!\n");
				exit(1);
			}
			els_clustershape_isLoaded = true;
		}
		return els_clustershape_;
	}
	const vector<float> &els_deta()
	{
		if (not els_deta_isLoaded) {
			if (els_deta_branch != 0) {
				els_deta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_deta_.begin(); i != els_deta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_deta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_deta_branch does not exist!\n");
				exit(1);
			}
			els_deta_isLoaded = true;
		}
		return els_deta_;
	}
	const vector<float> &els_dphi()
	{
		if (not els_dphi_isLoaded) {
			if (els_dphi_branch != 0) {
				els_dphi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dphi_.begin(); i != els_dphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_dphi_branch does not exist!\n");
				exit(1);
			}
			els_dphi_isLoaded = true;
		}
		return els_dphi_;
	}
	const vector<float> &els_ecal_iso()
	{
		if (not els_ecal_iso_isLoaded) {
			if (els_ecal_iso_branch != 0) {
				els_ecal_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ecal_iso_.begin(); i != els_ecal_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ecal_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_ecal_iso_branch does not exist!\n");
				exit(1);
			}
			els_ecal_iso_isLoaded = true;
		}
		return els_ecal_iso_;
	}
	const vector<float> &els_eminusp()
	{
		if (not els_eminusp_isLoaded) {
			if (els_eminusp_branch != 0) {
				els_eminusp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eminusp_.begin(); i != els_eminusp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eminusp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_eminusp_branch does not exist!\n");
				exit(1);
			}
			els_eminusp_isLoaded = true;
		}
		return els_eminusp_;
	}
	const vector<float> &els_eta()
	{
		if (not els_eta_isLoaded) {
			if (els_eta_branch != 0) {
				els_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eta_.begin(); i != els_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_eta_branch does not exist!\n");
				exit(1);
			}
			els_eta_isLoaded = true;
		}
		return els_eta_;
	}
	const vector<float> &els_hcal_iso()
	{
		if (not els_hcal_iso_isLoaded) {
			if (els_hcal_iso_branch != 0) {
				els_hcal_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_hcal_iso_.begin(); i != els_hcal_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_hcal_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_hcal_iso_branch does not exist!\n");
				exit(1);
			}
			els_hcal_iso_isLoaded = true;
		}
		return els_hcal_iso_;
	}
	const vector<float> &els_he()
	{
		if (not els_he_isLoaded) {
			if (els_he_branch != 0) {
				els_he_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_he_.begin(); i != els_he_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_he_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_he_branch does not exist!\n");
				exit(1);
			}
			els_he_isLoaded = true;
		}
		return els_he_;
	}
	const vector<float> &els_phi()
	{
		if (not els_phi_isLoaded) {
			if (els_phi_branch != 0) {
				els_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_phi_.begin(); i != els_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_phi_branch does not exist!\n");
				exit(1);
			}
			els_phi_isLoaded = true;
		}
		return els_phi_;
	}
	const vector<float> &els_pt()
	{
		if (not els_pt_isLoaded) {
			if (els_pt_branch != 0) {
				els_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pt_.begin(); i != els_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_pt_branch does not exist!\n");
				exit(1);
			}
			els_pt_isLoaded = true;
		}
		return els_pt_;
	}
	const vector<float> &els_track_iso()
	{
		if (not els_track_iso_isLoaded) {
			if (els_track_iso_branch != 0) {
				els_track_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_track_iso_.begin(); i != els_track_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_track_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_track_iso_branch does not exist!\n");
				exit(1);
			}
			els_track_iso_isLoaded = true;
		}
		return els_track_iso_;
	}
	const vector<float> &genels_eta()
	{
		if (not genels_eta_isLoaded) {
			if (genels_eta_branch != 0) {
				genels_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genels_eta_.begin(); i != genels_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genels_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_eta_branch does not exist!\n");
				exit(1);
			}
			genels_eta_isLoaded = true;
		}
		return genels_eta_;
	}
	const vector<float> &genels_phi()
	{
		if (not genels_phi_isLoaded) {
			if (genels_phi_branch != 0) {
				genels_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genels_phi_.begin(); i != genels_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genels_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_phi_branch does not exist!\n");
				exit(1);
			}
			genels_phi_isLoaded = true;
		}
		return genels_phi_;
	}
	const vector<float> &genels_pt()
	{
		if (not genels_pt_isLoaded) {
			if (genels_pt_branch != 0) {
				genels_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genels_pt_.begin(); i != genels_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genels_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_pt_branch does not exist!\n");
				exit(1);
			}
			genels_pt_isLoaded = true;
		}
		return genels_pt_;
	}
	const vector<float> &genjets_eta()
	{
		if (not genjets_eta_isLoaded) {
			if (genjets_eta_branch != 0) {
				genjets_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genjets_eta_.begin(); i != genjets_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genjets_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genjets_eta_branch does not exist!\n");
				exit(1);
			}
			genjets_eta_isLoaded = true;
		}
		return genjets_eta_;
	}
	const vector<float> &genjets_phi()
	{
		if (not genjets_phi_isLoaded) {
			if (genjets_phi_branch != 0) {
				genjets_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genjets_phi_.begin(); i != genjets_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genjets_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genjets_phi_branch does not exist!\n");
				exit(1);
			}
			genjets_phi_isLoaded = true;
		}
		return genjets_phi_;
	}
	const vector<float> &genjets_pt()
	{
		if (not genjets_pt_isLoaded) {
			if (genjets_pt_branch != 0) {
				genjets_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genjets_pt_.begin(); i != genjets_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genjets_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genjets_pt_branch does not exist!\n");
				exit(1);
			}
			genjets_pt_isLoaded = true;
		}
		return genjets_pt_;
	}
	const vector<float> &genmus_eta()
	{
		if (not genmus_eta_isLoaded) {
			if (genmus_eta_branch != 0) {
				genmus_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genmus_eta_.begin(); i != genmus_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genmus_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_eta_branch does not exist!\n");
				exit(1);
			}
			genmus_eta_isLoaded = true;
		}
		return genmus_eta_;
	}
	const vector<float> &genmus_phi()
	{
		if (not genmus_phi_isLoaded) {
			if (genmus_phi_branch != 0) {
				genmus_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genmus_phi_.begin(); i != genmus_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genmus_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_phi_branch does not exist!\n");
				exit(1);
			}
			genmus_phi_isLoaded = true;
		}
		return genmus_phi_;
	}
	const vector<float> &genmus_pt()
	{
		if (not genmus_pt_isLoaded) {
			if (genmus_pt_branch != 0) {
				genmus_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = genmus_pt_.begin(); i != genmus_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch genmus_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_pt_branch does not exist!\n");
				exit(1);
			}
			genmus_pt_isLoaded = true;
		}
		return genmus_pt_;
	}
	const vector<float> &gentop_eta()
	{
		if (not gentop_eta_isLoaded) {
			if (gentop_eta_branch != 0) {
				gentop_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = gentop_eta_.begin(); i != gentop_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch gentop_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gentop_eta_branch does not exist!\n");
				exit(1);
			}
			gentop_eta_isLoaded = true;
		}
		return gentop_eta_;
	}
	const vector<float> &gentop_phi()
	{
		if (not gentop_phi_isLoaded) {
			if (gentop_phi_branch != 0) {
				gentop_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = gentop_phi_.begin(); i != gentop_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch gentop_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gentop_phi_branch does not exist!\n");
				exit(1);
			}
			gentop_phi_isLoaded = true;
		}
		return gentop_phi_;
	}
	const vector<float> &gentop_pt()
	{
		if (not gentop_pt_isLoaded) {
			if (gentop_pt_branch != 0) {
				gentop_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = gentop_pt_.begin(); i != gentop_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch gentop_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gentop_pt_branch does not exist!\n");
				exit(1);
			}
			gentop_pt_isLoaded = true;
		}
		return gentop_pt_;
	}
	const vector<float> &mus_eta()
	{
		if (not mus_eta_isLoaded) {
			if (mus_eta_branch != 0) {
				mus_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_eta_.begin(); i != mus_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_eta_branch does not exist!\n");
				exit(1);
			}
			mus_eta_isLoaded = true;
		}
		return mus_eta_;
	}
	const vector<float> &mus_iso()
	{
		if (not mus_iso_isLoaded) {
			if (mus_iso_branch != 0) {
				mus_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso_.begin(); i != mus_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_iso_branch does not exist!\n");
				exit(1);
			}
			mus_iso_isLoaded = true;
		}
		return mus_iso_;
	}
	const vector<float> &mus_phi()
	{
		if (not mus_phi_isLoaded) {
			if (mus_phi_branch != 0) {
				mus_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_phi_.begin(); i != mus_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_phi_branch does not exist!\n");
				exit(1);
			}
			mus_phi_isLoaded = true;
		}
		return mus_phi_;
	}
	const vector<float> &mus_pt()
	{
		if (not mus_pt_isLoaded) {
			if (mus_pt_branch != 0) {
				mus_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pt_.begin(); i != mus_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pt_branch does not exist!\n");
				exit(1);
			}
			mus_pt_isLoaded = true;
		}
		return mus_pt_;
	}
	const vector<float> &pfjets_eta()
	{
		if (not pfjets_eta_isLoaded) {
			if (pfjets_eta_branch != 0) {
				pfjets_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pfjets_eta_.begin(); i != pfjets_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pfjets_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfjets_eta_branch does not exist!\n");
				exit(1);
			}
			pfjets_eta_isLoaded = true;
		}
		return pfjets_eta_;
	}
	const vector<float> &pfjets_phi()
	{
		if (not pfjets_phi_isLoaded) {
			if (pfjets_phi_branch != 0) {
				pfjets_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pfjets_phi_.begin(); i != pfjets_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pfjets_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfjets_phi_branch does not exist!\n");
				exit(1);
			}
			pfjets_phi_isLoaded = true;
		}
		return pfjets_phi_;
	}
	const vector<float> &pfjets_pt()
	{
		if (not pfjets_pt_isLoaded) {
			if (pfjets_pt_branch != 0) {
				pfjets_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pfjets_pt_.begin(); i != pfjets_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pfjets_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfjets_pt_branch does not exist!\n");
				exit(1);
			}
			pfjets_pt_isLoaded = true;
		}
		return pfjets_pt_;
	}
	const vector<float> &reco_els_eta()
	{
		if (not reco_els_eta_isLoaded) {
			if (reco_els_eta_branch != 0) {
				reco_els_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_els_eta_.begin(); i != reco_els_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_els_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_els_eta_branch does not exist!\n");
				exit(1);
			}
			reco_els_eta_isLoaded = true;
		}
		return reco_els_eta_;
	}
	const vector<float> &reco_els_phi()
	{
		if (not reco_els_phi_isLoaded) {
			if (reco_els_phi_branch != 0) {
				reco_els_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_els_phi_.begin(); i != reco_els_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_els_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_els_phi_branch does not exist!\n");
				exit(1);
			}
			reco_els_phi_isLoaded = true;
		}
		return reco_els_phi_;
	}
	const vector<float> &reco_els_pt()
	{
		if (not reco_els_pt_isLoaded) {
			if (reco_els_pt_branch != 0) {
				reco_els_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_els_pt_.begin(); i != reco_els_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_els_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_els_pt_branch does not exist!\n");
				exit(1);
			}
			reco_els_pt_isLoaded = true;
		}
		return reco_els_pt_;
	}
	const vector<float> &reco_genjet_eta()
	{
		if (not reco_genjet_eta_isLoaded) {
			if (reco_genjet_eta_branch != 0) {
				reco_genjet_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_genjet_eta_.begin(); i != reco_genjet_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_genjet_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_genjet_eta_branch does not exist!\n");
				exit(1);
			}
			reco_genjet_eta_isLoaded = true;
		}
		return reco_genjet_eta_;
	}
	const vector<float> &reco_genjet_phi()
	{
		if (not reco_genjet_phi_isLoaded) {
			if (reco_genjet_phi_branch != 0) {
				reco_genjet_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_genjet_phi_.begin(); i != reco_genjet_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_genjet_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_genjet_phi_branch does not exist!\n");
				exit(1);
			}
			reco_genjet_phi_isLoaded = true;
		}
		return reco_genjet_phi_;
	}
	const vector<float> &reco_genjet_pt()
	{
		if (not reco_genjet_pt_isLoaded) {
			if (reco_genjet_pt_branch != 0) {
				reco_genjet_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_genjet_pt_.begin(); i != reco_genjet_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_genjet_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_genjet_pt_branch does not exist!\n");
				exit(1);
			}
			reco_genjet_pt_isLoaded = true;
		}
		return reco_genjet_pt_;
	}
	const vector<float> &reco_jet_eta()
	{
		if (not reco_jet_eta_isLoaded) {
			if (reco_jet_eta_branch != 0) {
				reco_jet_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_jet_eta_.begin(); i != reco_jet_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_jet_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_jet_eta_branch does not exist!\n");
				exit(1);
			}
			reco_jet_eta_isLoaded = true;
		}
		return reco_jet_eta_;
	}
	const vector<float> &reco_jet_phi()
	{
		if (not reco_jet_phi_isLoaded) {
			if (reco_jet_phi_branch != 0) {
				reco_jet_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_jet_phi_.begin(); i != reco_jet_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_jet_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_jet_phi_branch does not exist!\n");
				exit(1);
			}
			reco_jet_phi_isLoaded = true;
		}
		return reco_jet_phi_;
	}
	const vector<float> &reco_jet_pt()
	{
		if (not reco_jet_pt_isLoaded) {
			if (reco_jet_pt_branch != 0) {
				reco_jet_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_jet_pt_.begin(); i != reco_jet_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_jet_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_jet_pt_branch does not exist!\n");
				exit(1);
			}
			reco_jet_pt_isLoaded = true;
		}
		return reco_jet_pt_;
	}
	const vector<float> &reco_mus_eta()
	{
		if (not reco_mus_eta_isLoaded) {
			if (reco_mus_eta_branch != 0) {
				reco_mus_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_mus_eta_.begin(); i != reco_mus_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_mus_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_mus_eta_branch does not exist!\n");
				exit(1);
			}
			reco_mus_eta_isLoaded = true;
		}
		return reco_mus_eta_;
	}
	const vector<float> &reco_mus_phi()
	{
		if (not reco_mus_phi_isLoaded) {
			if (reco_mus_phi_branch != 0) {
				reco_mus_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_mus_phi_.begin(); i != reco_mus_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_mus_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_mus_phi_branch does not exist!\n");
				exit(1);
			}
			reco_mus_phi_isLoaded = true;
		}
		return reco_mus_phi_;
	}
	const vector<float> &reco_mus_pt()
	{
		if (not reco_mus_pt_isLoaded) {
			if (reco_mus_pt_branch != 0) {
				reco_mus_pt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = reco_mus_pt_.begin(); i != reco_mus_pt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch reco_mus_pt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch reco_mus_pt_branch does not exist!\n");
				exit(1);
			}
			reco_mus_pt_isLoaded = true;
		}
		return reco_mus_pt_;
	}
	const vector<int> &genels_ggmom_id()
	{
		if (not genels_ggmom_id_isLoaded) {
			if (genels_ggmom_id_branch != 0) {
				genels_ggmom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_ggmom_id_branch does not exist!\n");
				exit(1);
			}
			genels_ggmom_id_isLoaded = true;
		}
		return genels_ggmom_id_;
	}
	const vector<int> &genels_gmom_id()
	{
		if (not genels_gmom_id_isLoaded) {
			if (genels_gmom_id_branch != 0) {
				genels_gmom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_gmom_id_branch does not exist!\n");
				exit(1);
			}
			genels_gmom_id_isLoaded = true;
		}
		return genels_gmom_id_;
	}
	const vector<int> &genels_mom_id()
	{
		if (not genels_mom_id_isLoaded) {
			if (genels_mom_id_branch != 0) {
				genels_mom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genels_mom_id_branch does not exist!\n");
				exit(1);
			}
			genels_mom_id_isLoaded = true;
		}
		return genels_mom_id_;
	}
	const vector<int> &genmus_ggmom_id()
	{
		if (not genmus_ggmom_id_isLoaded) {
			if (genmus_ggmom_id_branch != 0) {
				genmus_ggmom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_ggmom_id_branch does not exist!\n");
				exit(1);
			}
			genmus_ggmom_id_isLoaded = true;
		}
		return genmus_ggmom_id_;
	}
	const vector<int> &genmus_gmom_id()
	{
		if (not genmus_gmom_id_isLoaded) {
			if (genmus_gmom_id_branch != 0) {
				genmus_gmom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_gmom_id_branch does not exist!\n");
				exit(1);
			}
			genmus_gmom_id_isLoaded = true;
		}
		return genmus_gmom_id_;
	}
	const vector<int> &genmus_mom_id()
	{
		if (not genmus_mom_id_isLoaded) {
			if (genmus_mom_id_branch != 0) {
				genmus_mom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmus_mom_id_branch does not exist!\n");
				exit(1);
			}
			genmus_mom_id_isLoaded = true;
		}
		return genmus_mom_id_;
	}
	const vector<int> &gentop_id()
	{
		if (not gentop_id_isLoaded) {
			if (gentop_id_branch != 0) {
				gentop_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gentop_id_branch does not exist!\n");
				exit(1);
			}
			gentop_id_isLoaded = true;
		}
		return gentop_id_;
	}
	const vector<int> &gentop_mom_id()
	{
		if (not gentop_mom_id_isLoaded) {
			if (gentop_mom_id_branch != 0) {
				gentop_mom_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gentop_mom_id_branch does not exist!\n");
				exit(1);
			}
			gentop_mom_id_isLoaded = true;
		}
		return gentop_mom_id_;
	}
	unsigned int &ngenlep()
	{
		if (not ngenlep_isLoaded) {
			if (ngenlep_branch != 0) {
				ngenlep_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ngenlep_branch does not exist!\n");
				exit(1);
			}
			ngenlep_isLoaded = true;
		}
		return ngenlep_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern hlt_class hlt;
#endif

namespace baby {
	const float &calo_ht();
	const float &calo_mht_eta();
	const float &calo_mht_phi();
	const float &calo_mht_pt();
	const float &gen_ht();
	const float &gen_met();
	const float &gen_metcalo();
	const float &gen_metcalononprompt();
	const float &met_eta();
	const float &met_phi();
	const float &met_pt();
	const float &pf_ht();
	const float &pf_mht_eta();
	const float &pf_mht_phi();
	const float &pf_mht_pt();
	const float &wl1ht200();
	const vector<float> &bjets_csv();
	const vector<float> &bjets_eta();
	const vector<float> &bjets_phi();
	const vector<float> &bjets_pt();
	const vector<float> &ele_eta();
	const vector<float> &ele_phi();
	const vector<float> &ele_pt();
	const vector<float> &els_clustershape();
	const vector<float> &els_deta();
	const vector<float> &els_dphi();
	const vector<float> &els_ecal_iso();
	const vector<float> &els_eminusp();
	const vector<float> &els_eta();
	const vector<float> &els_hcal_iso();
	const vector<float> &els_he();
	const vector<float> &els_phi();
	const vector<float> &els_pt();
	const vector<float> &els_track_iso();
	const vector<float> &genels_eta();
	const vector<float> &genels_phi();
	const vector<float> &genels_pt();
	const vector<float> &genjets_eta();
	const vector<float> &genjets_phi();
	const vector<float> &genjets_pt();
	const vector<float> &genmus_eta();
	const vector<float> &genmus_phi();
	const vector<float> &genmus_pt();
	const vector<float> &gentop_eta();
	const vector<float> &gentop_phi();
	const vector<float> &gentop_pt();
	const vector<float> &mus_eta();
	const vector<float> &mus_iso();
	const vector<float> &mus_phi();
	const vector<float> &mus_pt();
	const vector<float> &pfjets_eta();
	const vector<float> &pfjets_phi();
	const vector<float> &pfjets_pt();
	const vector<float> &reco_els_eta();
	const vector<float> &reco_els_phi();
	const vector<float> &reco_els_pt();
	const vector<float> &reco_genjet_eta();
	const vector<float> &reco_genjet_phi();
	const vector<float> &reco_genjet_pt();
	const vector<float> &reco_jet_eta();
	const vector<float> &reco_jet_phi();
	const vector<float> &reco_jet_pt();
	const vector<float> &reco_mus_eta();
	const vector<float> &reco_mus_phi();
	const vector<float> &reco_mus_pt();
	const vector<int> &genels_ggmom_id();
	const vector<int> &genels_gmom_id();
	const vector<int> &genels_mom_id();
	const vector<int> &genmus_ggmom_id();
	const vector<int> &genmus_gmom_id();
	const vector<int> &genmus_mom_id();
	const vector<int> &gentop_id();
	const vector<int> &gentop_mom_id();
	const unsigned int &ngenlep();
}
#endif
