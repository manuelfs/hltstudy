#include "hlt_class.h"
hlt_class hlt;
namespace baby {
	const float &gen_ht() { return hlt.gen_ht(); }
	const float &gen_met() { return hlt.gen_met(); }
	const float &gen_metcalo() { return hlt.gen_metcalo(); }
	const float &gen_metcalononprompt() { return hlt.gen_metcalononprompt(); }
	const float &met_eta() { return hlt.met_eta(); }
	const float &met_phi() { return hlt.met_phi(); }
	const float &met_pt() { return hlt.met_pt(); }
	const float &pf_ht() { return hlt.pf_ht(); }
	const float &pf_mht_eta() { return hlt.pf_mht_eta(); }
	const float &pf_mht_phi() { return hlt.pf_mht_phi(); }
	const float &pf_mht_pt() { return hlt.pf_mht_pt(); }
	const float &wl1ht200() { return hlt.wl1ht200(); }
	const vector<float> &els_eta() { return hlt.els_eta(); }
	const vector<float> &els_phi() { return hlt.els_phi(); }
	const vector<float> &els_pt() { return hlt.els_pt(); }
	const vector<float> &genels_eta() { return hlt.genels_eta(); }
	const vector<float> &genels_phi() { return hlt.genels_phi(); }
	const vector<float> &genels_pt() { return hlt.genels_pt(); }
	const vector<float> &genjets_eta() { return hlt.genjets_eta(); }
	const vector<float> &genjets_phi() { return hlt.genjets_phi(); }
	const vector<float> &genjets_pt() { return hlt.genjets_pt(); }
	const vector<float> &genmus_eta() { return hlt.genmus_eta(); }
	const vector<float> &genmus_phi() { return hlt.genmus_phi(); }
	const vector<float> &genmus_pt() { return hlt.genmus_pt(); }
	const vector<float> &mus_eta() { return hlt.mus_eta(); }
	const vector<float> &mus_phi() { return hlt.mus_phi(); }
	const vector<float> &mus_pt() { return hlt.mus_pt(); }
	const vector<float> &pfjets_eta() { return hlt.pfjets_eta(); }
	const vector<float> &pfjets_phi() { return hlt.pfjets_phi(); }
	const vector<float> &pfjets_pt() { return hlt.pfjets_pt(); }
	const unsigned int &ngenlep() { return hlt.ngenlep(); }
}
