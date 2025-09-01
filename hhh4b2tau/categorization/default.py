# coding: utf-8

"""
Exemplary selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import
from columnflow.columnar_util import attach_coffea_behavior,  EMPTY_FLOAT, optional_column as optional

from hhh4b2tau.config.variables import build_higgs_reco, build_gen_higgs, build_higgs_reco_gen

ak = maybe_import("awkward")
np = maybe_import("numpy")


#
# lepton channels
#

@categorizer(uses={"channel_id"})
def cat_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.etau.id


@categorizer(uses={"channel_id"})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.mutau.id


@categorizer(uses={"channel_id"})
def cat_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.tautau.id


@categorizer(uses={"channel_id"})
def cat_ee(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.ee.id


@categorizer(uses={"channel_id"})
def cat_mumu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.mumu.id


@categorizer(uses={"channel_id"})
def cat_emu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.channel_id == self.config_inst.channels.n.emu.id


#
# QCD regions
#

@categorizer(uses={"leptons_os"})
def cat_os(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # oppositive sign leptons
    return events, events.leptons_os == 1


@categorizer(uses={"leptons_os"})
def cat_ss(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # same sign leptons
    return events, events.leptons_os == 0


@categorizer(uses={"tau2_isolated"})
def cat_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # isolated tau2
    return events, events.tau2_isolated == 1


@categorizer(uses={"tau2_isolated"})
def cat_noniso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # non-isolated tau2
    return events, events.tau2_isolated == 0


#
# kinematic regions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1

@categorizer(uses={"Jet.{pt,phi,eta,mass}"})
def cat_3j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # exactly three jets
    return events, ak.num(events.Jet.pt, axis=1) == 3

@categorizer(uses={"Jet.{pt,phi,eta,mass}"})
def cat_4j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # four or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 4

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_3b0j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # exactly 3 b jets and no additional non b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    threeb_mask = (ak.sum(btag_mask, axis=1) == 3)
    zeroj_mask = (ak.sum(~btag_mask, axis=1) == 0)
    return events, (threeb_mask & zeroj_mask)

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_3b1j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # exactly 3 b jets and one or more additional non b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    threeb_mask = (ak.sum(btag_mask, axis=1) == 3)
    onej_mask = (ak.sum(~btag_mask, axis=1) >= 1)
    return events, (threeb_mask & onej_mask)

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_1b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 1

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_2b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 2

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_3b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 3

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_4b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # four or more b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 4

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_5b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # four or more b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 5



@categorizer(uses={"Jet.{pt,phi,eta,mass}"})
def cat_j1_pt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # cut on leading jet pt
    return events, ak.fill_none((ak.firsts(events.Jet.pt[:, :1]) >= 100), False)

@categorizer(uses={"Jet.{pt,phi,eta,mass}"})
def cat_j1_pt_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # cut on leading jet pt
    return events, ak.fill_none((ak.firsts(events.Jet.pt[:, :1]) < 100), False)

# create standard categorizer to reduce combinatorics
@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}", "leptons_os", "channel_id", "tau2_isolated",})
def cat_standard(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (ak.sum((events.Jet.btagDeepFlavB >= btag_wp), axis=1) >= 2)
    iso_mask = (events.tau2_isolated == 1)
    os_mask = (events.leptons_os == 1)
    channel_mask = (events.channel_id == self.config_inst.channels.n.mutau.id)
    standard_mask = btag_mask & os_mask & iso_mask & channel_mask
    return events, standard_mask

@categorizer(uses={"chi2"})
def cat_chi2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ak.fill_none((events.chi2 <= 0.31), False)

@categorizer(uses={"chi2"})
def cat_chi2_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ak.fill_none((events.chi2 > 0.31), False)

# categorizer class to access created variables more conveniently
class _VarCut(Categorizer):
    def init_func(self):
        self.uses = {
            "{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}", 
            optional("BB1_idx"), optional("BB2_idx"), 
            "Jet.btagDeepFlavB",
        }

    def rebuild_higgs_reco_mask(self, events: ak.Array, obj: str, var: str, lower = None, higher = None):
        value = build_higgs_reco(events, obj, var)
        self.emptyf_mask = (value != EMPTY_FLOAT)
        self.lower_mask = np.ones(len(events), dtype=bool)
        self.higher_mask = np.ones(len(events), dtype=bool)
        if lower != None:
            self.lower_mask = (value <= lower)
        if higher != None:
            self.higher_mask = (value >= higher)
        self.final_mask = self.emptyf_mask & self.lower_mask & self.higher_mask


@_VarCut.categorizer()          
def cat_mh3(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", lower = 125)
    return events, self.final_mask

@_VarCut.categorizer()          
def cat_mh3_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", lower = 125)
    return events, ~self.final_mask


@_VarCut.categorizer()
def cat_leps_cos(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "cos", higher = -0.25)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_leps_cos_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "cos", higher = -0.25)
    return events, ~self.final_mask

@_VarCut.categorizer()
def cat_leps_dr(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = 2.0)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_leps_dr_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = 2.0)
    return events, ~self.final_mask


@_VarCut.categorizer()
def cat_bb1_dr(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "bb1", "dr", lower = 2.4)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_bb1_dr_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "bb1", "dr", lower = 2.4)
    return events, ~self.final_mask

@_VarCut.categorizer()
def cat_m3j2l(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "3j2l", "mass", higher = 450.0)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_m3j2l_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "3j2l", "mass", higher = 450.0)
    return events, ~self.final_mask

# combine most optimal combination
@_VarCut.categorizer()
def cat_var_cuts(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", lower = 125)
    h3_mass_mask = self.final_mask
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = 2.0)
    leps_dr_mask = self.final_mask
    self.rebuild_higgs_reco_mask(events, "bb1", "dr", lower = 2.4)
    jj1_dr_mask =  self.final_mask
    j1_pt_mask = ak.fill_none((ak.firsts(events.Jet.pt[:, :1]) >= 100), False)
    combo_mask = (h3_mass_mask & leps_dr_mask & jj1_dr_mask & j1_pt_mask)
    return events, combo_mask

@_VarCut.categorizer()
def cat_var_cuts_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", lower = 125)
    h3_mass_mask = self.final_mask
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = 2.0)
    leps_dr_mask = self.final_mask
    self.rebuild_higgs_reco_mask(events, "bb1", "dr", lower = 2.4)
    jj1_dr_mask =  self.final_mask
    j1_pt_mask = ak.fill_none((ak.firsts(events.Jet.pt[:, :1]) >= 100), False)
    combo_mask = (h3_mass_mask & leps_dr_mask & jj1_dr_mask & j1_pt_mask)
    return events, ~combo_mask
