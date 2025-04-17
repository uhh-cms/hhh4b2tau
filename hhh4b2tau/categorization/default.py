# coding: utf-8

"""
Exemplary selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import
from columnflow.columnar_util import attach_coffea_behavior,  EMPTY_FLOAT

from hhh4b2tau.config.variables import build_higgs_reco

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
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2

@categorizer(uses={"Jet.{pt,phi,eta,mass}"})
def cat_4j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # four or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 4

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_1btag(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # one or more b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 1

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_2btag(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # one or more b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 2

@categorizer(uses={"Jet.{pt,phi,eta,mass,btagDeepFlavB}"})
def cat_3btag(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # one or more b jets
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    btag_mask = (events.Jet.btagDeepFlavB >= btag_wp)
    return events, ak.sum(btag_mask, axis=1) >= 3

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

# categorizer class to access created variables more conveniently
class _VarCut(Categorizer):
    def init_func(self):
        self.uses = {"{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}", "BB{1,2}_idx"}

    def rebuild_higgs_reco_mask(self, events: ak.Array, obj: str, var: str, lower = None, higher = None):
        value = build_higgs_reco(events, obj, var)
        self.emptyf_mask = (value != EMPTY_FLOAT)
        self.lower_mask = np.ones(len(events), dtype=bool)
        self.higher_mask = np.ones(len(events), dtype=bool)
        if lower != None:
            self.lower_mask = (value < lower)
        if higher != None:
            self.higher_mask = (value >= higher)
        self.final_mask = self.emptyf_mask & self.lower_mask & self.higher_mask


@_VarCut.categorizer()          
def cat_h3_mass(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", lower = 125)
    return events, self.final_mask

@_VarCut.categorizer()          
def cat_h3_mass_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "h3", "mass", higher = 125)
    return events, self.final_mask


@_VarCut.categorizer()
def cat_leps_cos(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", higher = -0.25)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_leps_cos_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = -0.25)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_leps_dr(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", lower = 2.4)
    return events, self.final_mask

@_VarCut.categorizer()
def cat_leps_dr_orth(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    self.rebuild_higgs_reco_mask(events, "leps", "dr", higher = 2.4)
    return events, self.final_mask
