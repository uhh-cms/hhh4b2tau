# coding: utf-8

"""
Exemplary selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import


ak = maybe_import("awkward")


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
    # noon-isolated tau2
    return events, events.tau2_isolated == 0


#
# kinematic regions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1


@categorizer(uses={"Jet.pt"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2


    # using variables to make cuts but not saving them to events yet to avoid issues down stream
    # events = self[detector_variables](events, lepton_results, **kwargs)
    # events = self[detector_variables](events, **kwargs)
    # cos_bb1_mask = ak.fill_none(events.cos_bb1 > -0.25, False)
    # cos_bb2_mask = ak.fill_none(events.cos_bb2 > 0.5, False)
    # cos_tautau_mask = ak.fill_none(events.cos_tautau > 0.0, False)
    # cos_h12_mask = ak.fill_none(events.cos_h12 > -0.6, False)
    # cos_h13_mask = ak.fill_none(events.cos_h13 < 0.75, False)
    # cos_h23_mask = ak.fill_none(events.cos_h23 < 0.6, False)
    # delta_r_bb1_mask = ak.fill_none(events.delta_r_bb1 < 1.8, False)
    # delta_r_bb2_mask = ak.fill_none(events.delta_r_bb2 < 2.0, False)
    # delta_r_tautau_mask = ak.fill_none(events.delta_r_tautau < 2.6, False)
    # delta_r_h12_mask = ak.fill_none((events.delta_r_h12 < 3.6), False)
    # delta_r_h13_mask = ak.fill_none((events.delta_r_h13 < 3.6) & (events.delta_r_h13 > 2.0), False)
    # delta_r_h23_mask = ak.fill_none((events.delta_r_h23 < 3.4) & (events.delta_r_h23 > 2.0), False)
    # h3_mass_mask = ak.fill_none((events.h3_mass < 125.0) & (events.h3_mass > 50.0), False)
    # m_3b2tau_mask = ak.fill_none(events.m_3b2tau > 300.0, False)
    # m_3b2tau_pt_mask = ak.fill_none(events.m_3b2tau_pt > 300.0, False)
    # mhhh_mask = ak.fill_none(events.mhhh > 400.0, False)