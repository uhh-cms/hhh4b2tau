from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.production.processes import process_ids
from hhh4b2tau.production.processes import process_ids_genmatched_higgs

from hhh4b2tau.production.higgs_reco import higgs_reco_mass_diff, higgs_reco_chi2

ak = maybe_import("awkward")

@producer(
    uses={higgs_reco_mass_diff, process_ids_genmatched_higgs},
)
def produce_genmatched_procids_mass_diff(self, events: ak.Array, **kwargs):

    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[higgs_reco_mass_diff](events, **kwargs)

        events = self[process_ids_genmatched_higgs](events, **kwargs)

    return events

@produce_genmatched_procids_mass_diff.init
def produce_genmatched_procids_mass_diff_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column

@producer(
    uses={higgs_reco_chi2, process_ids_genmatched_higgs},
)
def produce_genmatched_procids_chi2(self, events: ak.Array, **kwargs):

    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[higgs_reco_chi2](events, **kwargs)

        events = self[process_ids_genmatched_higgs](events, **kwargs)

    return events

@produce_genmatched_procids_mass_diff.init
def produce_genmatched_procids_mass_diff_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column







