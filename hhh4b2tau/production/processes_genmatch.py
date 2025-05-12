from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.production.processes import process_ids
from columnflow.production.categories import category_ids
from hhh4b2tau.production.processes import process_ids_genmatched_higgs

from hhh4b2tau.production.higgs_reco import higgs_reco_mds, higgs_reco_chi2, higgs_reco_dhh

ak = maybe_import("awkward")

@producer(
    uses={higgs_reco_mds, process_ids_genmatched_higgs, category_ids,},
    produces={higgs_reco_mds, category_ids,},
)
def produce_genmatched_procids_mds(self, events: ak.Array, **kwargs):

    events = self[higgs_reco_mds](events, **kwargs)
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[process_ids_genmatched_higgs](events, **kwargs)

    return events

@produce_genmatched_procids_mds.init
def produce_genmatched_procids_mds_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column

@producer(
    uses={higgs_reco_chi2, process_ids_genmatched_higgs, category_ids,},
    produces={higgs_reco_chi2, category_ids,},
)
def produce_genmatched_procids_chi2(self, events: ak.Array, **kwargs):
    events = self[higgs_reco_chi2](events, **kwargs)
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[process_ids_genmatched_higgs](events, **kwargs)

    return events

@produce_genmatched_procids_chi2.init
def produce_genmatched_procids_chi2_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column


@producer(
    uses={higgs_reco_dhh, process_ids_genmatched_higgs, category_ids,},
    produces={higgs_reco_dhh, category_ids,},
)
def produce_genmatched_procids_dhh(self, events: ak.Array, **kwargs):
    events = self[higgs_reco_dhh](events, **kwargs)
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[process_ids_genmatched_higgs](events, **kwargs)

    return events

@produce_genmatched_procids_dhh.init
def produce_genmatched_procids_dhh_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column







