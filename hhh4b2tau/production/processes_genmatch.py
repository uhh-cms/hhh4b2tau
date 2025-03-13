from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.production.processes import process_ids

from hhh4b2tau.production.higgs_reco import higgs_reco_mass_diff

ak = maybe_import("awkward")

@producer(
    uses={higgs_reco_mass_diff, "Gen_Matched_H{1,2}_idx"},
)
def produce_genmatched_procids_mass_diff(self, events: ak.Array, **kwargs):

    if self.dataset_inst.has_tag("hhh"):
        # do stuff
        events = self[higgs_reco_mass_diff](events, **kwargs)

        from IPython import embed
        embed(header=f"in {self.__class__.__name__}")

    return events

@produce_genmatched_procids_mass_diff.init
def produce_genmatched_procids_mass_diff_init(self):
    self.out_column = process_ids.produces
    self.uses |= self.out_column
    self.produces |= self.out_column