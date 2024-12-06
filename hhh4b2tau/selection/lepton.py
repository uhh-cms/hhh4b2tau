# coding: utf-8

"""
Lepton selection methods.
"""

from __future__ import annotations

import law

from collections import defaultdict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.util import sorted_indices_from_mask
from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.util import maybe_import

from hhh4b2tau.production.example import cutflow_features

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={f"Tau.{var}" 
          for var in ("pt", "eta", "phi", "mass")},
)
def tau_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    tau_mask = (events.Tau.pt >= 20.0) & (abs(events.Tau.eta) < 2.5)
    tau_sel1 = ak.sum(tau_mask, axis=1) >= 1
    tau_sel2 = ak.sum(tau_mask, axis=1) >= 2


    return events, SelectionResult(
        steps={
            "one_tau": tau_sel1,
            "two_tau": tau_sel2
        },
        objects={
            "Tau": {
                "Tau": tau_mask,
            },
        },
    )