# coding: utf-8

"""
Producers that determine the generator-level particles related to a top quark decay.
"""

from __future__ import annotations
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import (
    set_ak_column, remove_ak_column, attach_behavior, EMPTY_FLOAT, get_ak_routes, remove_ak_column,
    optional_column as optional
)
from columnflow.types import Sequence
import numpy as np

ak = maybe_import("awkward")

class _GenPartMatchBase(Producer):

    def __init__(self, *args, **kwargs):
        # first, call the init function of the super class (Producer)
        super().__init__(*args, **kwargs)
        # define variables that are needed for the gen matching
        self.mothers: tuple[str] = tuple()
        self.children: tuple[str] = tuple()
        self.variables: tuple[str] = ('pt', 'eta', 'phi', 'mass', 'pdgId', )

    def init_func(self):
        # to perform a gen matching, we need information about the GenPartons
        # therefore, request all available information
        self.uses=(
            {"GenPart.*"}
        )
        # derived classes should produce the following output columns
        # note that this needs to be set explicitly when calling get_decay_idx!
        self.produces=(
            {
                optional(f"gen_{mother}_to_{child}.{var}")
                for mother in self.mothers
                for child in self.children
                for var in self.variables
            } |
            {   
                optional(f"gen_{child}.{var}")
                for child in self.children
                for var in self.variables
            } |
            {   
                optional(f"gen_tth_{child}.{var}")
                for child in self.children
                for var in self.variables
            } |
            {   
                optional(f"gen_tth_{mother}_to_{child}.{var}")
                for mother in self.mothers
                for child in self.children
                for var in self.variables
            }
        )

    def get_decay_idx(
        self,
        events: ak.Array,
        mother_id: int,
        children_id: int,
        children_output_name: str | None = None,
        mother_output_name: str | None = None,
        mother_gen_flags: list | None = None,
        children_gen_flags: list | None = None,
    ):
        abs_id = abs(events.GenPart.pdgId)
        try:
            if not mother_gen_flags:
                mother_gen_flags = ["isLastCopy", "fromHardProcess"]
            if not children_gen_flags:
                children_gen_flags = ["isFirstCopy", "fromHardProcess"]
            # obtain indices and GenParticles for mothers that match input 'mother_id'
            mask_mother_id = abs_id == mother_id
            mask_gen_stati = ak.mask(events.GenPart, mask_mother_id).hasFlags(*mother_gen_flags)
            mask_gen_stati = ak.fill_none(mask_gen_stati, False, axis=-1)
            mother_idx = ak.local_index(events.GenPart.pt, axis=-1)[mask_gen_stati]
            mothers = events.GenPart[mother_idx]

            # sort mothers for pt
            sorted_mothers_idx = ak.argsort(mothers.pt, axis=-1, ascending=False)
            mother_idx = mother_idx[sorted_mothers_idx]
            mothers = mothers[sorted_mothers_idx]

            # get corresponding decay products, aka children
            children = mothers.distinctChildren
            # make sure you only consider real children
            children_gen_stati_mask = children.hasFlags(*children_gen_flags)
            children = children[children_gen_stati_mask]
            abs_children_id = abs(children.pdgId)

            children_mask = abs_children_id == children_id
            any_relevant_children_mask = ak.any(children_mask, axis=-1)
            relevant_mother_idx = mother_idx[any_relevant_children_mask]
            relevant_mothers = events.GenPart[relevant_mother_idx]

            #update children masks
            children = children[any_relevant_children_mask]
            children_mask = children_mask[any_relevant_children_mask]
            
            relevant_children = children[children_mask]
            relevant_children_idx = ak.local_index(relevant_children.pt, axis=-1)

            sorted_children_idx = ak.argsort(relevant_children.pt, axis=-1, ascending=False)
            relevant_children_idx = relevant_children_idx[sorted_children_idx]

        except Exception as e:
            from IPython import embed
            embed(header=str(e))
            # raise e

        if children_output_name:
            for var in ('pt', 'eta', 'phi', 'mass', 'pdgId'):
                events = set_ak_column(events, f'{children_output_name}.{var}', getattr(relevant_children, var))
        if mother_output_name:
            for var in ('pt', 'eta', 'phi', 'mass', 'pdgId'):
                events = set_ak_column(events, f'{mother_output_name}.{var}', getattr(relevant_mothers, var))

        return events, relevant_mother_idx, relevant_children


@_GenPartMatchBase.producer(
    mothers = ('h', 'tau', ),
    children = ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
)
def gen_higgs_decay_products(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Creates a new ragged column "gen_higgs_decay" with one element per hard higgs boson. Each element is
    a GenParticleArray with five or more objects in a distinct order: higgs boson, bottom quark, tau lepton,
    W boson, down-type quark or charged lepton, up-type quark or neutrino, and any additional decay
    produces of the W boson (if any, then most likely photon radiations). Per event, the structure
    will be similar to:

    .. code-block:: python

        [
            # event 1
            [
                # top 1
                [t1, b1, W1, q1/l, q2/n(, additional_w_decay_products)],
                # top 2
                [...],
            ],
            # event 2
            ...
        ]
    """
    # from IPython import embed; embed()
    events, h_b_idx, h_b_particles = self.get_decay_idx(
        events,
        mother_id=25,
        children_id=5,
        children_output_name="gen_b",
        mother_output_name="gen_h_to_b",
    )
    events, h_tau_idx, h_tau_particles = self.get_decay_idx(
        events,
        mother_id=25,
        children_id=15,
        children_output_name="gen_tau",
        mother_output_name="gen_h_to_tau",
    )

    events, tau_taunu_idx, tau_taunu_particles = self.get_decay_idx(
        events,
        mother_id=15,
        children_id=16,
        children_output_name="gen_taunu",
        mother_output_name="gen_tau_to_taunu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_electron_idx, tau_electron_particles = self.get_decay_idx(
        events,
        mother_id=15,
        children_id=11,
        children_output_name="gen_electron",
        mother_output_name="gen_tau_to_electron",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_enu_idx, tau_enu_particles = self.get_decay_idx(
        events,
        mother_id=15,
        children_id=12,
        children_output_name="gen_enu",
        mother_output_name="gen_tau_to_enu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_muon_idx, tau_muon_particles = self.get_decay_idx(
        events,
        mother_id=15,
        children_id=13,
        children_output_name="gen_muon",
        mother_output_name="gen_tau_to_muon",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_munu_idx, tau_munu_particles = self.get_decay_idx(
        events,
        mother_id=15,
        children_id=14,
        children_output_name="gen_munu",
        mother_output_name="gen_tau_to_munu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )
    # from IPython import embed
    # embed(header="in gen_higgs_decay_products after W identification")
    return events
  
# Access decay products for ttH channel
@_GenPartMatchBase.producer(
        mothers=('w', 'h', 't'),
        children=('tau', 'taunu', 'b1', 'b2')
)
def gen_tth_decay_products(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events, h_b_idx, h_b_particles = self.get_decay_idx(
        events,
        mother_id=25,
        children_id=5,
        children_output_name="gen_tth_b1",
        mother_output_name="gen_tth_h_to_b",
    )

    events, t_b_idx, t_b_particles = self.get_decay_idx(
        events,
        mother_id=6,
        children_id=5,
        children_output_name="gen_tth_b2",
        mother_output_name="gen_tth_t_to_b",
    )

    events, w_tau_idx, w_tau_particles = self.get_decay_idx(
        events,
        mother_id=24,
        children_id=15,
        children_output_name="gen_tth_tau",
        mother_output_name="gen_tth_w_to_tau",
    )

    events, h_b_idx, h_b_particles = self.get_decay_idx(
        events,
        mother_id=24,
        children_id=16,
        children_output_name="gen_tth_taunu",
        mother_output_name="gen_tth_w_to_taunu",
    )
    # from IPython import embed; embed(header='in gen_tth')
    return events