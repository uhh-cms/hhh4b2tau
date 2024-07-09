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


@producer(
    uses=(
        {"GenPart.*"}
        ),
    # produces={"gen_higgs_decay.*", "z_pion_neg", "z_pion_pos", "z_kaon_neg", "z_kaon_pos", "pion_neg.*", "pion_pos", "pion_neg_E", "pion_pos_E"},
    produces=({

        optional(f"gen_{mother}_to_{child}.{var}")
        # optional(f"gen_{child}_from_{mother}.{var}")
        for mother in ('h', 'tau', )
        for child in ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )
    } |
    {   optional(f"gen_{child}_from_{mother}.{var}")
        for mother in ('h', 'tau', )
        for child in ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )}
    
    ),
)
def gen_higgs_decay_products(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Creates a new ragged column "gen_higgs_decay" with one element per hard higgs boson. Each element is
    a GenParticleArray with five or more objects in a distinct order: higgs boson, bottom quark, tau lepton,
    W boson, down-type quark or charged lepton, up-type quark or neutrino, and any additional decay
    produces of the W boson (if any, then most likly photon radiations). Per event, the structure
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
    # n = 19362

    # # find hard h boson
    
    # h = events.GenPart[abs_id == 25]
    # h = h[h.hasFlags("isFirstCopy", "fromHardProcess")]
    # h = ak.drop_none(h, behavior=h.behavior)

    # # distinct higgs boson children (b's and tau's)
    # h_children = h.distinctChildrenDeep
    # abs_children_id = abs(h_children.pdgId)
    abs_id = abs(events.GenPart.pdgId)
    def get_decay_idx(
        events: ak.Array,
        mother_id: int,
        children_id: int,
        children_output_name: str | None = None,
        mother_output_name: str | None = None,
        mother_gen_flags: list | None = None,
        children_gen_flags: list | None = None,
    ):
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

    events, h_b_idx, h_b_particles = get_decay_idx(
        events,
        mother_id=25,
        children_id=5,
        children_output_name="gen_b_from_h",
        mother_output_name="gen_h_to_b",
    )
    events, h_tau_idx, h_tau_particles = get_decay_idx(
        events,
        mother_id=25,
        children_id=15,
        children_output_name="gen_tau_from_h",
        mother_output_name="gen_h_to_tau",
    )
    # from IPython import embed; embed()
    events, tau_taunu_idx, tau_taunu_particles = get_decay_idx(
        events,
        mother_id=15,
        children_id=16,
        children_output_name="gen_taunu_from_tau",
        mother_output_name="gen_tau_to_taunu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_electron_idx, tau_electron_particles = get_decay_idx(
        events,
        mother_id=15,
        children_id=11,
        children_output_name="gen_electron_from_tau",
        mother_output_name="gen_tau_to_electron",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_enu_idx, tau_enu_particles = get_decay_idx(
        events,
        mother_id=15,
        children_id=12,
        children_output_name="gen_enu_from_tau",
        mother_output_name="gen_tau_to_enu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_muon_idx, tau_muon_particles = get_decay_idx(
        events,
        mother_id=15,
        children_id=13,
        children_output_name="gen_muon_from_tau",
        mother_output_name="gen_tau_to_muon",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )

    events, tau_munu_idx, tau_munu_particles = get_decay_idx(
        events,
        mother_id=15,
        children_id=14,
        children_output_name="gen_munu_from_tau",
        mother_output_name="gen_tau_to_munu",
        children_gen_flags=["isFirstCopy", "isPromptTauDecayProduct"],
    )
    
    
    # from IPython import embed
    # embed(header="after get_decay_idx")


    # # get b's
    # b = events.GenPart[abs_id == 5]
    # b = b[b.hasFlags("isFirstCopy", "fromHardProcess")]
    # # remove optional (first remove nones, then update the type)
    # b = ak.drop_none(b, behavior=b.behavior)

    # # get tau's
    # tau = events.GenPart[abs_id == 15]
    # tau = tau[tau.hasFlags("isLastCopy", "fromHardProcess")]
    # # remove optional
    # tau = ak.drop_none(tau, behavior=tau.behavior)
    
    # # distinct tau children
    # tau_children = tau.distinctChildren
    # abs_tau_children_id = abs(tau_children.pdgId)
    # tau_children = ak.drop_none(tau_children, behavior=tau_children.behavior)

    # # separate tau-positive, tau-negative
    # tau_pos = events.GenPart[events.GenPart.pdgId == -15]
    # tau_pos = tau_pos[tau_pos.hasFlags("isLastCopy", "fromHardProcess")]
    # tau_pos = ak.drop_none(tau_pos, behavior=tau_pos.behavior)
    
    # tau_neg = events.GenPart[events.GenPart.pdgId == 15]
    # tau_neg = tau_neg[tau_neg.hasFlags("isLastCopy", "fromHardProcess")]
    # tau_neg = ak.drop_none(tau_neg, behavior=tau_neg.behavior)
    
    # tau_children_pos = tau_pos.distinctChildren
    # abs_tau_children_pos_id = abs(tau_children_pos.pdgId)
    # tau_children_pos = ak.drop_none(tau_children_pos, behavior=tau_children_pos.behavior)
    
    # tau_children_neg = tau_neg.distinctChildren
    # abs_tau_children_neg_id = abs(tau_children_neg.pdgId)
    # tau_children_neg = ak.drop_none(tau_children_neg, behavior=tau_children_neg.behavior)

    # h_tau = tau.parent.parent
    # h_tau = ak.drop_none(h_tau, behavior=h_tau.behavior)
   
    # h_b = b.parent
    # h_b = ak.drop_none(h_b, behavior=h_b.behavior)
    
    # # # get nu's
    # # nu = ak.firsts(tau_children[abs_tau_children_id == 16], axis=2)
    # # # remove optional
    # # nu = ak.drop_none(nu, behavior=nu.behavior)

    # # get nu's
    # nu = events.GenPart[abs_id == 16]
    # nu = nu[nu.hasFlags("isFirstCopy", "isDirectTauDecayProduct")]
    # # remove optional (first remove nones, then update the type)
    # nu = ak.drop_none(nu, behavior=nu.behavior)

    # # from IPython import embed
    # # embed(header="in gen_higgs_decay_products after neutrino identification")
    
    # # # get w's
    # # w = ak.firsts(tau_children[abs_tau_children_id == 24], axis=2)
    # # # remove optional
    # # w = ak.drop_none(w, behavior=w.behavior)


    # # get w's
    # w = events.GenPart[abs_id == 24]
    # w = w[w.hasFlags("isFirstCopy", "fromHardProcess")]
    # # remove optional (first remove nones, then update the type)
    # w = ak.drop_none(w, behavior=w.behavior)

    # from IPython import embed
    # embed(header="in gen_higgs_decay_products after W identification")
    


#     # decays might also be effective into a variable number of mesons -> add them
#     w_children = tau_children[abs_tau_children_id != 16]
#     # remove optional
#     w_children = ak.drop_none(w_children, behavior=w_children.behavior)
    
#     # temporarily here: sum children and make them GenParticles
#     w_sum = w_children.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_sum = set_ak_column(w_sum, c, getattr(w_sum, c))
#     for c in ["t", "x", "y", "z"]:
#         w_sum = remove_ak_column(w_sum, c)
#     tau_sign = tau.pdgId // abs(tau.pdgId)
#     w_sum = set_ak_column(w_sum, "pdgId", tau.pdgId - 39 * tau_sign)
#     w_sum = attach_behavior(w_sum, "GenParticle")


#    # visible decay of w
#     w_visible = tau_children[(abs_tau_children_id != 12) & (abs_tau_children_id != 14) & (abs_tau_children_id != 16)]
#     # remove optional
#     w_visible = ak.drop_none(w_visible, behavior=w_visible.behavior)

#     w_visible_sum = w_visible.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_visible_sum = set_ak_column(w_visible_sum, c, getattr(w_visible_sum, c))
#     for c in ["t", "x", "y", "z"]:
#         w_visible_sum = remove_ak_column(w_visible_sum, c)
#     tau_sign = tau.pdgId // abs(tau.pdgId)
#     w_visible_sum = set_ak_column(w_visible_sum, "pdgId", tau.pdgId - 39 * tau_sign)
#     w_visible_sum = attach_behavior(w_visible_sum, "GenParticle")



#     w_vis_had_ne = tau_children_neg[
#         ((abs(tau_children_neg.pdgId) < 11) | (abs(tau_children_neg.pdgId) > 16)) &
#         (tau_children_neg.pdgId != 22) &
#         ~tau_children_neg.hasFlags("isPrompt")
#     ]
#     # visible hadronic w decay    
#     w_vis_had_neg = w_vis_had_ne.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_vis_had_neg = set_ak_column(w_vis_had_neg, c, getattr(w_vis_had_neg, c))
#     for c in ["t", "x", "y", "z"]:
#         w_vis_had_neg = remove_ak_column(w_vis_had_neg, c)
#     tau_neg_sign = tau_neg.pdgId // abs(tau_neg.pdgId)
#     w_vis_had_neg = set_ak_column(w_vis_had_neg, "pdgId", tau_neg.pdgId - 39 * tau_neg_sign)
#     w_vis_had_neg = attach_behavior(w_vis_had_neg, "GenParticle")


#     w_vis_had_po = tau_children_pos[
#         ((abs(tau_children_pos.pdgId) < 11) | (abs(tau_children_pos.pdgId) > 16)) &
#         (tau_children_pos.pdgId != 22) &
#         ~tau_children_pos.hasFlags("isPrompt")
#     ]
#     w_vis_had_pos = w_vis_had_po.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_vis_had_pos = set_ak_column(w_vis_had_pos, c, getattr(w_vis_had_pos, c))
#     for c in ["t", "x", "y", "z"]:
#         w_vis_had_pos = remove_ak_column(w_vis_had_pos, c)
#     tau_pos_sign = tau_pos.pdgId // abs(tau_pos.pdgId)
#     w_vis_had_pos = set_ak_column(w_vis_had_pos, "pdgId", tau_pos.pdgId - 39 * tau_pos_sign)
#     w_vis_had_pos = attach_behavior(w_vis_had_pos, "GenParticle")    

#     # ------------------------------
#     # visible leptonic w decay
#     w_vis_lep_po = tau_children_pos[
#         ((abs(tau_children_pos.pdgId) == 11) | (abs(tau_children_pos.pdgId) == 13) | (abs(tau_children_pos.pdgId) == 15)) &
#         ~tau_children_pos.hasFlags("isPrompt")
#     ]
    
#     w_vis_lep_pos = w_vis_lep_po.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_vis_lep_pos = set_ak_column(w_vis_lep_pos, c, getattr(w_vis_lep_pos, c))
#     for c in ["t", "x", "y", "z"]:
#         w_vis_lep_pos = remove_ak_column(w_vis_lep_pos, c)
#     w_vis_lep_pos = set_ak_column(w_vis_lep_pos, "pdgId", tau_pos.pdgId - 39 * tau_pos_sign)
#     w_vis_lep_pos = attach_behavior(w_vis_lep_pos, "GenParticle")


#     w_vis_lep_ne = tau_children_neg[
#         ((abs(tau_children_neg.pdgId) == 11) | (abs(tau_children_neg.pdgId) == 13) | (abs(tau_children_neg.pdgId) == 15)) &
#         ~tau_children_neg.hasFlags("isPrompt")
#     ]
    
#     w_vis_lep_neg = w_vis_lep_ne.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         w_vis_lep_neg = set_ak_column(w_vis_lep_neg, c, getattr(w_vis_lep_neg, c))
#     for c in ["t", "x", "y", "z"]:
#         w_vis_lep_neg = remove_ak_column(w_vis_lep_neg, c)
#     w_vis_lep_neg = set_ak_column(w_vis_lep_neg, "pdgId", tau_neg.pdgId - 39 * tau_neg_sign)
#     w_vis_lep_neg = attach_behavior(w_vis_lep_neg, "GenParticle")
    
#     # pion plus 
#     pi_pos_neg = tau_children_neg[(tau_children_neg.pdgId == 211) & ~tau_children_neg.hasFlags("isPrompt")] 
#     pion_pos_neg = pi_pos_neg.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_pos_neg = set_ak_column(pion_pos_neg, c, getattr(pion_pos_neg, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_pos_neg = remove_ak_column(pion_pos_neg, c)
#     pion_pos_neg = set_ak_column(pion_pos_neg, "pdgId", 211)
#     pion_pos_neg = attach_behavior(pion_pos_neg, "GenParticle")

#     pi_pos_pos = tau_children_pos[(tau_children_pos.pdgId == 211) & ~tau_children_pos.hasFlags("isPrompt")] 
#     pion_pos_pos = pi_pos_pos.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_pos_pos = set_ak_column(pion_pos_pos, c, getattr(pion_pos_pos, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_pos_pos = remove_ak_column(pion_pos_pos, c)
#     pion_pos_pos = set_ak_column(pion_pos_pos, "pdgId", 211)
#     pion_pos_pos = attach_behavior(pion_pos_pos, "GenParticle")
    
#     # pion minus
#     pi_neg_neg = tau_children_neg[(tau_children_neg.pdgId == -211) & ~tau_children_neg.hasFlags("isPrompt")] 
#     pion_neg_neg = pi_neg_neg.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_neg_neg = set_ak_column(pion_neg_neg, c, getattr(pion_neg_neg, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_neg_neg = remove_ak_column(pion_neg_neg, c)
#     pion_neg_neg = set_ak_column(pion_neg_neg, "pdgId", -211)
#     pion_neg_neg = attach_behavior(pion_neg_neg, "GenParticle")

#     pi_neg_pos = tau_children_pos[(tau_children_pos.pdgId == -211) & ~tau_children_pos.hasFlags("isPrompt")] 
#     pion_neg_pos = pi_neg_pos.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_neg_pos = set_ak_column(pion_neg_pos, c, getattr(pion_neg_pos, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_neg_pos = remove_ak_column(pion_neg_pos, c)
#     pion_neg_pos = set_ak_column(pion_neg_pos, "pdgId", -211)
#     pion_neg_pos = attach_behavior(pion_neg_pos, "GenParticle")
    

#     pion = tau_children_pos[(tau_children_pos.pdgId == 111) & ~tau_children_pos.hasFlags("isPrompt")]
#     pion_zero_pos = pion.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_zero_pos = set_ak_column(pion_zero_pos, c, getattr(pion_zero_pos, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_zero_pos = remove_ak_column(pion_zero_pos, c)
#     pion_zero_pos = set_ak_column(pion_zero_pos, "pdgId", 111)
#     pion_zero_pos = attach_behavior(pion_zero_pos, "GenParticle")

#     pion1 = tau_children_neg[(tau_children_neg.pdgId == 111) & ~tau_children_neg.hasFlags("isPrompt")]
#     pion_zero_neg = pion1.sum(axis=-1)
#     for c in ["pt", "eta", "phi", "mass"]:
#         pion_zero_neg = set_ak_column(pion_zero_neg, c, getattr(pion_zero_neg, c))
#     for c in ["t", "x", "y", "z"]:
#         pion_zero_neg = remove_ak_column(pion_zero_neg, c)
#     pion_zero_neg = set_ak_column(pion_zero_neg, "pdgId", 111)
#     pion_zero_neg = attach_behavior(pion_zero_neg, "GenParticle")
    


#     tau_neg_2c = tau_children_neg[ak.num(tau_children_neg, axis=2) == 2]
#     pion_neg = ak.flatten(tau_neg_2c[tau_neg_2c.pdgId == -211], axis=2)
#     # from IPython import embed
#     # embed()
#     kaon_neg = ak.flatten(tau_neg_2c[tau_neg_2c.pdgId == -321], axis=2)

#     tau_pos_2c = tau_children_pos[ak.num(tau_children_pos, axis=2) == 2]
#     pion_pos = ak.flatten(tau_pos_2c[tau_pos_2c.pdgId == 211], axis=2)
#     kaon_pos = ak.flatten(tau_pos_2c[tau_pos_2c.pdgId == 321], axis=2)

    return events
