# coding: utf-8

"""
Tau scale factor production.
"""

import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view, layout_ak_array


ak = maybe_import("awkward")
np = maybe_import("numpy")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses={
        # custom columns created upstream, probably by a selector
        "single_triggered", "cross_triggered",
        # nano columns
        "Tau.pt", "Tau.eta", "Tau.genPartFlav", "Tau.decayMode",
    },
    produces={
        "tau_weight",
    } | {
        f"tau_weight_{unc}_{direction}"
        for direction in ["up", "down"]
        for unc in [
            "jet_dm0", "jet_dm1", "jet_dm10", "e_barrel", "e_endcap",
            "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
        ]
    },
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_tau_file=(lambda self, external_files: external_files.tau_sf),
    # function to determine the tau tagger name
    get_tau_tagger=(lambda self: self.config_inst.x.tau_tagger),
)
def tau_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Producer for tau ID weights. Requires an external file in the config under ``tau_sf``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "tau_sf": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/TAU/2017_UL/tau.json.gz",  # noqa
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files.

    The name of the tagger should be given as an auxiliary entry in the config:

    .. code-block:: python

        cfg.x.tau_tagger = "DeepTau2017v2p1"

    It is used to extract correction set names such as "DeepTau2017v2p1VSjet". *get_tau_tagger* can
    be adapted in a subclass in case it is stored differently in the config.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """
    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Tau.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayMode, axis=1)
    match = flat_np_view(events.Tau.genPartFlav, axis=1)

    # define channel / trigger dependent masks
    single_triggered = events.single_triggered
    cross_triggered = events.cross_triggered
    dm_mask = (
        (events.Tau.decayMode == 0) |
        (events.Tau.decayMode == 1) |
        (events.Tau.decayMode == 10) |
        (events.Tau.decayMode == 11)
    )

    #
    # compute nominal ID weights
    #

    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    wp_config = self.config_inst.x.tau_trigger_working_points

    # helpers to create corrector arguments
    if self.id_vs_jet_corrector.version == 0:
        # pt, dm, genmatch, jet wp, syst, sf type
        tau_args = lambda mask, syst: (pt[mask], dm[mask], match[mask], wp_config.id_vs_jet_v0, syst, "dm")
    elif self.id_vs_jet_corrector.version in (1, 2, 3):
        # pt, dm, genmatch, jet wp, e wp, syst, sf type
        tau_args = lambda mask, syst: (pt[mask], dm[mask], match[mask], *wp_config.id_vs_jet_gv0, syst, "dm")
    else:
        raise NotImplementedError

    if self.id_vs_e_corrector.version == 0:
        e_args = lambda mask, wp, syst: (abseta[mask], match[mask], wp, syst)
    elif self.id_vs_e_corrector.version in (1,):
        e_args = lambda mask, wp, syst: (abseta[mask], dm[mask], match[mask], wp, syst)
    else:
        raise NotImplementedError

    mu_args = lambda mask, wp, syst: (abseta[mask], match[mask], wp, syst)

    # genuine taus
    tau_mask = flat_np_view(dm_mask & (events.Tau.genPartFlav == 5), axis=1)
    sf_nom[tau_mask] = self.id_vs_jet_corrector(*tau_args(tau_mask, "nom"))

    # electrons faking taus
    e_mask = ((events.Tau.genPartFlav == 1) | (events.Tau.genPartFlav == 3))
    if self.config_inst.campaign.x.run == 3:
        e_mask = e_mask & (events.Tau.decayMode != 5) & (events.Tau.decayMode != 6)
    e_single_mask = flat_np_view((e_mask & single_triggered), axis=1)
    e_cross_mask = flat_np_view((e_mask & cross_triggered), axis=1)
    sf_nom[e_single_mask] = self.id_vs_e_corrector(*e_args(e_single_mask, wp_config.id_vs_e_single, "nom"))
    sf_nom[e_cross_mask] = self.id_vs_e_corrector(*e_args(e_cross_mask, wp_config.id_vs_e_cross, "nom"))

    # muons faking taus
    mu_mask = ((events.Tau.genPartFlav == 2) | (events.Tau.genPartFlav == 4))
    mu_single_mask = flat_np_view((mu_mask & single_triggered), axis=1)
    mu_cross_mask = flat_np_view((mu_mask & cross_triggered), axis=1)
    sf_nom[mu_single_mask] = self.id_vs_mu_corrector(*mu_args(mu_single_mask, wp_config.id_vs_mu_single, "nom"))
    sf_nom[mu_cross_mask] = self.id_vs_mu_corrector(*mu_args(mu_cross_mask, wp_config.id_vs_mu_cross, "nom"))

    # create and store weights
    events = set_ak_column_f32(events, "tau_weight", reduce_mul(sf_nom))

    #
    # compute varied ID weights
    #

    for direction in ["up", "down"]:
        # genuine taus -> split into decay modes
        sf_tau_dm0 = sf_nom.copy()
        sf_tau_dm1 = sf_nom.copy()
        sf_tau_dm10 = sf_nom.copy()
        tau_dm0_mask = tau_mask & (dm == 0)
        tau_dm1_mask = tau_mask & (dm == 1)
        tau_dm10_mask = tau_mask & ((dm == 10) | (dm == 11))
        sf_tau_dm0[tau_dm0_mask] = self.id_vs_jet_corrector(*tau_args(tau_dm0_mask, direction))
        sf_tau_dm1[tau_dm1_mask] = self.id_vs_jet_corrector(*tau_args(tau_dm1_mask, direction))
        sf_tau_dm10[tau_dm10_mask] = self.id_vs_jet_corrector(*tau_args(tau_dm10_mask, direction))
        events = set_ak_column_f32(events, f"tau_weight_jet_dm0_{direction}", reduce_mul(sf_tau_dm0))
        events = set_ak_column_f32(events, f"tau_weight_jet_dm1_{direction}", reduce_mul(sf_tau_dm1))
        events = set_ak_column_f32(events, f"tau_weight_jet_dm10_{direction}", reduce_mul(sf_tau_dm10))

        # electron fakes -> split into 2 eta regions
        for region, region_mask in [
            ("barrel", (abseta < 1.5)),
            ("endcap", (abseta >= 1.5)),
        ]:
            sf_e = sf_nom.copy()
            e_single_region_mask = e_single_mask & region_mask
            e_cross_region_mask = e_cross_mask & region_mask
            sf_e[e_single_region_mask] = self.id_vs_e_corrector(
                *e_args(e_single_region_mask, wp_config.id_vs_e_single, direction),
            )
            sf_e[e_cross_region_mask] = self.id_vs_e_corrector(
                *e_args(e_cross_region_mask, wp_config.id_vs_e_cross, direction),
            )
            events = set_ak_column_f32(events, f"tau_weight_e_{region}_{direction}", reduce_mul(sf_e))

        # muon fakes -> split into 5 eta regions
        for region, region_mask in [
            ("0p0To0p4", (abseta < 0.4)),
            ("0p4To0p8", ((abseta >= 0.4) & (abseta < 0.8))),
            ("0p8To1p2", ((abseta >= 0.8) & (abseta < 1.2))),
            ("1p2To1p7", ((abseta >= 1.2) & (abseta < 1.7))),
            ("1p7To2p3", (abseta >= 1.7)),
        ]:
            sf_mu = sf_nom.copy()
            mu_single_region_mask = mu_single_mask & region_mask
            mu_cross_region_mask = mu_cross_mask & region_mask
            sf_mu[mu_single_region_mask] = self.id_vs_mu_corrector(
                *mu_args(mu_single_region_mask, wp_config.id_vs_mu_single, direction),
            )
            sf_mu[mu_cross_region_mask] = self.id_vs_mu_corrector(
                *mu_args(mu_cross_region_mask, wp_config.id_vs_mu_cross, direction),
            )
            events = set_ak_column_f32(events, f"tau_weight_mu_{region}_{direction}", reduce_mul(sf_mu))

    return events


@tau_weights.requires
def tau_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@tau_weights.setup
def tau_weights_setup(self: Producer, reqs: dict, inputs: dict, reader_targets: InsertableDict) -> None:
    bundle = reqs["external_files"]

    # create the trigger and id correctors
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_tau_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    tagger_name = self.get_tau_tagger()
    self.id_vs_jet_corrector = correction_set[f"{tagger_name}VSjet"]
    self.id_vs_e_corrector = correction_set[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector = correction_set[f"{tagger_name}VSmu"]

    # check versions
    assert self.id_vs_jet_corrector.version in (0, 1, 2, 3)
    assert self.id_vs_e_corrector.version in (0, 1)
    assert self.id_vs_mu_corrector.version in (0, 1)


@producer(
    uses={
        "channel_id", "single_triggered", "cross_triggered",
        "Tau.pt", "Tau.decayMode",
    },
    produces={
        "tau_trigger_weight",
    } | {
        f"tau_trigger_weight_{ch}_{direction}"
        for direction in ["up", "down"]
        for ch in ["etau", "mutau", "tautau"]  # TODO: add tautauvbf when existing
    },
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_tau_file=(lambda self, external_files: external_files.tau_sf),
)
def trigger_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Producer for trigger scale factors derived by the TAU POG. Requires an external file in the
    config under ``tau_sf``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "tau_sf": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/TAU/2017_UL/tau.json.gz",  # noqa
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files. A correction set named ``"tau_trigger"`` is extracted from it.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """
    # get channels from the config
    ch_etau = self.config_inst.get_channel("etau")
    ch_mutau = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Tau.pt), axis=1, mask_identity=False)

    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    dm = flat_np_view(events.Tau.decayMode, axis=1)

    #
    # compute nominal trigger weight
    #

    # define channel / trigger dependent masks
    channel_id = events.channel_id
    single_triggered = events.single_triggered
    dm_mask = (
        (events.Tau.decayMode == 0) |
        (events.Tau.decayMode == 1) |
        (events.Tau.decayMode == 10) |
        (events.Tau.decayMode == 11)
    )
    tautau_mask = flat_np_view(
        dm_mask & (events.Tau.pt >= 40.0) & (channel_id == ch_tautau.id),
        axis=1,
    )
    # not existing yet
    # tautauvbf_mask = flat_np_view(dm_mask & (channel_id == ch_tautau.id), axis=1)
    etau_mask = flat_np_view(
        dm_mask & (channel_id == ch_etau.id) & single_triggered & (events.Tau.pt >= 25.0),
        axis=1,
    )
    mutau_mask = flat_np_view(
        dm_mask & (channel_id == ch_mutau.id) & single_triggered & (events.Tau.pt >= 25.0),
        axis=1,
    )

    # start with flat ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    wp_config = self.config_inst.x.tau_trigger_working_points
    eval_args = lambda mask, ch, syst: (pt[mask], dm[mask], ch, wp_config.trigger_corr, "sf", syst)
    sf_nom[etau_mask] = self.trigger_corrector(*eval_args(etau_mask, "etau", "nom"))
    sf_nom[mutau_mask] = self.trigger_corrector(*eval_args(mutau_mask, "mutau", "nom"))
    sf_nom[tautau_mask] = self.trigger_corrector(*eval_args(tautau_mask, "ditau", "nom"))

    # create and store weights
    events = set_ak_column_f32(events, "tau_trigger_weight", reduce_mul(sf_nom))

    #
    # compute varied trigger weights
    #

    for direction in ["up", "down"]:
        for ch, ch_corr, mask in [
            ("etau", "etau", etau_mask),
            ("mutau", "mutau", mutau_mask),
            ("tautau", "ditau", tautau_mask),
            # ("tautauvbf", "ditauvbf", tautauvbf_mask),
        ]:
            sf_unc = sf_nom.copy()
            sf_unc[mask] = self.trigger_corrector(*eval_args(mask, ch_corr, direction))
            events = set_ak_column_f32(events, f"tau_trigger_weight_{ch}_{direction}", reduce_mul(sf_unc))

    return events


@trigger_weights.requires
def trigger_weights_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@trigger_weights.setup
def trigger_weights_setup(self: Producer, reqs: dict, inputs: dict, reader_targets: InsertableDict) -> None:
    bundle = reqs["external_files"]

    # create the trigger and id correctors
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    # load the correction set
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_tau_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    self.trigger_corrector = correction_set["tau_trigger"]

    # check versions
    assert self.trigger_corrector.version in [0, 1]
