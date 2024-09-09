# coding: utf-8

"""
Configuration of the HHH4b2tau analysis.
"""

from __future__ import annotations

import functools
import importlib

import law
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import


ak = maybe_import("awkward")


#
# the main analysis object
#

analysis_hhh4b2tau = ana = od.Analysis(
    name="analysis_hhh4b2tau",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = [
    "$CF_BASE/sandboxes/cf.sh",
    "$CF_BASE/sandboxes/venv_columnar.sh",
]

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    # "$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}


#
# setup configs
#

# an example config is setup below, based on cms NanoAOD v9 for Run2 2017, focussing on
# ttbar and single top MCs, plus single muon data
# update this config or add additional ones to accomodate the needs of your analysis
from hhh4b2tau.config.config_run3_2022 import add_config

def add_lazy_config(
    campaign_module: str,
    campaign_attr: str,
    config_name: str,
    config_id: int,
    limit_dataset_files: int | None = None,
    **kwargs,
):
    def create_factory(
        config_id: int,
        config_name_postfix: str = "",
        limit_dataset_files: int | None = None,
    ):
        def factory(configs: od.UniqueObjectIndex):
            # import the campaign
            mod = importlib.import_module(campaign_module)
            campaign = getattr(mod, campaign_attr)

            return add_config(
                ana,
                campaign.copy(),
                config_name=config_name + config_name_postfix,
                config_id=config_id,
                limit_dataset_files=limit_dataset_files,
                **kwargs,
            )
        return factory

    ana.configs.add_lazy_factory(config_name, create_factory(config_id))
    ana.configs.add_lazy_factory(f"{config_name}_limited", create_factory(config_id + 200, "_limited", limit_dataset_files))


# add_config(
#     analysis=ana,
#     campaign=campaign_run3_2022_preEE_nano_v12.copy(),
#     config_name=campaign_run3_2022_preEE_nano_v12.name,
#     config_id=30,  # 3 2022 0 (run3 year 2022 preEE) 
# )

# add_config(
#     analysis=ana,
#     campaign=campaign_run3_2022_preEE_nano_v12.copy(),
#     config_name=f"{campaign_run3_2022_preEE_nano_v12.name}_limited",
#     config_id=31,  # 3 2022 0 (run3 year 2022 preEE limited)
#     limit_dataset_files=2, 
# )
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_preEE_nano_uhh_v12",
    campaign_attr="campaign_run3_2022_preEE_nano_uhh_v12",
    config_name="run3_2022_preEE",
    config_id=5,
    limit_dataset_files=2,
)

add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_preEE_nano_v12",
    campaign_attr="campaign_run3_2022_preEE_nano_v12",
    config_name="run3_2022_preEE_nano_v12",
    config_id=6,
    limit_dataset_files=2,
)
