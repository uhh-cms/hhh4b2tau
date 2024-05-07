# coding: utf-8

"""
Configuration of the HHH4b2tau analysis.
"""

import functools

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

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    "$CF_BASE/sandboxes/cmssw_default.sh",
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

from cmsdb.campaigns.run3_2022_preEE_nano_v12 import campaign_run3_2022_preEE_nano_v12

from hhh4b2tau.config.config_run3_2022 import add_config

add_config(
    analysis=ana,
    campaign=campaign_run3_2022_preEE_nano_v12.copy(),
    config_name=campaign_run3_2022_preEE_nano_v12.name,
    config_id=30,  # 3 2022 0 (run3 year 2022 preEE) 
)

add_config(
    analysis=ana,
    campaign=campaign_run3_2022_preEE_nano_v12.copy(),
    config_name=f"{campaign_run3_2022_preEE_nano_v12.name}_limited",
    config_id=31,  # 3 2022 0 (run3 year 2022 preEE limited)
    limit_dataset_files=2, 
)
