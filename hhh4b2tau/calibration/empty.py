from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.met import met_phi
from columnflow.calibration.cms.jets import jec, jer
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column



np = maybe_import("numpy")
ak = maybe_import("awkward")


# derive calibrators to add settings
jec_nominal = jec.derive("jec_nominal", cls_dict={"uncertainty_sources": [], "data_only": True})
jec_full = jec.derive("jec_nominal", cls_dict={"mc_only": True})
@calibrator(
    uses=set(),
    produces=set(),
)
def empty(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    return events