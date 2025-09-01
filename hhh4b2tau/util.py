# coding: utf-8

"""
Helpful utils.
"""

from __future__ import annotations

__all__ = []

import law

from columnflow.types import Any
from columnflow.columnar_util import ArrayFunction, deferred_column
from columnflow.util import maybe_import

from typing import Hashable, Iterable, Callable

np = maybe_import("numpy")
ak = maybe_import("awkward")

_logger = law.logger.get_logger(__name__)

@deferred_column
def IF_NANO_V9(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V11(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 11 else None


@deferred_column
def IF_NANO_V12(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 12 else None

@deferred_column
def IF_NANO_V14(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 14 else None

@deferred_column
def IF_NANO_GE_V10(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None

@deferred_column
def IF_RUN_2(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.run == 2 else None


@deferred_column
def IF_RUN_3(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.run == 3 else None


@deferred_column
def IF_RUN_3_2022(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if (func.config_inst.campaign.x.run == 3 and func.config_inst.campaign.x.year == 2022) else None


@deferred_column
def IF_DATASET_HAS_LHE_WEIGHTS(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return None if func.dataset_inst.has_tag("no_lhe_weights") else self.get()


@deferred_column
def IF_DATASET_IS_DY(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return self.get() if func.dataset_inst.has_tag("is_dy") else None

@deferred_column
def IF_DATASET_IS_W_LNU(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return self.get() if func.dataset_inst.has_tag("w_lnu") else None

@deferred_column
def IF_DATASET_IS_TT(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return self.get() if func.dataset_inst.has_tag("is_ttbar") else None

def hash_events(arr: np.ndarray) -> np.ndarray:
    """
    Helper function to create a hash value from the event, run and luminosityBlock columns.
    The values are padded to specific lengths and concatenated to a single integer.
    """
    import awkward as ak

    def assert_value(arr: np.ndarray, field: str, max_value: int) -> None:
        """
        Helper function to check if a column does not exceed a maximum value.
        """
        digits = len(str(arr[field].to_numpy().max()))
        assert digits <= max_value, f"{field} digit count is {digits} and exceed max value {max_value}"

    max_digits_run = 6
    max_digits_luminosityBlock = 6
    max_digits_event = 8
    assert_value(arr, "run", max_digits_run)
    assert_value(arr, "luminosityBlock", max_digits_luminosityBlock)
    assert_value(arr, "event", max_digits_event)

    max_digits_hash = max_digits_event + max_digits_luminosityBlock + max_digits_run
    assert max_digits_hash <= 20, "sum of digits exceeds int64"

    # upcast to int64 to avoid overflow
    return (
        ak.values_astype(arr.run, np.int64) * 10**(max_digits_luminosityBlock + max_digits_event) +
        ak.values_astype(arr.luminosityBlock, np.int64) * 10**max_digits_event +
        ak.values_astype(arr.event, np.int64)
    )

def round_sig(
    value: int | float | np.number,
    sig: int = 4,
    convert: Callable | None = None,
) -> int | float | np.number:
    """
    Helper function to round number *value* on *sig* significant digits and
    optionally transform output to type *convert*
    """
    if not np.isfinite(value):
        # cannot round infinite
        _logger.warning("cannot round infinite number")
        return value

    from math import floor, log10

    def try_rounding(_value):
        try:
            n_digits = sig - int(floor(log10(abs(_value)))) - 1
            if convert in (int, np.int8, np.int16, np.int32, np.int64):
                # do not round on decimals when converting to integer
                n_digits = min(n_digits, 0)
            return round(_value, n_digits)
        except Exception:
            _logger.warning(f"Cannot round number {value} to {sig} significant digits. Number will not be rounded")
            return value

    # round first to not lose information from type conversion
    rounded_value = try_rounding(value)

    # convert number if "convert" is given
    if convert not in (None, False):
        try:
            rounded_value = convert(rounded_value)
        except Exception:
            _logger.warning(f"Cannot convert {rounded_value} to {convert.__name__}")
            return rounded_value

        # some types need rounding again after converting (e.g. np.float32 to float)
        rounded_value = try_rounding(rounded_value)

    return rounded_value
