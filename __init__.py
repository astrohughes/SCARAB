"""
astroframes — fast, correct aerospace reference frames and time systems.

Rust core with Python bindings via PyO3/maturin.

Usage:
    from astroframes import Epoch, TimeScale, Frame, CartesianState

    t = Epoch.from_utc(2024, 6, 15, 12, 0, 0.0)
    state = CartesianState([6778.0, 0, 0], [0, 7.5, 0], t, Frame.GCRF)
"""
from ._astroframes import (
    CartesianState,
    Epoch,
    Frame,
    TimeScale,
)

__all__ = [
    "Epoch",
    "TimeScale",
    "Frame",
    "CartesianState",
]

__version__ = "0.1.0"
