"""magnetisme — modélisation numérique de matériaux magnétiques (Weiss, Ising 2D, Onsager, Wolff)."""
from .ising import HAS_NUMBA, Ising2D, mesurer, thermodynamique  # noqa: F401
from .theorie import TC_ONSAGER  # noqa: F401

__version__ = "2.0.0"
