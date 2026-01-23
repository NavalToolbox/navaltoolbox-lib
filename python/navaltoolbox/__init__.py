# Copyright (C) 2026 Antoine ANCEAU
#
# This file is part of navaltoolbox.
#
# navaltoolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
NavalToolbox - High-performance naval architecture library.

This library provides tools for hydrostatics, stability analysis,
and tank management in naval architecture applications.

Classes:
    Hull: A hull geometry loaded from an STL file.
    Vessel: A vessel containing hulls, tanks, and silhouettes.
    Silhouette: A 2D profile for wind heeling calculations.
    Tank: A tank with fluid management capabilities.
    HydrostaticsCalculator: Calculator for hydrostatic properties.
    StabilityCalculator: Calculator for stability curves (GZ).
    DownfloodingOpening: Openings for downflooding analysis.
    OpeningType: Types of downflooding openings.

Example:
    >>> from navaltoolbox import Hull, Vessel, HydrostaticsCalculator
    >>> hull = Hull("ship.stl")
    >>> vessel = Vessel(hull)
    >>> calc = HydrostaticsCalculator(vessel)
    >>> state = calc.from_draft(5.0)
    >>> print(f"Displacement: {state.displacement:.0f} kg")
"""

# Re-export all types from the native module
from .navaltoolbox import (
    Hull,
    Vessel,
    Silhouette,
    OpeningType,
    DownfloodingOpening,
    HydrostaticState,
    HydrostaticsCalculator,
    StabilityPoint,
    StabilityCurve,
    StabilityCalculator,
    Tank,
    CompleteStabilityResult,
    CriterionResult,
    CriteriaResult,
    CriteriaContext,
    ScriptEngine,
)

__all__ = [
    "Hull",
    "Vessel",
    "Silhouette",
    "OpeningType",
    "DownfloodingOpening",
    "HydrostaticState",
    "HydrostaticsCalculator",
    "StabilityPoint",
    "StabilityCurve",
    "StabilityCalculator",
    "Tank",
    "CompleteStabilityResult",
    "CriterionResult",
    "CriteriaResult",
    "CriteriaContext",
    "ScriptEngine",
]
