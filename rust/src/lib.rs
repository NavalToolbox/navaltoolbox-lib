// Copyright (C) 2026 Antoine ANCEAU
//
// This file is part of navaltoolbox.
//
// navaltoolbox is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

//! # NavalToolbox
//!
//! High-performance naval architecture library for hydrostatics, stability, and tank calculations.
//!
//! ## Features
//! - Hull geometry loading (STL, VTK)
//! - Multi-hull vessel support
//! - Hydrostatic calculations (volume, CoB, waterplane properties, metacentric heights)
//! - Stability analysis (KN curves, GZ curves)
//! - Tank management with free surface effects
//!
//! ## Usage
//! ```rust,ignore
//! use navaltoolbox::{Hull, Vessel, HydrostaticsCalculator};
//!
//! let hull = Hull::from_stl("ship.stl")?;
//! let vessel = Vessel::new(hull);
//! let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
//! let state = calc.from_draft(5.0, 0.0, 0.0)?;
//! println!("Volume: {} m³", state.volume);
//! ```
//!
//! ## Disclaimer
//!
//! NavalToolbox has been developed with care to ensure that all models and methods are correct.
//! However, **results must not be considered as a guarantee of performance**. The author cannot
//! be held responsible for any inaccuracies in the calculations or for any consequences arising
//! from the use of this software. Users are advised to independently verify critical calculations
//! and to use this software as a tool to support, not replace, professional engineering judgment.

pub mod downflooding;
pub mod hull;
pub mod hydrostatics;
pub mod mesh;
pub mod scripting;
pub mod silhouette;
pub mod stability;
pub mod tanks;
pub mod vessel;

// Re-export main types
pub use downflooding::{DownfloodingOpening, OpeningGeometry, OpeningType};
pub use hull::Hull;
pub use hydrostatics::{HydrostaticState, HydrostaticsCalculator};
pub use scripting::{
    CriteriaContext, CriteriaResult, CriteriaStatus, CriterionResult, PlotData, PlotElement,
    ScriptEngine, ScriptError,
};
pub use silhouette::Silhouette;
pub use stability::{
    CompleteStabilityResult, StabilityCalculator, StabilityCurve, StabilityPoint, WindHeelingData,
};
pub use tanks::{Tank, TankState};
pub use vessel::Vessel;

// ============================================================================
// Python Bindings
// ============================================================================

#[cfg(feature = "python")]
mod python;

#[cfg(feature = "python")]
pub use python::*;
