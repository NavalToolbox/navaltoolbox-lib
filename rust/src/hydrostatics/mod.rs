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

//! Hydrostatics module.
//!
//! Provides hydrostatic calculations for vessels.

mod calculator;
mod dataclasses;
pub use crate::mesh::{calculate_waterplane_properties, WaterplaneProperties};
pub use calculator::HydrostaticsCalculator;
pub use dataclasses::{HydrostaticState, TankOptions};

pub(crate) use calculator::{
    build_contact_face_refs, calculate_mesh_area, compute_contact_area_from_precomputed,
    detect_contact_area,
};
