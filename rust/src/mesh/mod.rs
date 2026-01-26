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

//! Mesh utilities module.
//!
//! Provides mesh loading, transformation, and clipping functionality.

mod clipper;
mod loader;
mod transform;

pub use clipper::{clip_at_waterline, clip_by_axis_aligned_plane, Axis};
pub use loader::{load_stl, load_vtk};
pub use transform::{get_bounds, transform_mesh, transform_point};
