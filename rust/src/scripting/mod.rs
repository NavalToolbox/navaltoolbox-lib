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

//! Scripting module for stability criteria verification.
//!
//! Provides Rhai scripting support for implementing portable stability criteria
//! that work from both Python bindings and pure Rust applications.

mod context;
mod engine;
mod result;

pub use context::CriteriaContext;
pub use engine::{ScriptEngine, ScriptError};
pub use result::{CriteriaResult, CriteriaStatus, CriterionResult, PlotData, PlotElement};
