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
Type stubs for navaltoolbox - High-performance naval architecture library.

This module provides type hints for IDE autocompletion and static type checking.
"""

from typing import List, Tuple

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
]


class Hull:
    """A hull geometry loaded from an STL file.
    
    The Hull class represents a 3D mesh geometry that can be used for
    hydrostatic and stability calculations. It supports loading from STL files,
    transformations, and export operations.
    
    Example:
        >>> hull = Hull("path/to/hull.stl")
        >>> print(hull.num_triangles())
        >>> hull.scale(0.001)  # Convert from mm to meters
    """
    
    def __init__(self, file_path: str) -> None:
        """Load a hull from an STL file.
        
        Args:
            file_path: Path to the STL file (ASCII or binary format).
        
        Raises:
            IOError: If the file cannot be read or parsed.
        """
        ...
    
    def get_bounds(self) -> Tuple[float, float, float, float, float, float]:
        """Returns the bounding box of the hull.
        
        Returns:
            A tuple (xmin, xmax, ymin, ymax, zmin, zmax) in meters.
        """
        ...
    
    def num_triangles(self) -> int:
        """Returns the number of triangles in the mesh."""
        ...
    
    def num_vertices(self) -> int:
        """Returns the number of vertices in the mesh."""
        ...
    
    def transform(
        self,
        translation: Tuple[float, float, float],
        rotation: Tuple[float, float, float],
        pivot: Tuple[float, float, float],
    ) -> None:
        """Applies a transformation to the hull.
        
        Args:
            translation: Translation vector (dx, dy, dz) in meters.
            rotation: Rotation angles (roll, pitch, yaw) in degrees.
            pivot: Pivot point (x, y, z) for rotation.
        """
        ...
    
    def scale(self, factor: float) -> None:
        """Scales the hull uniformly.
        
        Args:
            factor: Scale factor (e.g., 0.001 to convert mm to m).
        """
        ...
    
    def scale_xyz(self, sx: float, sy: float, sz: float) -> None:
        """Scales the hull non-uniformly along each axis.
        
        Args:
            sx: Scale factor along X axis.
            sy: Scale factor along Y axis.
            sz: Scale factor along Z axis.
        """
        ...
    
    def export_stl(self, file_path: str) -> None:
        """Exports the hull to an STL file.
        
        Args:
            file_path: Output path for the STL file.
        
        Raises:
            IOError: If the file cannot be written.
        """
        ...


class Vessel:
    """A vessel containing one or more hulls, tanks, and silhouettes.
    
    The Vessel class is the main container for naval architecture calculations.
    It can hold multiple hulls (for multihull vessels), tanks for loading
    conditions, and silhouettes for wind heeling calculations.
    
    Example:
        >>> hull = Hull("hull.stl")
        >>> vessel = Vessel(hull)
        >>> print(f"LBP: {vessel.lbp:.2f} m")
    """
    
    def __init__(self, hull: Hull) -> None:
        """Create a vessel from a hull.
        
        Args:
            hull: The main hull geometry.
        """
        ...
    
    def get_bounds(self) -> Tuple[float, float, float, float, float, float]:
        """Returns the bounding box of all hulls.
        
        Returns:
            A tuple (xmin, xmax, ymin, ymax, zmin, zmax) in meters.
        """
        ...
    
    @property
    def ap(self) -> float:
        """Returns the Aft Perpendicular position (X coordinate)."""
        ...
    
    @ap.setter
    def ap(self, value: float) -> None:
        """Sets the Aft Perpendicular position (X coordinate)."""
        ...
    
    @property
    def fp(self) -> float:
        """Returns the Forward Perpendicular position (X coordinate)."""
        ...
    
    @fp.setter
    def fp(self, value: float) -> None:
        """Sets the Forward Perpendicular position (X coordinate)."""
        ...
    
    @property
    def lbp(self) -> float:
        """Returns the Length Between Perpendiculars in meters."""
        ...
    
    def num_hulls(self) -> int:
        """Returns the number of hulls."""
        ...
    
    def num_tanks(self) -> int:
        """Returns the number of tanks."""
        ...
    
    def add_tank(self, tank: "Tank") -> None:
        """Add a tank to the vessel.
        
        Args:
            tank: The tank to add.
        """
        ...
    
    def get_total_tanks_mass(self) -> float:
        """Returns the total tanks mass in kg."""
        ...
    
    def get_tanks_center_of_gravity(self) -> Tuple[float, float, float]:
        """Returns the tanks center of gravity [x, y, z] in meters."""
        ...
    
    def add_silhouette(self, silhouette: "Silhouette") -> None:
        """Add a silhouette profile to the vessel.
        
        Args:
            silhouette: The silhouette to add.
        """
        ...
    
    def num_silhouettes(self) -> int:
        """Returns the number of silhouettes."""
        ...
    
    def has_silhouettes(self) -> bool:
        """Returns true if there are any silhouettes."""
        ...
    
    def clear_silhouettes(self) -> None:
        """Removes all silhouettes."""
        ...
    
    def get_total_emerged_area(self, waterline_z: float) -> float:
        """Returns the total emerged area from all silhouettes.
        
        Args:
            waterline_z: Waterline height in meters.
        
        Returns:
            Emerged lateral area in m².
        """
        ...
    
    def get_combined_emerged_centroid(self, waterline_z: float) -> Tuple[float, float]:
        """Returns the combined emerged centroid [x, z].
        
        Args:
            waterline_z: Waterline height in meters.
        
        Returns:
            Centroid coordinates (x, z) in meters.
        """
        ...


class Silhouette:
    """A 2D silhouette profile in the X-Z plane for wind heeling calculations.
    
    Silhouettes represent the lateral projected area of a vessel's superstructure,
    used for calculating wind heeling moments in stability analysis.
    
    Example:
        >>> silhouette = Silhouette("superstructure.dxf")
        >>> print(f"Lateral area: {silhouette.get_area():.2f} m²")
    """
    
    def __init__(self, file_path: str) -> None:
        """Load a silhouette from a DXF file.
        
        Args:
            file_path: Path to the DXF file.
        
        Raises:
            IOError: If the file cannot be read or parsed.
        """
        ...
    
    @staticmethod
    def from_vtk(file_path: str) -> "Silhouette":
        """Load a silhouette from a VTK file (.vtk or .vtp polyline).
        
        Args:
            file_path: Path to the VTK file.
        
        Returns:
            A new Silhouette instance.
        
        Raises:
            IOError: If the file cannot be read or parsed.
        """
        ...
    
    @staticmethod
    def from_points(points: List[Tuple[float, float]], name: str) -> "Silhouette":
        """Create a silhouette from a list of points.
        
        Args:
            points: List of (x, z) coordinates defining the contour.
            name: Name identifier for the silhouette.
        
        Returns:
            A new Silhouette instance.
        """
        ...
    
    @property
    def name(self) -> str:
        """Returns the silhouette name."""
        ...
    
    def num_points(self) -> int:
        """Returns the number of points in the contour."""
        ...
    
    def is_closed(self) -> bool:
        """Returns true if the contour is closed."""
        ...
    
    def get_points(self) -> List[Tuple[float, float, float]]:
        """Returns the points as a list of tuples [(x, y, z), ...]."""
        ...
    
    def get_area(self) -> float:
        """Returns the total lateral area in m²."""
        ...
    
    def get_centroid(self) -> Tuple[float, float]:
        """Returns the centroid [x, z] in meters."""
        ...
    
    def get_bounds(self) -> Tuple[float, float, float, float]:
        """Returns the bounding box (x_min, x_max, z_min, z_max)."""
        ...
    
    def get_emerged_area(self, waterline_z: float) -> float:
        """Returns the emerged area above waterline in m².
        
        Args:
            waterline_z: Waterline height in meters.
        """
        ...
    
    def get_emerged_centroid(self, waterline_z: float) -> Tuple[float, float]:
        """Returns the centroid of emerged area [x, z].
        
        Args:
            waterline_z: Waterline height in meters.
        """
        ...


class OpeningType:
    """Type of opening that can cause downflooding.
    
    Use the static methods to create specific opening types:
    - OpeningType.vent()
    - OpeningType.air_pipe()
    - OpeningType.hatch()
    - OpeningType.door()
    - OpeningType.window()
    - OpeningType.other("custom_name")
    """
    
    @staticmethod
    def vent() -> "OpeningType":
        """Creates a vent opening type."""
        ...
    
    @staticmethod
    def air_pipe() -> "OpeningType":
        """Creates an air pipe opening type."""
        ...
    
    @staticmethod
    def hatch() -> "OpeningType":
        """Creates a hatch opening type."""
        ...
    
    @staticmethod
    def door() -> "OpeningType":
        """Creates a door opening type."""
        ...
    
    @staticmethod
    def window() -> "OpeningType":
        """Creates a window opening type."""
        ...
    
    @staticmethod
    def other(name: str) -> "OpeningType":
        """Creates a custom opening type.
        
        Args:
            name: Custom name for the opening type.
        """
        ...


class DownfloodingOpening:
    """A downflooding opening point or contour.
    
    Downflooding openings represent locations where water can enter the vessel
    when heeled or trimmed. They are used in IMO intact stability calculations
    to determine flooding angles.
    
    Example:
        >>> opening = DownfloodingOpening.from_point(
        ...     "engine_vent", (50.0, 2.5, 8.0), OpeningType.vent()
        ... )
        >>> if opening.is_submerged(heel=30.0, trim=0.0, pivot=(50, 0, 5), waterline_z=3.0):
        ...     print("Opening flooded!")
    """
    
    @staticmethod
    def from_point(
        name: str,
        position: Tuple[float, float, float],
        opening_type: OpeningType,
    ) -> "DownfloodingOpening":
        """Create a downflooding opening from a single point.
        
        Args:
            name: Identifier for the opening.
            position: Position (x, y, z) in meters.
            opening_type: Type of opening.
        
        Returns:
            A new DownfloodingOpening instance.
        """
        ...
    
    @staticmethod
    def from_contour(
        name: str,
        points: List[Tuple[float, float, float]],
        opening_type: OpeningType,
    ) -> "DownfloodingOpening":
        """Create a downflooding opening from a contour (polyline).
        
        Args:
            name: Identifier for the opening.
            points: List of (x, y, z) coordinates defining the contour.
            opening_type: Type of opening.
        
        Returns:
            A new DownfloodingOpening instance.
        """
        ...
    
    @property
    def name(self) -> str:
        """Returns the opening name."""
        ...
    
    @property
    def is_active(self) -> bool:
        """Check if opening is active (considered in calculations)."""
        ...
    
    def set_active(self, active: bool) -> None:
        """Set opening active state.
        
        Args:
            active: True to include in calculations, False to ignore.
        """
        ...
    
    def num_points(self) -> int:
        """Get number of points defining the opening."""
        ...
    
    def get_points(self) -> List[Tuple[float, float, float]]:
        """Get all points as [(x, y, z), ...]."""
        ...
    
    def is_submerged(
        self,
        heel: float,
        trim: float,
        pivot: Tuple[float, float, float],
        waterline_z: float,
    ) -> bool:
        """Check if the opening is submerged at given heel/trim/draft.
        
        Args:
            heel: Heel angle in degrees (positive = starboard down).
            trim: Trim angle in degrees (positive = stern down).
            pivot: Rotation pivot point (x, y, z) in meters.
            waterline_z: Waterline height in meters.
        
        Returns:
            True if any point of the opening is below the waterline.
        """
        ...


class HydrostaticState:
    """Result of hydrostatic calculations.
    
    Contains the hydrostatic properties at a specific floating condition.
    
    Attributes:
        draft: Draft at midship in meters.
        trim: Trim angle in degrees.
        heel: Heel angle in degrees.
        volume: Submerged volume in m³.
        displacement: Displacement mass in kg.
        cob: Center of buoyancy as tuple (lcb, tcb, vcb).
        cog: Center of gravity as tuple (lcg, tcg, vcg) if specified, None otherwise.
        lcb: Longitudinal center of buoyancy (X) in meters.
        tcb: Transverse center of buoyancy (Y) in meters.
        vcb: Vertical center of buoyancy (Z) in meters.
        lcg: Longitudinal center of gravity (X) in meters, or None.
        tcg: Transverse center of gravity (Y) in meters, or None.
        vcg: Vertical center of gravity (Z) in meters, or None.
        waterplane_area: Waterplane area in m².
        lcf: Longitudinal center of floatation (X) in meters.
        bmt: Transverse metacentric radius in meters.
        bml: Longitudinal metacentric radius in meters.
        gmt: Transverse metacentric height with FSC in meters, or None.
        gml: Longitudinal metacentric height with FSC in meters, or None.
        gmt_dry: Transverse metacentric height without FSC in meters, or None.
        gml_dry: Longitudinal metacentric height without FSC in meters, or None.
        lwl: Waterline length in meters.
        bwl: Waterline breadth in meters.
        los: Length overall submerged in meters.
        wetted_surface_area: Wetted surface area in m².
        midship_area: Midship section area in m².
        cm: Midship coefficient.
        cb: Block coefficient.
        cp: Prismatic coefficient.
        free_surface_correction_t: Transverse FSC in meters.
        free_surface_correction_l: Longitudinal FSC in meters.
        stiffness_matrix: 6x6 hydrostatic stiffness matrix (flattened).
    """
    
    draft: float
    trim: float
    heel: float
    volume: float
    displacement: float
    
    @property
    def cob(self) -> Tuple[float, float, float]:
        """Center of buoyancy (lcb, tcb, vcb) in meters."""
        ...
    
    @property
    def cog(self) -> Tuple[float, float, float] | None:
        """Center of gravity (lcg, tcg, vcg) if specified, None otherwise."""
        ...
    
    @property
    def lcb(self) -> float:
        """Longitudinal center of buoyancy (X) in meters."""
        ...
    
    @property
    def tcb(self) -> float:
        """Transverse center of buoyancy (Y) in meters."""
        ...
    
    @property
    def vcb(self) -> float:
        """Vertical center of buoyancy (Z) in meters."""
        ...
    
    @property
    def lcg(self) -> float | None:
        """Longitudinal center of gravity (X) in meters, or None."""
        ...
    
    @property
    def tcg(self) -> float | None:
        """Transverse center of gravity (Y) in meters, or None."""
        ...
    
    @property
    def vcg(self) -> float | None:
        """Vertical center of gravity (Z) in meters, or None."""
        ...
    
    waterplane_area: float
    lcf: float
    bmt: float
    bml: float
    
    @property
    def gmt(self) -> float | None:
        """Transverse metacentric height with FSC in meters, or None if VCG not specified."""
        ...
    
    @property
    def gml(self) -> float | None:
        """Longitudinal metacentric height with FSC in meters, or None if VCG not specified."""
        ...
    
    @property
    def gmt_dry(self) -> float | None:
        """Transverse metacentric height without FSC in meters, or None."""
        ...
    
    @property
    def gml_dry(self) -> float | None:
        """Longitudinal metacentric height without FSC in meters, or None."""
        ...
    
    lwl: float
    bwl: float
    los: float
    wetted_surface_area: float
    midship_area: float
    cm: float
    cb: float
    cp: float
    free_surface_correction_t: float
    free_surface_correction_l: float
    stiffness_matrix: List[float]


class HydrostaticsCalculator:
    """Calculator for hydrostatic properties.
    
    Performs hydrostatic calculations on a vessel, including volume,
    center of buoyancy, and draft finding.
    
    Example:
        >>> hull = Hull("hull.stl")
        >>> vessel = Vessel(hull)
        >>> calc = HydrostaticsCalculator(vessel, water_density=1025.0)
        >>> state = calc.from_draft(draft=5.0)
        >>> print(f"Displacement: {state.displacement:.0f} kg")
    """
    
    def __init__(self, vessel: Vessel, water_density: float = 1025.0) -> None:
        """Create a hydrostatics calculator for a vessel.
        
        Args:
            vessel: The vessel to analyze.
            water_density: Water density in kg/m³ (default: seawater 1025).
        """
        ...
    
    def from_draft(
        self,
        draft: float,
        trim: float = 0.0,
        heel: float = 0.0,
        vcg: float | None = None,
    ) -> HydrostaticState:
        """Calculate hydrostatics at a given draft, trim, and heel.
        
        Args:
            draft: Draft at midship in meters.
            trim: Trim angle in degrees (default: 0).
            heel: Heel angle in degrees (default: 0).
            vcg: Optional vertical center of gravity in meters for GMT/GML calculation.
        
        Returns:
            HydrostaticState with calculated properties.
        
        Raises:
            ValueError: If no submerged volume at this draft.
        """
        ...
    
    def from_displacement(
        self,
        displacement_mass: float,
        vcg: float | None = None,
        cog: Tuple[float, float, float] | None = None,
        trim: float | None = None,
        heel: float | None = None,
    ) -> HydrostaticState:
        """Calculate hydrostatics for a given displacement with optional constraints.
        
        Args:
            displacement_mass: Target displacement in kg.
            vcg: Optional vertical center of gravity (m) for GM calculations.
            cog: Optional (lcg, tcg, vcg) tuple in meters for full COG specification.
                 (overrides vcg if both are provided)
            trim: Optional trim angle in degrees.
            heel: Optional heel angle in degrees.
        
        Returns:
            Complete HydrostaticState.
        
        Raises:
            ValueError: If constraints are invalid or unsatisfiable.
        
        Examples:
            >>> # Basic: find draft for displacement
            >>> state = calc.from_displacement(8635000.0)
            
            >>> # With VCG only: compute GMT/GML
            >>> state = calc.from_displacement(8635000.0, vcg=7.555)
            
            >>> # With full COG: for trim optimization
            >>> state = calc.from_displacement(8635000.0, cog=(71.67, 0.0, 7.555))
        """
        ...
    
    @property
    def water_density(self) -> float:
        """Returns the water density in kg/m³."""
        ...


class StabilityPoint:
    """A point on a stability curve.
    
    Represents the stability properties at a specific heel angle.
    
    Attributes:
        heel: Heel angle in degrees.
        draft: Draft at this heel angle in meters.
        trim: Trim angle at this heel in degrees.
        gz: Righting arm (GZ) in meters.
        is_flooding: True if any downflooding opening is submerged.
        flooded_openings: List of names of submerged openings.
    """
    
    heel: float
    draft: float
    trim: float
    gz: float
    is_flooding: bool
    flooded_openings: List[str]


class StabilityCurve:
    """A complete GZ stability curve.
    
    Contains the full righting arm curve for a specific loading condition.
    
    Example:
        >>> curve = calc.calculate_gz_curve(
        ...     displacement_mass=10000,
        ...     cog=(50.0, 0.0, 5.0),
        ...     heels=[0, 10, 20, 30, 40, 50, 60]
        ... )
        >>> for heel, gz in zip(curve.heels(), curve.values()):
        ...     print(f"Heel: {heel}°, GZ: {gz:.3f} m")
    """
    
    @property
    def displacement(self) -> float:
        """Returns the displacement in kg."""
        ...
    
    def heels(self) -> List[float]:
        """Returns the heel angles in degrees."""
        ...
    
    def values(self) -> List[float]:
        """Returns the GZ values in meters."""
        ...
    
    def points(self) -> List[Tuple[float, float, float, float]]:
        """Returns the points as a list of tuples (heel, draft, trim, gz)."""
        ...


class StabilityCalculator:
    """Calculator for stability curves (GZ).
    
    Performs intact stability calculations, generating GZ curves for
    specified loading conditions.
    
    Example:
        >>> hull = Hull("hull.stl")
        >>> vessel = Vessel(hull)
        >>> calc = StabilityCalculator(vessel, water_density=1025.0)
        >>> curve = calc.calculate_gz_curve(
        ...     displacement_mass=50000,
        ...     cog=(45.0, 0.0, 6.5),
        ...     heels=list(range(0, 91, 5))
        ... )
    """
    
    def __init__(self, vessel: Vessel, water_density: float = 1025.0) -> None:
        """Create a stability calculator for a vessel.
        
        Args:
            vessel: The vessel to analyze.
            water_density: Water density in kg/m³ (default: seawater 1025).
        """
        ...
    
    def calculate_gz_curve(
        self,
        displacement_mass: float,
        cog: Tuple[float, float, float],
        heels: List[float],
    ) -> StabilityCurve:
        """Calculate the GZ curve for a given loading condition.
        
        Args:
            displacement_mass: Displacement in kg.
            cog: Center of gravity (x, y, z) in meters.
            heels: List of heel angles in degrees to calculate.
        
        Returns:
            StabilityCurve with GZ values at each heel angle.
        """
        ...


class Tank:
    """A tank with fluid management capabilities.
    
    Represents a fluid tank on a vessel, with support for fill level
    management and free surface effect calculations.
    
    Example:
        >>> tank = Tank.from_box(
        ...     name="FO_1P",
        ...     x_min=20.0, x_max=25.0,
        ...     y_min=-5.0, y_max=0.0,
        ...     z_min=0.0, z_max=3.0,
        ...     fluid_density=850.0  # Fuel oil
        ... )
        >>> tank.fill_percent = 80.0
        >>> print(f"Mass: {tank.fluid_mass:.0f} kg")
    """
    
    @staticmethod
    def from_box(
        name: str,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
        z_min: float,
        z_max: float,
        fluid_density: float,
    ) -> "Tank":
        """Create a box-shaped tank.
        
        Args:
            name: Tank identifier (e.g., "FO_1P", "FW_2S").
            x_min: Minimum X coordinate in meters.
            x_max: Maximum X coordinate in meters.
            y_min: Minimum Y coordinate in meters (negative = port).
            y_max: Maximum Y coordinate in meters (positive = starboard).
            z_min: Minimum Z coordinate in meters (bottom).
            z_max: Maximum Z coordinate in meters (top).
            fluid_density: Fluid density in kg/m³.
        
        Returns:
            A new Tank instance.
        """
        ...
    
    @property
    def name(self) -> str:
        """Returns the tank name."""
        ...
    
    @property
    def total_volume(self) -> float:
        """Returns the total volume in m³."""
        ...
    
    @property
    def fill_level(self) -> float:
        """Returns the fill level as a fraction (0-1)."""
        ...
    
    @fill_level.setter
    def fill_level(self, level: float) -> None:
        """Sets the fill level as a fraction (0-1)."""
        ...
    
    @property
    def fill_percent(self) -> float:
        """Returns the fill level as a percentage (0-100)."""
        ...
    
    @fill_percent.setter
    def fill_percent(self, percent: float) -> None:
        """Sets the fill level as a percentage (0-100)."""
        ...
    
    @property
    def fill_volume(self) -> float:
        """Returns the filled volume in m³."""
        ...
    
    @property
    def fluid_mass(self) -> float:
        """Returns the fluid mass in kg."""
        ...
    
    @property
    def center_of_gravity(self) -> Tuple[float, float, float]:
        """Returns the center of gravity [x, y, z] in meters."""
        ...
    
    @property
    def free_surface_moment_t(self) -> float:
        """Returns the transverse free surface moment in m⁴."""
        ...
    
    @property
    def free_surface_moment_l(self) -> float:
        """Returns the longitudinal free surface moment in m⁴."""
        ...
