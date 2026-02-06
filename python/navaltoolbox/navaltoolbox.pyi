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
    "Appendage",
    "DeckEdge",
    "DeckEdgeSide",
    "OpeningType",
    "DownfloodingOpening",
    "HydrostaticState",
    "HydrostaticsCalculator",
    "StabilityPoint",
    "StabilityCurve",
    "StabilityCalculator",
    "Tank",
    "WindHeelingData",
    "CompleteStabilityResult",
    "CriterionResult",
    "CriteriaResult",
    "CriteriaContext",
    "ScriptEngine",
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
    
    @staticmethod
    def from_box(
        length: float,
        breadth: float,
        depth: float,
    ) -> "Hull":
        """Create a box hull.
        
        Args:
            length: Length of the box in meters.
            breadth: Breadth of the box in meters.
            depth: Depth of the box in meters.
        
        Returns:
            A new Hull instance.
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
    
    def simplify(self, target_count: int) -> None:
        """Simplifies the hull mesh to a target number of triangles.
        
        Args:
            target_count: Target number of triangles for the simplified mesh.
        """
        ...
    
    def to_simplified(self, target_count: int) -> "Hull":
        """Returns a simplified copy of the hull.
        
        Args:
            target_count: Target number of triangles for the simplified mesh.
        
        Returns:
            A new simplified Hull instance.
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
    
    def get_vertices(self) -> List[Tuple[float, float, float]]:
        """Returns vertices as list of tuples (x, y, z)."""
        ...
    
    def get_faces(self) -> List[Tuple[int, int, int]]:
        """Returns faces as list of tuples (i, j, k)."""
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
    
    def add_opening(self, opening: "DownfloodingOpening") -> None:
        """Add a downflooding opening to the vessel.
        
        Args:
            opening: The opening to add.
        """
        ...

    def num_openings(self) -> int:
        """Returns the number of downflooding openings."""
        ...
    
    def clear_openings(self) -> None:
        """Removes all downflooding openings."""
        ...

    def get_hulls(self) -> List["Hull"]:
        """Get all hulls.

        Returns:
            List of Hull objects.
        """
        ...

    def get_tanks(self) -> List["Tank"]:
        """Get all tanks.

        Returns:
            List of Tank objects.
        """
        ...

    def get_silhouettes(self) -> List["Silhouette"]:
        """Get all silhouettes.

        Returns:
            List of Silhouette objects.
        """
        ...

    def get_openings(self) -> List["DownfloodingOpening"]:
        """Get all downflooding openings.

        Returns:
            List of DownfloodingOpening objects.
        """
        ...
    
    # Appendage methods
    
    def add_appendage(self, appendage: "Appendage") -> None:
        """Add an appendage to the vessel."""
        ...
        
    def num_appendages(self) -> int:
        """Returns the number of appendages."""
        ...
        
    def clear_appendages(self) -> None:
        """Removes all appendages."""
        ...
        
    def get_appendages(self) -> List["Appendage"]:
        """Get all appendages."""
        ...
        
    def get_total_appendage_volume(self) -> float:
        """Returns the total appendage volume in m³."""
        ...
        
    def get_total_appendage_wetted_surface(self) -> float:
        """Returns the total appendage wetted surface in m²."""
        ...
        
    # Deck Edge methods
    
    def add_deck_edge(self, deck_edge: "DeckEdge") -> None:
        """Add a deck edge to the vessel."""
        ...
        
    def num_deck_edges(self) -> int:
        """Returns the number of deck edges."""
        ...
        
    def has_deck_edges(self) -> bool:
        """Returns true if any deck edges are defined."""
        ...
        
    def clear_deck_edges(self) -> None:
        """Removes all deck edges."""
        ...
        
    def get_deck_edges(self) -> List["DeckEdge"]:
        """Get all deck edges."""
        ...
        
    def get_min_freeboard(self, heel: float, trim: float, waterline_z: float) -> float | None:
        """Calculate minimum freeboard across all deck edges."""
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
        """Load a silhouette from a file (DXF, VTK, or VTP).
        
        Args:
            file_path: Path to the geometry file.
        
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



class Appendage:
    """An appendage (additional volume element) attached to the vessel.
    
    Appendages represent volume contributions from items like keels, rudders,
    bulbous bows, etc. that are not part of the main hull geometry.
    """
    
    @staticmethod
    def from_point(name: str, center: Tuple[float, float, float], volume: float) -> "Appendage":
        """Create an appendage from a point (fixed volume at position)."""
        ...
        
    @staticmethod
    def from_file(name: str, file_path: str) -> "Appendage":
        """Create an appendage from an STL or VTK file."""
        ...
        
    @staticmethod
    def from_box(name: str, xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float) -> "Appendage":
        """Create an appendage from a box (parallelepiped)."""
        ...
        
    @staticmethod
    def from_cube(name: str, center: Tuple[float, float, float], volume: float) -> "Appendage":
        """Create an appendage from a cube (center and volume)."""
        ...
        
    @staticmethod
    def from_sphere(name: str, center: Tuple[float, float, float], volume: float) -> "Appendage":
        """Create an appendage from a sphere (center and volume)."""
        ...
    
    @property
    def name(self) -> str:
        """The appendage name."""
        ...
    
    @name.setter
    def name(self, value: str) -> None:
        ...
        
    @property
    def volume(self) -> float:
        """Volume in m³."""
        ...
        
    @property
    def center(self) -> Tuple[float, float, float]:
        """Center of volume (x, y, z) in meters."""
        ...
        
    @property
    def wetted_surface(self) -> float | None:
        """Wetted surface area in m², or None if not set."""
        ...
        
    @wetted_surface.setter
    def wetted_surface(self, value: float | None) -> None:
        ...
        
    def geometry_type(self) -> str:
        """Returns the geometry type (Point, Mesh, Box, etc.)."""
        ...
    
    def get_mesh_data(self) -> Tuple[List[Tuple[float, float, float]], List[Tuple[int, int, int]]] | None:
        """Returns mesh data (vertices, faces) if geometry is a mesh.
        
        Returns:
            Tuple of (vertices, faces), or None if not a mesh.
        """
        ...
    
    @property
    def bounds(self) -> Tuple[float, float, float, float, float, float] | None:
        """Returns bounds (xmin, xmax, ymin, ymax, zmin, zmax)."""
        ...


class DeckEdgeSide:
    """Side of the deck edge (Port, Starboard, or Both)."""
    
    @staticmethod
    def port() -> "DeckEdgeSide": ...
    
    @staticmethod
    def starboard() -> "DeckEdgeSide": ...
    
    @staticmethod
    def both() -> "DeckEdgeSide": ...


class DeckEdge:
    """A deck edge contour (livet) for freeboard calculation."""
    
    @staticmethod
    def from_points(name: str, points: List[Tuple[float, float, float]], side: DeckEdgeSide) -> "DeckEdge":
        """Create a deck edge from a list of 3D points."""
        ...
        
    @staticmethod
    def from_file(name: str, file_path: str) -> "DeckEdge":
        """Load a deck edge from a DXF or VTK file."""
        ...
    
    @property
    def name(self) -> str:
        """The deck edge name."""
        ...
    
    @name.setter
    def name(self, value: str) -> None:
        ...
        
    def num_points(self) -> int:
        """Returns the number of points."""
        ...
        
    def get_points(self) -> List[Tuple[float, float, float]]:
        """Returns points as list of (x, y, z) tuples."""
        ...
        
    def get_side(self) -> str:
        """Returns the side as a string."""
        ...
        
    def get_freeboard(self, heel: float, trim: float, pivot: Tuple[float, float, float], waterline_z: float) -> float:
        """Calculate freeboard at given conditions."""
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
    
    @staticmethod
    def from_file(file_path: str, default_type: "OpeningType", name: str | None = None) -> List["DownfloodingOpening"]:
        """Load openings from a file (DXF or VTK).
        
        If name is provided:
        - For a single opening, it sets the name.
        - For multiple openings, it sets names as "{name}_{i+1}".
        
        Args:
            file_path: Path to the file.
            default_type: Default OpeningType for loaded openings.
            name: Optional base name for loaded openings.
            
        Returns:
            List of DownfloodingOpening instances.
        
        Raises:
            IOError: If the file cannot be read or parsed.
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
        draft_ap: Draft at Aft Perpendicular in meters.
        draft_fp: Draft at Forward Perpendicular in meters.
        draft_mp: Draft at Midship Perpendicular in meters.
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
        sectional_areas: List of (x, area) pairs representing the sectional area curve.
        freeboard: Minimum freeboard in meters, or None.
    """
    
    draft: float
    trim: float
    heel: float
    draft_ap: float
    draft_fp: float
    draft_mp: float
    volume: float
    displacement: float
    
    sectional_areas: List[Tuple[float, float]]
    freeboard: float | None
    
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
        num_stations: int | None = None,
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
    
    def from_drafts(
        self,
        draft_ap: float,
        draft_fp: float,
        heel: float = 0.0,
        vcg: float | None = None,
    ) -> HydrostaticState:
        """Calculate hydrostatics from drafts at Aft and Forward Perpendiculars.
        
        Args:
            draft_ap: Draft at Aft Perpendicular in meters.
            draft_fp: Draft at Forward Perpendicular in meters.
            heel: Heel angle in degrees (default: 0).
            vcg: Optional vertical center of gravity in meters for GMT/GML calculation.
        
        Returns:
            HydrostaticState with calculated properties.
        
        Raises:
            ValueError: If no submerged volume at these drafts.
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
        >>> curve = calc.gz_curve(
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
        
    def get_stability_points(self) -> List["StabilityPoint"]:
        """Returns the points as a list of StabilityPoint objects.
        
        This allows access to detailed information like is_flooding.
        """
        ...



class WindHeelingData:
    """Wind heeling data from silhouette calculation."""
    
    @property
    def emerged_area(self) -> float: ...
    @property
    def emerged_centroid(self) -> Tuple[float, float]: ...
    @property
    def wind_lever_arm(self) -> float: ...
    @property
    def waterline_z(self) -> float: ...


class CompleteStabilityResult:
    """Result of complete stability analysis."""
    
    @property
    def hydrostatics(self) -> HydrostaticState:
        """Hydrostatic state at equilibrium."""
        ...
    
    @property
    def gz_curve(self) -> StabilityCurve:
        """GZ stability curve."""
        ...
    
    @property
    def wind_data(self) -> WindHeelingData | None:
        """Wind heeling data (if silhouettes are defined)."""
        ...
    
    @property
    def displacement(self) -> float:
        """Displacement mass in kg."""
        ...
    
    @property
    def cog(self) -> Tuple[float, float, float]:
        """Center of gravity (LCG, TCG, VCG)."""
        ...
    
    @property
    def gm0(self) -> float | None:
        """Initial GM (fluid)."""
        ...
    
    @property
    def gm0_dry(self) -> float | None:
        """Initial GM (dry/solid)."""
        ...
    
    @property
    def max_gz(self) -> float | None:
        """Maximum GZ value."""
        ...
    
    @property
    def heel_at_max_gz(self) -> float | None:
        """Heel angle (degrees) at maximum GZ."""
        ...
    
    def has_wind_data(self) -> bool:
        """Returns true if wind data is available."""
        ...


class StabilityCalculator:
    """Calculator for stability curves (GZ).
    
    Performs intact stability calculations, generating GZ curves for
    specified loading conditions.
    
    Example:
        >>> hull = Hull("hull.stl")
        >>> vessel = Vessel(hull)
        >>> calc = StabilityCalculator(vessel, water_density=1025.0)
        >>> curve = calc.gz_curve(
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
    
    def gz_curve(
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
    
    def kn_curve(
        self,
        displacements: List[float],
        heels: List[float],
        lcg: float = 0.0,
        tcg: float = 0.0,
    ) -> List[StabilityCurve]:
        """Calculate KN curves (Righting Lever from Keel) for multiple displacements.
        
        This calculates stability curves assuming VCG = 0.
        Returns one curve per displacement.
        
        Args:
            displacements: List of displacements in kg.
            heels: List of heel angles in degrees.
            lcg: Longitudinal Center of Gravity in meters (default 0.0).
            tcg: Transverse Center of Gravity in meters (default 0.0).
        
        Returns:
            List[StabilityCurve]: One curve per displacement.
        """
        ...
    
    def complete_stability(
        self,
        displacement_mass: float,
        cog: Tuple[float, float, float],
        heels: List[float],
    ) -> CompleteStabilityResult:
        """Calculate complete stability analysis for a loading condition.
        
        Combines hydrostatic calculations, GZ curve, and wind heeling data
        (if silhouettes are available) for a single loading condition.
        
        Args:
            displacement_mass: Target displacement in kg.
            cog: Center of gravity (lcg, tcg, vcg) tuple.
            heels: List of heel angles for GZ curve in degrees.
        
        Returns:
            CompleteStabilityResult with hydrostatics, GZ curve, and wind data.
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
    
    def __init__(
        self,
        file_path: str,
        fluid_density: float = 1025.0,
        name: str | None = None,
    ) -> None:
        """Create a Tank from a file (STL or VTK).
        
        Args:
            file_path: Path to the geometry file.
            fluid_density: Fluid density in kg/m³.
            name: Optional name for the tank.
        """
        ...

    @staticmethod
    def from_box_hull_intersection(
        hull: Hull,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
        z_min: float,
        z_max: float,
        fluid_density: float = 1025.0,
        name: str = "HullTank",
    ) -> "Tank":
        """Create a Tank as the intersection of a box with a hull geometry.
        
        Args:
            hull: Hull object to intersect with.
            x_min: Minimum X coordinate in meters.
            x_max: Maximum X coordinate in meters.
            y_min: Minimum Y coordinate in meters.
            y_max: Maximum Y coordinate in meters.
            z_min: Minimum Z coordinate in meters.
            z_max: Maximum Z coordinate in meters.
            fluid_density: Fluid density in kg/m³.
            name: Tank name.
        
        Returns:
            A new Tank instance.
        """
        ...

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
            y_min: Minimum Y coordinate in meters (negative = starboard).
            y_max: Maximum Y coordinate in meters (positive = port).
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
    
    def get_vertices(self) -> List[Tuple[float, float, float]]:
        """Returns tank container vertices [(x,y,z)]."""
        ...
    
    def get_faces(self) -> List[Tuple[int, int, int]]:
        """Returns tank container faces [(i,j,k)]."""
        ...
    
    def get_fluid_vertices(self, heel: float = 0.0, trim: float = 0.0) -> List[Tuple[float, float, float]]:
        """Returns fluid mesh vertices [(x,y,z)] or empty list."""
        ...
    
    def get_fluid_faces(self, heel: float = 0.0, trim: float = 0.0) -> List[Tuple[int, int, int]]:
        """Returns fluid mesh faces [(i,j,k)] or empty list."""
        ...


class CriterionResult:
    """Result of a single criterion check."""
    
    @property
    def name(self) -> str: ...
    @property
    def description(self) -> str: ...
    @property
    def required_value(self) -> float: ...
    @property
    def actual_value(self) -> float: ...
    @property
    def unit(self) -> str: ...
    @property
    def status(self) -> str:
        """Pass/Fail status ('PASS', 'FAIL', 'WARN', 'N/A')."""
        ...
    @property
    def margin(self) -> float:
        """Margin between actual and required value."""
        ...
    @property
    def notes(self) -> str | None: ...
    @property
    def plot_id(self) -> str | None: ...


class CriteriaResult:
    """Result of a criteria verification script."""
    
    @property
    def regulation_name(self) -> str: ...
    @property
    def regulation_reference(self) -> str: ...
    @property
    def vessel_name(self) -> str: ...
    @property
    def loading_condition(self) -> str: ...
    @property
    def displacement(self) -> float: ...
    @property
    def overall_pass(self) -> bool: ...
    @property
    def pass_count(self) -> int: ...
    @property
    def fail_count(self) -> int: ...
    @property
    def notes(self) -> str: ...
    @property
    def criteria(self) -> List[CriterionResult]: ...
    @property
    def plots(self) -> List[str]:
        """List of JSON strings representing plot data."""
        ...


class CriteriaContext:
    """Context for Rhai scripts, wrapping stability results.
    
    This context provides access to stability data within Rhai scripts.
    The following methods are available when writing Rhai scripts:
    
    **GZ Curve Methods:**
    - `get_heels()` -> list of heel angles
    - `get_gz_values()` -> list of GZ values
    - `area_under_curve(from, to)` -> area in m·rad
    - `gz_at_angle(angle)` -> GZ value at angle
    - `find_max_gz()` -> map with `angle` and `value`
    - `find_angle_of_vanishing_stability()` -> angle or ()
    - `get_first_flooding_angle()` -> angle or ()
    - `find_equilibrium_angle(heeling_arm)` -> angle or ()
    - `find_second_intercept(heeling_arm)` -> angle or ()
    - `get_limiting_angle(default)` -> limiting angle
    
    **Hydrostatic Properties:**
    - `get_gm0()` -> initial GM in meters
    - `get_gm0_dry()` -> GM without FSC
    - `get_draft()` -> draft in meters
    - `get_trim()` -> trim angle in degrees
    - `get_displacement()` -> displacement in kg
    - `get_cog()` -> [x, y, z] center of gravity
    
    **Form Coefficients:**
    - `get_cb()` -> block coefficient
    - `get_cm()` -> midship coefficient  
    - `get_cp()` -> prismatic coefficient
    - `get_lwl()` -> waterline length
    - `get_bwl()` -> waterline breadth
    - `get_vcb()` -> vertical center of buoyancy
    
    **Wind Data:**
    - `has_wind_data()` -> bool
    - `get_emerged_area()` -> area in m²
    - `get_wind_lever_arm()` -> arm in meters
    - `calculate_wind_heeling_lever(pressure)` -> lever at given pressure
    
    **Parameters:**
    - `get_param(key)` -> value or ()
    - `has_param(key)` -> bool
    
    **Metadata:**
    - `get_vessel_name()` -> string
    - `get_loading_condition()` -> string
    """
    
    @staticmethod
    def from_result(
        result: CompleteStabilityResult,
        vessel_name: str,
        loading_condition: str,
    ) -> "CriteriaContext":
        """Create a context from a CompleteStabilityResult.
        
        Args:
            result: The stability calculation result.
            vessel_name: Name of the vessel.
            loading_condition: Description of loading condition.
        """
        ...
    
    def set_param(self, key: str, value: str | float | bool) -> None:
        """Set a parameter accessible to the script.
        
        Args:
            key: Parameter name.
            value: Parameter value (str, float, or bool).
        """
        ...
    
    def get_first_flooding_angle(self) -> float | None:
        """Get the first angle where downflooding occurs, or None."""
        ...
    
    def find_equilibrium_angle(self, heeling_arm: float) -> float | None:
        """Find the first stable equilibrium angle (where GZ = heeling_arm)."""
        ...
    
    def find_second_intercept(self, heeling_arm: float) -> float | None:
        """Find the second intercept angle (unstable equilibrium)."""
        ...


class ScriptEngine:
    """Rhai script execution engine."""
    
    def __init__(self) -> None: ...
    
    def run_script_file(self, path: str, context: CriteriaContext) -> CriteriaResult:
        """Run a Rhai script from file.
        
        Args:
            path: Path to .rhai script file.
            context: Data context for the script.
        
        Returns:
            Verification result.
        
        Raises:
            ValueError: If script execution fails.
        """
        ...
    
    def run_script(self, script: str, context: CriteriaContext) -> CriteriaResult:
        """Run a Rhai script from string.
        
        Args:
            script: Rhai script content.
            context: Data context for the script.
        """
        ...
