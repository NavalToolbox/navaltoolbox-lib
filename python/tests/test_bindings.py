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
Tests for navaltoolbox Python bindings.

These tests mirror the Rust integration tests to verify
Python bindings work correctly.
"""

import pytest
from pathlib import Path
import tempfile


class TestHullLoading:
    """Tests for Hull STL file loading."""

    def test_load_normal_stl(self):
        """Test loading a normal STL file (DTMB 5415)."""
        from navaltoolbox import Hull

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )

        if not stl_path.exists():
            pytest.skip(f"DTMB 5415 STL not found at {stl_path}")

        hull = Hull(str(stl_path))
        num_triangles = hull.num_triangles()
        bounds = hull.get_bounds()

        # Verify hull loaded successfully
        assert num_triangles > 1000, (
            f"Expected > 1000 triangles, got {num_triangles}"
        )

        # Verify bounds are reasonable for DTMB 5415
        loa = bounds[1] - bounds[0]
        boa = bounds[3] - bounds[2]

        assert 140.0 < loa < 160.0, f"LOA should be ~142-153m, got {loa:.2f}m"
        assert 15.0 < boa < 25.0, f"BOA should be ~19-20m, got {boa:.2f}m"

    def test_load_box_stl(self):
        """Test loading a box STL file."""
        from navaltoolbox import Hull

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "box_10x10.stl"
        )

        if not stl_path.exists():
            pytest.skip(f"Box STL not found at {stl_path}")

        hull = Hull(str(stl_path))
        num_triangles = hull.num_triangles()

        # Box should have a reasonable number of triangles
        assert num_triangles > 0, "Box hull should have triangles"

    def test_empty_stl_file_error(self):
        """Test that loading an empty file raises an appropriate error."""
        from navaltoolbox import Hull

        with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
            temp_path = f.name

        try:
            # Attempt to load empty file should raise an error
            with pytest.raises(Exception) as exc_info:
                _ = Hull(temp_path)

            error_msg = str(exc_info.value).lower()
            # Check that error message is meaningful
            is_meaningful = (
                "empty" in error_msg
                or "stl" in error_msg
                or "failed" in error_msg
            )
            assert is_meaningful, (
                f"Error message should mention empty/STL/failed, "
                f"got: {exc_info.value}"
            )

        finally:
            Path(temp_path).unlink(missing_ok=True)


class TestTank:
    """Tests for Tank class."""

    def test_box_tank_volume(self):
        """Test box tank volume calculation."""
        from navaltoolbox import Tank

        # 10 x 5 x 2 = 100 m³
        tank = Tank.from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0)

        expected_volume = 100.0
        error = abs(tank.total_volume - expected_volume)

        assert error < 1.0, (
            f"Tank volume error: expected {expected_volume}, "
            f"got {tank.total_volume}"
        )

    def test_tank_fill_level(self):
        """Test tank fill level management."""
        from navaltoolbox import Tank

        tank = Tank.from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0)

        # Set 50% fill
        tank.fill_percent = 50.0

        assert abs(tank.fill_level - 0.5) < 1e-6
        assert abs(tank.fill_percent - 50.0) < 1e-6

        # Fill volume should be ~50 m³
        assert (
            abs(tank.fill_volume - 50.0) < 1.0
        ), f"Fill volume should be ~50, got {tank.fill_volume}"

        # Fluid mass should be ~50000 kg
        assert (
            abs(tank.fluid_mass - 50000.0) < 100.0
        ), f"Fluid mass should be ~50000, got {tank.fluid_mass}"

    def test_tank_center_of_gravity(self):
        """Test tank center of gravity calculation."""
        from navaltoolbox import Tank

        tank = Tank.from_box("Test", 0.0, 10.0, -2.5, 2.5, 0.0, 2.0, 1000.0)
        tank.fill_percent = 50.0

        cog = tank.center_of_gravity

        # For a box tank at 50% fill:
        # x should be at center (5.0)
        # y should be at center (0.0)
        # z should be at half of fill height
        assert abs(cog[0] - 5.0) < 0.1, f"CoG x should be ~5.0, got {cog[0]}"
        assert abs(cog[1]) < 0.1, f"CoG y should be ~0.0, got {cog[1]}"
        assert cog[2] < 1.0, f"CoG z should be < 1.0 (half fill), got {cog[2]}"

    def test_tank_free_surface_moment(self):
        """Test free surface moment calculation."""
        from navaltoolbox import Tank

        tank = Tank.from_box("Test", 0.0, 10.0, -2.5, 2.5, 0.0, 2.0, 1000.0)
        tank.fill_percent = 50.0

        fsm_t = tank.free_surface_moment_t
        fsm_l = tank.free_surface_moment_l

        assert fsm_t > 0.0, "FSM_t should be > 0 at partial fill"
        assert fsm_l > 0.0, "FSM_l should be > 0 at partial fill"

        # For box tank:
        # I_t = L * B³ / 12 = 10 * 5³ / 12 ≈ 104.17
        # I_l = B * L³ / 12 = 5 * 10³ / 12 ≈ 416.67
        assert abs(fsm_t - 104.17) < 5.0, (
            f"FSM_t should be ~104.17, got {fsm_t}"
        )
        assert abs(fsm_l - 416.67) < 20.0, (
            f"FSM_l should be ~416.67, got {fsm_l}"
        )

    def test_tank_no_free_surface_when_full(self):
        """Test no free surface when tank is full."""
        from navaltoolbox import Tank

        tank = Tank.from_box("Test", 0.0, 10.0, -2.5, 2.5, 0.0, 2.0, 1000.0)
        tank.fill_percent = 100.0

        assert tank.free_surface_moment_t == 0.0, "FSM_t should be 0 when full"
        assert tank.free_surface_moment_l == 0.0, "FSM_l should be 0 when full"


class TestVessel:
    """Tests for Vessel class."""

    def test_vessel_from_tank_manager(self):
        """Test vessel tank management."""
        from navaltoolbox import Vessel, Tank, Hull

        # We need a hull to create a vessel - use the DTMB5415
        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )
        if not stl_path.exists():
            pytest.skip("DTMB5415 STL file not found")

        hull = Hull(str(stl_path))
        vessel = Vessel(hull)

        assert vessel.num_tanks() == 0

        # Add a tank
        tank = Tank.from_box("Fuel1", 50.0, 60.0, -5.0, 5.0, 0.0, 3.0, 850.0)
        vessel.add_tank(tank)

        assert vessel.num_tanks() == 1


class TestDTMB5415:
    """DTMB5415 validation tests."""

    # Reference GZ data from Ariffin (2017) PhD Thesis Figure 3.15
    REFERENCE_DATA = {
        0: 0.000,
        5: 0.171,
        10: 0.339,
        15: 0.505,
        20: 0.674,
        25: 0.848,
        30: 0.993,
        35: 1.069,
        40: 1.077,
        45: 1.025,
        50: 0.924,
        55: 0.789,
        60: 0.625,
    }

    # Loading condition
    DISPLACEMENT = 8635000.0  # kg
    LCG = 71.670
    TCG = 0.0
    VCG = 7.555

    @pytest.fixture
    def dtmb5415_vessel(self):
        """Load DTMB5415 hull and create vessel."""
        from navaltoolbox import Hull, Vessel

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )
        if not stl_path.exists():
            pytest.skip("DTMB5415 STL file not found")

        hull = Hull(str(stl_path))
        return Vessel(hull)

    def test_dtmb5415_hull_loads(self):
        """Test DTMB5415 hull loads correctly."""
        from navaltoolbox import Hull

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )
        if not stl_path.exists():
            pytest.skip("DTMB5415 STL file not found")

        hull = Hull(str(stl_path))

        assert hull.num_triangles() > 1000, "Hull should have many triangles"

        bounds = hull.get_bounds()
        loa = bounds[1] - bounds[0]
        boa = bounds[3] - bounds[2]

        print(f"DTMB5415: LOA={loa:.2f}m, BOA={boa:.2f}m")

        assert 150.0 < loa < 160.0, f"LOA should be ~153m, got {loa}"
        assert 18.0 < boa < 22.0, f"BOA should be ~20.5m, got {boa}"

    def test_dtmb5415_volume_at_reference_draft(self, dtmb5415_vessel):
        """Test volume at reference draft."""
        from navaltoolbox import HydrostaticsCalculator

        calc = HydrostaticsCalculator(dtmb5415_vessel, 1025.0)
        state = calc.from_draft(6.15, 0.0, 0.0, self.VCG)

        print(f"DTMB5415 at T=6.15m: Volume={state.volume:.1f}m³")

        # Reference: ~8424 m³
        assert 7500.0 < state.volume < 9000.0, (
            f"Volume should be ~8424m³, got {state.volume:.1f}"
        )

    def test_dtmb5415_gz_curve_accuracy(self, dtmb5415_vessel):
        """Test GZ curve accuracy against reference data."""
        from navaltoolbox import StabilityCalculator

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)

        heels = list(self.REFERENCE_DATA.keys())
        cog = (self.LCG, self.TCG, self.VCG)
        curve = calc.gz_curve(self.DISPLACEMENT, cog, heels)

        tolerance = 0.10  # 10cm

        passed = 0
        total = len(heels)

        print("\nDTMB5415 GZ Validation:")
        print("Heel      Calc GZ    Ref GZ     Error")
        print("-" * 45)

        for heel, gz in zip(curve.heels(), curve.values()):
            ref_gz = self.REFERENCE_DATA.get(int(heel), 0)
            error = abs(gz - ref_gz)

            status = "✓" if error < tolerance else "✗"
            if error < tolerance:
                passed += 1

            print(
                f"{heel:5.1f}°    {gz:7.3f}m   "
                f"{ref_gz:7.3f}m   {error:.3f}m {status}"
            )

        pass_rate = passed / total
        print(f"\nPassed: {passed}/{total} ({pass_rate * 100:.0f}%)")

        assert pass_rate >= 0.7, (
            f"At least 70% should pass, got {pass_rate * 100:.0f}%"
        )

    def test_dtmb5415_max_gz_location(self, dtmb5415_vessel):
        """Test maximum GZ location."""
        from navaltoolbox import StabilityCalculator

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)

        heels = list(range(0, 65, 5))
        cog = (self.LCG, self.TCG, self.VCG)
        curve = calc.gz_curve(self.DISPLACEMENT, cog, heels)

        # Find max GZ
        gz_values = curve.values()
        max_gz = max(gz_values)
        max_idx = gz_values.index(max_gz)
        max_heel = heels[max_idx]

        print(
            f"DTMB5415 Max GZ: {max_gz:.3f}m at {max_heel}° "
            f"(Ref: ~1.08m at ~38°)"
        )

        # Max should be between 30-50°
        assert 30.0 <= max_heel <= 50.0, f"Max GZ at wrong angle: {max_heel}°"

        # Max GZ should be ~1.0m
        assert 0.7 < max_gz < 1.5, f"Max GZ should be ~1.0m, got {max_gz:.3f}m"

    def test_kn_curve_binding(self, dtmb5415_vessel):
        """Test KN curve binding with multi-displacement support."""
        from navaltoolbox import StabilityCalculator

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)
        heels = [0, 10, 20]

        # Test with multiple displacements
        displacements = [self.DISPLACEMENT, self.DISPLACEMENT * 0.9]
        curves = calc.kn_curve(
            displacements, heels, lcg=self.LCG, tcg=self.TCG
        )

        assert len(curves) == 2, "Should return one curve per displacement"

        # Verify consistent results with GZ(VCG=0)
        kn_curve_0 = curves[0]

        # Compare with gz_curve(VCG=0)
        gz_curve_0 = calc.gz_curve(
            self.DISPLACEMENT, (self.LCG, self.TCG, 0.0), heels
        )

        for kn_val, gz_val in zip(kn_curve_0.values(), gz_curve_0.values()):
            assert abs(kn_val - gz_val) < 1e-6, "KN should equal GZ(VCG=0)"


class TestWallSidedFormula:
    """Wall-sided GZ formula validation (box barge)."""

    def test_wall_sided_gz(self):
        """
        Validate GZ against exact Wall-Sided Formula.

        Box: L=10, B=10, D=10.
        Draft T = 5.0m.
        KB = T/2 = 2.50m.
        BM = B² / 12T = 100 / 60 = 1.6667m.
        KM = KB + BM = 4.1667m.
        KG = 2.0m.
        GM = KM - KG = 2.1667m.

        Formula: GZ = (GM + 0.5 * BM * tan²(φ)) * sin(φ)
        """
        import math
        from navaltoolbox import Hull, Vessel, StabilityCalculator

        # Create box hull using from_box
        hull = Hull.from_box(10.0, 10.0, 10.0)
        vessel = Vessel(hull)
        stab_calc = StabilityCalculator(vessel, water_density=1025.0)

        # Parameters
        L, B = 10.0, 10.0  # D not used after box creation
        T = 5.0  # Draft
        KG = 2.0  # VCG

        # Calculated values
        KB = T / 2.0  # 2.5
        BM = (B ** 2) / (12.0 * T)  # 1.6667
        GM = KB + BM - KG  # 2.1667

        # Displacement for T=5m
        volume = L * B * T
        displacement = volume * 1025.0

        # COG: centered at (L/2, 0, KG)
        cog = (L / 2.0, 0.0, KG)

        # Calculate GZ curve
        heels = [0.0, 10.0, 20.0, 30.0]
        curve = stab_calc.gz_curve(displacement, cog, heels)

        # Validate against wall-sided formula
        # points() returns [(heel, draft, trim, gz), ...]
        for (heel, draft, trim, gz) in curve.points():
            heel_rad = math.radians(heel)
            tan_phi = math.tan(heel_rad)
            sin_phi = math.sin(heel_rad)

            # Wall-sided formula: GZ = (GM + 0.5 * BM * tan²φ) * sin(φ)
            gz_theoretical = (GM + 0.5 * BM * tan_phi ** 2) * sin_phi

            # Allow some tolerance (numerical integration vs analytical)
            tolerance = 0.05  # 5cm
            error = abs(gz - gz_theoretical)

            assert error < tolerance, (
                f"At heel={heel}°: "
                f"calculated GZ={gz:.4f}m, "
                f"theoretical={gz_theoretical:.4f}m, "
                f"error={error:.4f}m"
            )


class TestWaterplaneCalculations:
    """Waterplane property validation tests."""

    def test_dtmb5415_metacentric_height(self):
        """
        Validate GMT and KMt against DTMB 5415 reference data.

        Reference (Semarak paper, Table 3):
        - Draft: 6.15 m
        - VCG: 7.555 m
        - KMt: 9.493 m
        - GMT: 1.95 m (KMt - VCG)

        Tolerance: ±0.1 m
        """
        from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )

        if not stl_path.exists():
            pytest.skip(f"DTMB 5415 STL not found at {stl_path}")

        hull = Hull(str(stl_path))
        vessel = Vessel(hull)
        calc = HydrostaticsCalculator(vessel, 1025.0)

        # Calculate at reference draft with VCG
        state = calc.from_draft(6.15, vcg=7.555)

        # Calculate KMt
        kmt = state.vcb + state.bmt
        gmt = state.gmt

        # Reference values
        ref_kmt = 9.493
        ref_gmt = 1.95
        tolerance = 0.1

        # Validate KMt
        kmt_error = abs(kmt - ref_kmt)
        assert kmt_error < tolerance, (
            f"KMt error {kmt_error:.3f}m exceeds tolerance {tolerance}m. "
            f"Calculated: {kmt:.3f}m, Reference: {ref_kmt}m"
        )

        # Validate GMT
        gmt_error = abs(gmt - ref_gmt)
        assert gmt_error < tolerance, (
            f"GMT error {gmt_error:.3f}m exceeds tolerance {tolerance}m. "
            f"Calculated: {gmt:.3f}m, Reference: {ref_gmt}m"
        )

        print("\n✅ DTMB 5415 Waterplane Validation:")
        print(f"   KMt: {kmt:.3f}m (ref: {ref_kmt}m, error: {kmt_error:.3f}m)")
        print(f"   GMT: {gmt:.3f}m (ref: {ref_gmt}m, error: {gmt_error:.3f}m)")
        print(f"   Awp: {state.waterplane_area:.1f}m²")

    def test_waterplane_properties_exposed(self):
        """Verify waterplane properties are accessible via bindings."""
        from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )

        if not stl_path.exists():
            pytest.skip("DTMB 5415 STL not found")

        hull = Hull(str(stl_path))
        vessel = Vessel(hull)
        calc = HydrostaticsCalculator(vessel, 1025.0)

        state = calc.from_draft(6.15, vcg=7.555)

        # Check that all waterplane properties exist and are reasonable
        assert hasattr(state, "waterplane_area")
        assert hasattr(state, "lcf")
        assert hasattr(state, "bmt")
        assert hasattr(state, "bml")
        assert hasattr(state, "lwl")
        assert hasattr(state, "bwl")

        # Values should be positive and reasonable
        assert state.waterplane_area > 1000.0, "Awp should be > 1000 m²"
        assert state.bmt > 0.0, "BMt should be positive"
        assert state.bml > 0.0, "BMl should be positive"
        assert state.lwl > 100.0, "LWL should be > 100m"
        assert state.bwl > 10.0, "BWL should be > 10m"

    def test_gmt_gml_optional_without_vcg(self):
        """Verify GMT/GML are None when VCG not provided."""
        from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "dtmb5415.stl"
        )

        if not stl_path.exists():
            pytest.skip("DTMB 5415 STL not found")

        hull = Hull(str(stl_path))
        vessel = Vessel(hull)
        calc = HydrostaticsCalculator(vessel, 1025.0)

        # Calculate without VCG
        state = calc.from_draft(6.15)

        # GMT/GML should be None
        assert state.gmt is None, "GMT should be None when VCG not provided"
        assert state.gml is None, "GML should be None when VCG not provided"
        assert state.cog is None, "COG should be None when VCG not provided"

        # But waterplane properties should still be calculated
        assert state.waterplane_area > 0.0
        assert state.bmt > 0.0
        assert state.bml > 0.0


class TestFreeSurfaceCorrection:
    """Free surface correction validation tests."""

    def test_calc_fsc_box_tank(self):
        """
        Validate FSC calculation for a box tank.

        Tank: 5m x 5m x 2m
        FSM_t = L * B³ / 12 = 5 * 125 / 12 = 52.083 m⁴
        Displacement (Mass) calculated from box hull.
        FSC = FSM / Displacement
        """
        from navaltoolbox import Hull, Vessel, Tank, HydrostaticsCalculator

        # 1. Setup Vessel
        stl_path = (
            Path(__file__).parent.parent.parent
            / "rust"
            / "tests"
            / "data"
            / "box_10x10.stl"
        )
        if not stl_path.exists():
            # Fallback to DTMB if box not available
            stl_path = (
                Path(__file__).parent.parent.parent
                / "rust"
                / "tests"
                / "data"
                / "dtmb5415.stl"
            )

        if not stl_path.exists():
            pytest.skip("No suitable STL found for FSC test")

        hull = Hull(str(stl_path))
        vessel = Vessel(hull)

        # 2. Add Tank
        tank = Tank.from_box(
            name="TestTank",
            x_min=0.0,
            x_max=5.0,
            y_min=-2.5,
            y_max=2.5,
            z_min=0.0,
            z_max=2.0,
            fluid_density=1000.0,
        )
        tank.fill_percent = 50.0
        vessel.add_tank(tank)

        expected_fsm = 5.0 * 5.0**3 / 12.0  # 52.083
        assert abs(tank.free_surface_moment_t - expected_fsm) < 0.1

        # 3. Calculate
        calc = HydrostaticsCalculator(vessel, 1025.0)
        state = calc.from_draft(5.0, vcg=5.0)

        # 4. Validate
        gmt_wet = state.gmt
        gmt_dry = state.gmt_dry

        # Check dry calculation
        kmt = state.vcb + state.bmt
        expected_dry = kmt - 5.0
        assert abs(gmt_dry - expected_dry) < 1e-6, "GMT dry mismatch"

        # Check wet calculation
        # FSC = (FSM * fluid_density) / displacement
        fluid_density = 1000.0  # Tank fluid density in kg/m³
        expected_fsc = (
            tank.free_surface_moment_t * fluid_density / state.displacement
        )
        actual_fsc = gmt_dry - gmt_wet

        assert abs(actual_fsc - expected_fsc) < 1e-6, (
            f"FSC mismatch: expected {expected_fsc}, got {actual_fsc}"
        )
        assert gmt_wet < gmt_dry, "Wet GMT should be less than dry GMT"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
