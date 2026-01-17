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

These tests mirror the Rust integration tests to verify Python bindings work correctly.
"""

import pytest
import math
import os
from pathlib import Path


class TestTank:
    """Tests for Tank class."""

    def test_box_tank_volume(self):
        """Test box tank volume calculation."""
        from navaltoolbox import Tank

        # 10 x 5 x 2 = 100 m³
        tank = Tank.from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0)

        expected_volume = 100.0
        error = abs(tank.total_volume - expected_volume)

        assert error < 1.0, f"Tank volume error: expected {expected_volume}, got {tank.total_volume}"

    def test_tank_fill_level(self):
        """Test tank fill level management."""
        from navaltoolbox import Tank

        tank = Tank.from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0)

        # Set 50% fill
        tank.fill_percent = 50.0

        assert abs(tank.fill_level - 0.5) < 1e-6
        assert abs(tank.fill_percent - 50.0) < 1e-6

        # Fill volume should be ~50 m³
        assert abs(tank.fill_volume - 50.0) < 1.0, f"Fill volume should be ~50, got {tank.fill_volume}"

        # Fluid mass should be ~50000 kg
        assert abs(tank.fluid_mass - 50000.0) < 100.0, f"Fluid mass should be ~50000, got {tank.fluid_mass}"

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
        assert abs(fsm_t - 104.17) < 5.0, f"FSM_t should be ~104.17, got {fsm_t}"
        assert abs(fsm_l - 416.67) < 20.0, f"FSM_l should be ~416.67, got {fsm_l}"

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
        stl_path = Path(__file__).parent.parent.parent / "rust" / "tests" / "data" / "dtmb5415.stl"
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

        stl_path = Path(__file__).parent.parent.parent / "rust" / "tests" / "data" / "dtmb5415.stl"
        if not stl_path.exists():
            pytest.skip("DTMB5415 STL file not found")

        hull = Hull(str(stl_path))
        return Vessel(hull)

    def test_dtmb5415_hull_loads(self):
        """Test DTMB5415 hull loads correctly."""
        from navaltoolbox import Hull

        stl_path = Path(__file__).parent.parent.parent / "rust" / "tests" / "data" / "dtmb5415.stl"
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
        state = calc.calculate_at_draft(6.15, 0.0, 0.0, self.VCG)

        print(f"DTMB5415 at T=6.15m: Volume={state.volume:.1f}m³")

        # Reference: ~8424 m³
        assert 7500.0 < state.volume < 9000.0, f"Volume should be ~8424m³, got {state.volume:.1f}"

    def test_dtmb5415_gz_curve_accuracy(self, dtmb5415_vessel):
        """Test GZ curve accuracy against reference data."""
        from navaltoolbox import StabilityCalculator

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)

        heels = list(self.REFERENCE_DATA.keys())
        curve = calc.calculate_gz_curve(
            self.DISPLACEMENT,
            (self.LCG, self.TCG, self.VCG),
            heels
        )

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

            print(f"{heel:5.1f}°    {gz:7.3f}m   {ref_gz:7.3f}m   {error:.3f}m {status}")

        pass_rate = passed / total
        print(f"\nPassed: {passed}/{total} ({pass_rate * 100:.0f}%)")

        assert pass_rate >= 0.7, f"At least 70% should pass, got {pass_rate * 100:.0f}%"

    def test_dtmb5415_max_gz_location(self, dtmb5415_vessel):
        """Test maximum GZ location."""
        from navaltoolbox import StabilityCalculator

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)

        heels = list(range(0, 65, 5))
        curve = calc.calculate_gz_curve(
            self.DISPLACEMENT,
            (self.LCG, self.TCG, self.VCG),
            heels
        )

        # Find max GZ
        gz_values = curve.values()
        max_gz = max(gz_values)
        max_idx = gz_values.index(max_gz)
        max_heel = heels[max_idx]

        print(f"DTMB5415 Max GZ: {max_gz:.3f}m at {max_heel}° (Ref: ~1.08m at ~38°)")

        # Max should be between 30-50°
        assert 30.0 <= max_heel <= 50.0, f"Max GZ at wrong angle: {max_heel}°"

        # Max GZ should be ~1.0m
        assert 0.7 < max_gz < 1.5, f"Max GZ should be ~1.0m, got {max_gz:.3f}m"


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
        # This test requires creating a box mesh programmatically
        # Since we can't do that easily in Python, we skip this test
        # The equivalent test is run in Rust
        pytest.skip("Box mesh creation not available in Python bindings")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
