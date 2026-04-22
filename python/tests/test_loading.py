"""Tests for LoadingCondition feature."""

import json
import tempfile
from pathlib import Path


class TestMassCategory:
    """Tests for MassCategory enum."""

    def test_create_categories(self):
        from navaltoolbox import MassCategory
        ls = MassCategory.lightship()
        dw = MassCategory.deadweight()
        ot = MassCategory.other()
        assert repr(ls)  # Should not raise
        assert repr(dw)
        assert repr(ot)


class TestMassItem:
    """Tests for MassItem class."""

    def test_create_with_defaults(self):
        from navaltoolbox import MassItem
        item = MassItem("Test", 1000.0, (10.0, 0.0, 5.0))
        assert item.name == "Test"
        assert item.mass == 1000.0
        assert item.cog == (10.0, 0.0, 5.0)

    def test_create_with_category(self):
        from navaltoolbox import MassItem, MassCategory
        item = MassItem(
            "Lightship",
            5_000_000.0,
            (45.0, 0.0, 4.5),
            MassCategory.lightship(),
        )
        assert item.name == "Lightship"
        assert item.mass == 5_000_000.0


class TestLoadingCondition:
    """Tests for LoadingCondition class."""

    def test_create(self):
        from navaltoolbox import LoadingCondition
        lc = LoadingCondition("Test Loading")
        assert lc.name == "Test Loading"
        assert lc.num_masses() == 0
        assert lc.num_tank_overrides() == 0

    def test_add_masses(self):
        from navaltoolbox import LoadingCondition, MassCategory
        lc = LoadingCondition("Test")
        lc.add_mass_simple(
            "Lightship", 5_000_000, (45.0, 0.0, 4.5), MassCategory.lightship()
        )
        lc.add_mass_simple("Crew", 3_000, (35.0, 0.0, 8.0))
        assert lc.num_masses() == 2
        masses = lc.get_masses()
        assert len(masses) == 2
        assert masses[0].name == "Lightship"
        assert masses[1].name == "Crew"

    def test_remove_mass(self):
        from navaltoolbox import LoadingCondition
        lc = LoadingCondition("Test")
        lc.add_mass_simple("A", 1000, (0.0, 0.0, 0.0))
        lc.add_mass_simple("B", 2000, (0.0, 0.0, 0.0))
        assert lc.num_masses() == 2
        assert lc.remove_mass("A") is True
        assert lc.num_masses() == 1
        assert lc.remove_mass("X") is False

    def test_tank_fill_overrides(self):
        from navaltoolbox import LoadingCondition
        lc = LoadingCondition("Test")
        lc.set_tank_fill_percent("FO_1P", 95.0)
        lc.set_tank_fill("FW_1", 0.5)
        assert lc.num_tank_overrides() == 2
        fills = lc.get_tank_fills()
        assert abs(fills["FO_1P"] - 0.95) < 1e-6
        assert abs(fills["FW_1"] - 0.5) < 1e-6

    def test_apply_and_resolve(self):
        """Test apply/resolve with a vessel that has tanks."""
        from navaltoolbox import LoadingCondition, Hull, Vessel, Tank

        hull = Hull.from_box(length=100.0, breadth=20.0, depth=10.0)
        vessel = Vessel(hull)

        # Add a tank: 10x10x2 = 200 m³, density 1000 kg/m³
        tank = Tank.from_box("FO_1P", 20.0, 30.0, -5.0, 5.0, 0.0, 2.0, 1000.0)
        tank.fill_percent = 0.0  # Start empty
        vessel.add_tank(tank)

        lc = LoadingCondition("Test")
        lc.add_mass_simple("Lightship", 100_000, (50.0, 0.0, 5.0))
        lc.set_tank_fill_percent("FO_1P", 50.0)

        # Apply and resolve
        lc.apply(vessel)

        # Tank should now be 50% filled
        tanks = vessel.get_tanks()
        assert abs(tanks[0].fill_percent - 50.0) < 1e-6

        disp, cog = lc.resolve(vessel)
        # Total = 100_000 + tank_mass (~100_000 kg for 100m³ at 1000 kg/m³)
        assert disp > 100_000  # Must include tank mass
        assert cog[0] > 0  # LCG should be positive

        # Test item-only resolution (excluding tanks)
        item_disp, item_cog = lc.resolve_items()
        assert item_disp == 100_000
        assert item_cog[0] == 50.0

    def test_unaffected_tanks_keep_fill(self):
        """Tanks not in the override list must keep their current fill."""
        from navaltoolbox import LoadingCondition, Hull, Vessel, Tank

        hull = Hull.from_box(length=100.0, breadth=20.0, depth=10.0)
        vessel = Vessel(hull)

        tank1 = Tank.from_box("FO_1P", 20.0, 30.0, -5.0, 0.0, 0.0, 2.0, 1000.0)
        tank1.fill_percent = 80.0
        vessel.add_tank(tank1)

        tank2 = Tank.from_box("FW_1", 40.0, 50.0, -5.0, 5.0, 0.0, 2.0, 1000.0)
        tank2.fill_percent = 60.0
        vessel.add_tank(tank2)

        # Only override FO_1P, FW_1 should stay at 60%
        lc = LoadingCondition("Test")
        lc.add_mass_simple("Lightship", 100_000, (50.0, 0.0, 5.0))
        lc.set_tank_fill_percent("FO_1P", 10.0)
        lc.apply(vessel)

        tanks = vessel.get_tanks()
        assert abs(tanks[0].fill_percent - 10.0) < 1e-6  # FO_1P changed
        assert abs(tanks[1].fill_percent - 60.0) < 1e-6  # FW_1 unchanged

    def test_json_serialization(self):
        """Test JSON round-trip."""
        from navaltoolbox import LoadingCondition, MassCategory

        lc = LoadingCondition("Departure")
        lc.add_mass_simple(
            "Lightship", 5_000_000, (45.0, 0.0, 4.5), MassCategory.lightship()
        )
        lc.add_mass_simple("Crew", 3_000, (35.0, 0.0, 8.0))
        lc.set_tank_fill_percent("FO_1P", 95.0)
        lc.set_tank_fill("FW_1", 0.5)

        # Serialize
        json_str = lc.to_json()
        data = json.loads(json_str)
        assert data["name"] == "Departure"
        assert len(data["masses"]) == 2
        assert "FO_1P" in data["tank_fills"]

        # Deserialize
        lc2 = LoadingCondition.from_json(json_str)
        assert lc2.name == "Departure"
        assert lc2.num_masses() == 2
        assert lc2.num_tank_overrides() == 2

    def test_json_file_io(self):
        """Test JSON file save/load."""
        from navaltoolbox import LoadingCondition

        lc = LoadingCondition("Test File")
        lc.add_mass_simple("Test", 1000.0, (10.0, 0.0, 5.0))

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name

        try:
            lc.save_json(path)
            lc2 = LoadingCondition.load_json(path)
            assert lc2.name == "Test File"
            assert lc2.num_masses() == 1
        finally:
            Path(path).unlink(missing_ok=True)

    def test_copy(self):
        """Test copy with optional new name."""
        from navaltoolbox import LoadingCondition

        lc = LoadingCondition("Original")
        lc.add_mass_simple("A", 1000, (0.0, 0.0, 0.0))
        lc.set_tank_fill_percent("FO_1P", 95.0)

        # Copy without name -> keeps original
        lc2 = lc.copy()
        assert lc2.name == "Original"
        assert lc2.num_masses() == 1

        # Copy with new name
        lc3 = lc.copy("Arrival")
        assert lc3.name == "Arrival"
        assert lc3.num_masses() == 1
        assert lc3.num_tank_overrides() == 1

    def test_set_name(self):
        from navaltoolbox import LoadingCondition
        lc = LoadingCondition("Original")
        lc.name = "Modified"
        assert lc.name == "Modified"

    def test_from_loading_restores_fills(self):
        """Test that from_loading methods restore the tank fills."""
        from navaltoolbox import LoadingCondition, Hull, Vessel, Tank, HydrostaticsCalculator

        hull = Hull.from_box(length=100.0, breadth=20.0, depth=10.0)
        vessel = Vessel(hull)
        tank = Tank.from_box("FO_1P", 20.0, 30.0, -5.0, 5.0, 0.0, 2.0, 1000.0)
        tank.fill_percent = 25.0
        vessel.add_tank(tank)

        lc = LoadingCondition("Test")
        lc.add_mass_simple("Lightship", 100_000, (50.0, 0.0, 5.0))
        lc.set_tank_fill_percent("FO_1P", 100.0)

        calc = HydrostaticsCalculator(vessel)
        _ = calc.from_loading(lc)

        # Tank should be restored to 25.0%
        assert abs(vessel.get_tanks()[0].fill_percent - 25.0) < 1e-6

    def test_gz_curve_from_loading(self):
        """Test gz_curve_from_loading matches manual resolve_items workflow."""
        from navaltoolbox import LoadingCondition, Hull, Vessel, Tank, StabilityCalculator

        hull = Hull.from_box(length=100.0, breadth=20.0, depth=10.0)
        vessel = Vessel(hull)
        tank = Tank.from_box("FO_1P", 20.0, 30.0, -5.0, 5.0, 0.0, 2.0, 1000.0)
        vessel.add_tank(tank)

        lc = LoadingCondition("Test")
        lc.add_mass_simple("Lightship", 3_000_000, (50.0, 0.0, 5.0))
        lc.set_tank_fill_percent("FO_1P", 50.0)

        calc = StabilityCalculator(vessel)
        heels = [0.0, 10.0, 20.0]
        
        curve = calc.gz_curve_from_loading(lc, heels)
        assert len(curve.values()) == 3

    def test_from_loading_hydrostatics(self):
        """Test from_loading matches manual resolve workflow for hydrostatics."""
        from navaltoolbox import LoadingCondition, Hull, Vessel, Tank, HydrostaticsCalculator

        hull = Hull.from_box(length=100.0, breadth=20.0, depth=10.0)
        vessel = Vessel(hull)
        tank = Tank.from_box("FO_1P", 20.0, 30.0, -5.0, 5.0, 0.0, 2.0, 1000.0)
        vessel.add_tank(tank)

        lc = LoadingCondition("Test")
        lc.add_mass_simple("Lightship", 5_000_000, (50.0, 0.0, 5.0))
        lc.set_tank_fill_percent("FO_1P", 80.0)

        calc = HydrostaticsCalculator(vessel)
        state1 = calc.from_loading(lc)
        
        # Manual apply and resolve
        lc.apply(vessel)
        disp, cog = lc.resolve(vessel)
        state2 = calc.from_displacement(displacement_mass=disp, cog=cog)

        assert abs(state1.draft - state2.draft) < 1e-6
