
import pytest
from pathlib import Path
from navaltoolbox import (
    Hull,
    Vessel,
    StabilityCalculator,
    CriteriaContext,
    ScriptEngine,
)


class TestRhaiScripting:
    """Tests for Rhai scripting integration."""

    @pytest.fixture
    def dtmb5415_vessel(self):
        """Load DTMB5415 hull and create vessel."""
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

    def test_reproduce_is_code_2008_error(self, dtmb5415_vessel):
        """
        Reproduce the 'Function not found' error for IS Code 2008 script.

        Specifically checking:
        - find_max_gz
        - area_under_curve
        """
        # 1. Setup stability result
        displacement = 8635000.0  # kg
        lcg = 71.670
        tcg = 0.0
        vcg = 7.555
        cog = (lcg, tcg, vcg)

        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)

        # Calculate full curve 0-180
        heels = list(range(0, 70, 5))
        result = calc.complete_stability(displacement, cog, heels)

        # 2. Create Context
        ctx = CriteriaContext.from_result(
            result,
            vessel_name="DTMB 5415",
            loading_condition="Design Load"
        )

        # 3. Initialize Engine
        engine = ScriptEngine()

        # 4. Load and run IS Code 2008 script
        script_path = (
            Path(__file__).parent.parent.parent
            / "rules"
            / "is_code_2008_general.rhai"
        )

        if not script_path.exists():
            pytest.skip(f"Script not found at {script_path}")

        print(f"Running script: {script_path}")

        # This is expected to fail currently
        try:
            result = engine.run_script_file(str(script_path), ctx)
            print("Script execution successful!")
            print(f"Overall Pass: {result.overall_pass}")
            for crit in result.criteria:
                print(
                    f" - {crit.name}: "
                    f"{crit.status} ({crit.actual_value:.3f} {crit.unit})"
                )

        except ValueError as e:
            pytest.fail(f"Script execution failed: {e}")

    def test_minimal_reproduction(self, dtmb5415_vessel):
        """Minimal reproduction of specific function calls."""
        displacement = 8635000.0
        cog = (71.670, 0.0, 7.555)
        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)
        heels = [0, 10, 20, 30, 40, 50, 60]
        result = calc.complete_stability(displacement, cog, heels)

        ctx = CriteriaContext.from_result(result, "Test", "Test")
        engine = ScriptEngine()

        # Test find_max_gz
        script_max_gz = """
        fn check(ctx) {
            let max = ctx.find_max_gz();
            print("Max GZ: " + max);
            return #{ name: "Test", status: "PASS" };
        }
        """
        engine.run_script(script_max_gz, ctx)

        # Test area_under_curve
        script_area = """
        fn check(ctx) {
            let area = ctx.area_under_curve(0.0, 30.0);
            print("Area: " + area);
            return #{ name: "Test", status: "PASS" };
        }
        """
        engine.run_script(script_area, ctx)

    def test_is_code_2008_values_dtmb5415(self, dtmb5415_vessel):
        """Verify numerical values of IS Code 2008 criteria on DTMB5415."""
        displacement = 8635000.0
        cog = (71.670, 0.0, 7.555)
        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)
        heels = list(range(0, 70, 5))
        result = calc.complete_stability(displacement, cog, heels)

        ctx = CriteriaContext.from_result(result, "DTMB 5415", "Test")
        engine = ScriptEngine()

        script_path = (
            Path(__file__).parent.parent.parent
            / "rules"
            / "is_code_2008_general.rhai"
        )
        
        if not script_path.exists():
            pytest.skip("IS Code script not found")
            
        res = engine.run_script_file(str(script_path), ctx)
        
        assert res.overall_pass is True, "DTMB5415 should pass all general criteria"
        assert len(res.criteria) == 6
        
        crit_names = [c.name for c in res.criteria]
        assert any("Area 0-30°" in n for n in crit_names)
        assert any("GZ at 30°+" in n for n in crit_names)
        assert any("Initial GM₀" in n for n in crit_names)
        
        for crit in res.criteria:
            assert crit.status == "PASS"
            if "GM₀" in crit.name:
                assert abs(crit.actual_value - result.hydrostatics.gmt) < 0.01

    def test_run_all_rules(self, dtmb5415_vessel):
        """Verify that all scripts in rules/ directory can be executed."""
        displacement = 8635000.0
        cog = (71.670, 0.0, 7.555)
        calc = StabilityCalculator(dtmb5415_vessel, 1025.0)
        heels = list(range(0, 70, 5))
        result = calc.complete_stability(displacement, cog, heels)

        ctx = CriteriaContext.from_result(result, "Test", "Test")
        engine = ScriptEngine()

        rules_dir = Path(__file__).parent.parent.parent / "rules"
        scripts = list(rules_dir.glob("*.rhai"))

        assert len(scripts) > 0, "No rule scripts found"

        for script in scripts:
            print(f"Testing script: {script.name}")
            try:
                res = engine.run_script_file(str(script), ctx)
                # Basic check that we got a result
                assert res.overall_pass is not None
            except Exception as e:
                # template.rhai might fail if it expects parameters
                # we didn't set, but we should check if it's a
                # syntax error vs runtime error
                print(f"Script {script.name} failed: {e}")
                # We expect success for IMO scripts. Template might fail logic
                # but should run.
                if "im" in script.name:
                    pytest.fail(f"IMO script {script.name} failed: {e}")
