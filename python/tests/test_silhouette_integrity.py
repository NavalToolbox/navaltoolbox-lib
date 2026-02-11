import unittest
import tempfile
import os
import navaltoolbox as nt
import logging

# Configure logging to capture warnings during tests
logging.basicConfig(level=logging.WARNING)


class TestSilhouetteIntegrity(unittest.TestCase):

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if os.path.exists(f):
                os.unlink(f)

    def create_temp_file(self, content, suffix):
        tf = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
        tf.write(content)
        tf.close()
        self.temp_files.append(tf.name)
        return tf.name

    def create_dxf(self, points, closed=True):
        content = [
            "0", "SECTION", "2", "ENTITIES",
            "0", "LWPOLYLINE", "8", "0",
            "90", f"{len(points)}",
            "70", f"{1 if closed else 0}"
        ]
        for x, y in points:
            content.extend(["10", f"{x}", "20", f"{y}"])
        content.extend(["0", "ENDSEC", "0", "EOF"])
        return self.create_temp_file("\n".join(content), ".dxf")

    def test_csv_loading_xz(self):
        """Test loading X, Z format from CSV."""
        content = "0.0, 0.0\n10.0, 0.0\n10.0, 5.0\n0.0, 5.0"
        filepath = self.create_temp_file(content, ".csv")
        s = nt.Silhouette(filepath)
        self.assertAlmostEqual(s.get_area(), 50.0, places=2)
        self.assertTrue(s.is_closed())
        self.assertEqual(s.num_points(), 5)  # 4 + auto-close

    def test_txt_loading_xyz(self):
        """Test loading X, Y, Z format from TXT (Y ignored)."""
        content = "0.0 1.0 0.0\n10.0 2.0 0.0\n10.0 3.0 5.0\n0.0 4.0 5.0"
        filepath = self.create_temp_file(content, ".txt")
        s = nt.Silhouette(filepath)
        self.assertAlmostEqual(s.get_area(), 50.0, places=2)
        # Should project to X-Z plane, ignoring Y
        bounds = s.get_bounds()  # x_min, x_max, z_min, z_max
        self.assertAlmostEqual(bounds[1], 10.0)
        self.assertAlmostEqual(bounds[3], 5.0)

    def test_csv_comments_and_headers(self):
        """Test CSV with comments and headers."""
        content = "# Comment\nX, Z\n0.0, 0.0\n5.0, 5.0\n10.0, 0.0"
        filepath = self.create_temp_file(content, ".csv")
        s = nt.Silhouette(filepath)
        self.assertEqual(s.num_points(), 4)  # 3 + auto-close
        self.assertAlmostEqual(s.get_area(), 25.0, places=2)

    def test_dxf_closed_loading(self):
        """Test loading a closed DXF silhouette."""
        points = [(0, 0), (10, 0), (10, 5), (0, 5)]
        filepath = self.create_dxf(points, closed=True)
        s = nt.Silhouette(filepath)
        self.assertTrue(s.is_closed())
        self.assertAlmostEqual(s.get_area(), 50.0, places=2)

    def test_dxf_open_loading(self):
        """Test loading an open DXF silhouette."""
        points = [(0, 0), (10, 0), (10, 5)]
        filepath = self.create_dxf(points, closed=False)
        s = nt.Silhouette(filepath)
        self.assertFalse(s.is_closed())
        # Area calculation for open loop might differ based on shoelace
        # implementation checking it loads is key here

    def test_validation_zero_area(self):
        """Test that zero-area silhouette creation doesn't crash."""
        points = [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]
        s = nt.Silhouette.from_points(points, "ZeroArea")
        self.assertAlmostEqual(s.get_area(), 0.0)

    def test_stale_calculator_fix(self):
        """Verify that adding silhouette after calculator creation works."""
        hull = nt.Hull.from_box(20.0, 5.0, 2.0)
        vessel = nt.Vessel(hull)
        calc = nt.StabilityCalculator(vessel, 1025.0)

        # Add silhouette AFTER calculator creation
        s = nt.Silhouette.from_points([
            (0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0), (0.0, 0.0)
        ], "Test")
        vessel.add_silhouette(s)

        disp = 20.0 * 5.0 * 1.0 * 1025.0
        res = calc.complete_stability(disp, (10.0, 0.0, 1.0), [0.0])
        self.assertTrue(res.has_wind_data(),
                        "Calculator failed to see new silhouette")
        self.assertGreater(res.wind_data.emerged_area, 0.0)

    def test_stability_without_silhouette(self):
        """Verify behavior when no silhouette is present."""
        hull = nt.Hull.from_box(20.0, 5.0, 2.0)
        vessel = nt.Vessel(hull)
        calc = nt.StabilityCalculator(vessel, 1025.0)

        disp = 20.0 * 5.0 * 1.0 * 1025.0
        res = calc.complete_stability(disp, (10.0, 0.0, 1.0), [0.0])
        self.assertFalse(res.has_wind_data())


if __name__ == '__main__':
    unittest.main()
