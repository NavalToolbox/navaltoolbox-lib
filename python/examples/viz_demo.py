
import os
import sys
from pathlib import Path

# Add parent directory to path to allow running from examples folder without installing
sys.path.append(str(Path(__file__).parent.parent))

from navaltoolbox import Vessel, Hull, Tank, Silhouette, DownfloodingOpening, OpeningType
from navaltoolbox.visualization import plot_hydrostatic_condition

def main():
    # 1. Resolve path to DTMB STL
    # Assuming standard repo structure: python/examples/viz_demo.py -> rust/tests/data/dtmb5415.stl
    current_dir = Path(__file__).parent
    repo_root = current_dir.parent.parent
    stl_path = repo_root / "rust" / "tests" / "data" / "dtmb5415.stl"

    if not stl_path.exists():
        print(f"Error: STL file not found at {stl_path}")
        # Fallback for different CWD
        stl_path = Path("rust/tests/data/dtmb5415.stl")
        if not stl_path.exists():
             print("Could not locate dtmb5415.stl. Please check working directory.")
             return

    print(f"Loading hull from {stl_path}...")
    hull = Hull(str(stl_path))
    
    # Scale to valid size if it's the model scale version? 
    # Usually DTMB5415 LBP is ~142m full scale or ~5.72m model scale.
    # Let's check bounds.
    bounds = hull.get_bounds()
    length = bounds[1] - bounds[0]
    print(f"Hull Length: {length:.2f} m")

    if length < 10.0:
        print("Scaling up to ship scale (approx 1:24.8)...")
        hull.scale(24.8) # Approx scale for 5.72m model to 142m ship

    # 2. Setup Vessel
    vessel = Vessel(hull)
    vessel.ap = 0
    vessel.fp = 142

    print("Adding components...")

    # 3. Add Tanks
    # Fwd Fuel Tank
    tank1 = Tank.from_box(
        "Fwd Fuel",
        x_min=110.0, x_max=120.0,
        y_min=-4.0, y_max=4.0,
        z_min=4.0, z_max=8.0,
        fluid_density=850.0
    )
    tank1.fill_percent = 85.0
    vessel.add_tank(tank1)

    # Aft Ballast
    tank2 = Tank.from_box_hull_intersection(
        hull=hull,
        x_min=10.0, x_max=20.0,
        y_min=-8.0, y_max=8.0,
        z_min=0.0, z_max=10.0,
        fluid_density=1025.0,
        name="Aft Ballast"
    )
    tank2.fill_percent = 50.0
    vessel.add_tank(tank2)

    # 4. Add Silhouette (Fake superstructure)
    # Assuming deck at Z=10m roughly
    deck_z = bounds[5]
    if length < 10.0:
        deck_z *= 24.8

    points = [
        (10.0, deck_z),
        (80.0, deck_z),
        (80.0, deck_z + 8.0), # Bridge
        (60.0, deck_z + 8.0),
        (60.0, deck_z + 15.0), # Mast
        (55.0, deck_z + 15.0),
        (55.0, deck_z + 8.0),
        (10.0, deck_z + 4.0), # Poop deck
        (10.0, deck_z), # Close silhouette
    ]
    sil = Silhouette.from_points(points, "Profile")
    vessel.add_silhouette(sil)

    # 5. Add Openings
    vent = DownfloodingOpening.from_point(
        "Vent Stbd",
        (40.0, 6.0, deck_z + 1.0),
        OpeningType.vent()
    )
    vessel.add_opening(vent)

    from navaltoolbox.visualization import plot_vessel_3d

    # ... (Keep existing setup) ...

    # 6. Visualize Upright
    print("Generating 3D plot (Upright)...")
    fig = plot_vessel_3d(
        vessel, 
        title="DTMB 5415 - Visualization Demo",
        opacity_hull=0.5
    )

    output_file = "dtmb_demo.html"
    fig.write_html(output_file)
    print(f"Success! Visualization saved to: {os.path.abspath(output_file)}")
    
    # 7. Visualize Inclined (Hydrostatic Condition)
    print("Generating 3D plot (Inclined)...")
    # Simulate a condition: Draft 6.15m, Heel 15 deg, Trim 1 deg
    fig_hydro = plot_hydrostatic_condition(
        vessel,
        draft=6.15,
        heel=15.0,
        trim=1.0,
        title="DTMB 5415 - Hydrostatic (T=6.15m, Heel=15°, Trim=1°, COG shown)",
        opacity_hull=0.5,
        cog=(70.0, 0.0, 7.5), # Example COG
        show_axes=False
    )
    
    output_hydro = "dtmb_hydro.html"
    fig_hydro.write_html(output_hydro)
    print(f"Success! Hydrostatic visualization saved to: {os.path.abspath(output_hydro)}")

    os.system(f"open {output_file}")
    os.system(f"open {output_hydro}")

if __name__ == "__main__":
    main()
