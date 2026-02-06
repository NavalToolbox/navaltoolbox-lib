Coordinate System
=================

NavalToolbox uses a **right-handed coordinate system** consistent with standard naval architecture conventions.

Axes Definition
---------------

.. list-table::
   :widths: 15 40 25
   :header-rows: 1

   * - Axis
     - Direction
     - Sign Convention
   * - **X**
     - Longitudinal (length)
     - Positive towards bow (forward)
   * - **Y**
     - Transverse (beam)
     - Positive towards port (left when facing bow)
   * - **Z**
     - Vertical (height)
     - Positive upward

.. note::
   The origin is typically at the intersection of the aft perpendicular (AP), 
   centerline, and baseline (keel), but may vary depending on the hull geometry file.

Rotations
---------

Rotations follow the **right-hand rule**: curl your fingers around the axis in the 
direction of positive rotation, and your thumb points in the positive axis direction.

.. list-table::
   :widths: 15 30 35
   :header-rows: 1

   * - Rotation
     - Axis
     - Positive Direction
   * - **Heel (φ)**
     - Around X-axis
     - Port UP / Starboard down
   * - **Trim (θ)**
     - Around Y-axis
     - Bow down / Stern up

Physical Interpretation
-----------------------

**Heel:**

- **Positive heel** (φ > 0): Port side goes UP, starboard side goes DOWN
- **Negative heel** (φ < 0): Port side goes DOWN, starboard side goes UP

**Trim:**

- **Positive trim** (θ > 0): Bow goes DOWN, stern goes UP
- **Negative trim** (θ < 0): Bow goes UP, stern goes DOWN (stern trim)

Equilibrium Angles
------------------

When calculating equilibrium from a center of gravity (COG) offset:

- **TcG > 0** (weight on port): Results in **negative heel** (port goes down)
- **TcG < 0** (weight on starboard): Results in **positive heel** (starboard goes down)
- **LcG > LCB** (weight forward of buoyancy): Results in **positive trim** (bow down)
- **LcG < LCB** (weight aft of buoyancy): Results in **negative trim** (stern down)

GZ Curve Convention
-------------------

The righting arm GZ is calculated as:

.. math::
   GZ = G_y^{world} - B_y^{world}

Where :math:`G_y^{world}` and :math:`B_y^{world}` are the Y-coordinates of the 
center of gravity and center of buoyancy in the world frame (waterplane horizontal).

- **Positive GZ**: Restoring moment (ship tends to return upright)
- **Negative GZ**: Capsizing moment (ship tends to heel further)

At small heels, this simplifies to: :math:`GZ \approx GM \cdot \sin(\phi)`
