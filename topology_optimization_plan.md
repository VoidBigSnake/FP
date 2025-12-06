# Stator Topology Optimization Plan (15° Sector)

## Goal
Use the existing immune algorithm (`clonalg.m`) to optimize stator tooth/slot topology within a single 15° sector, then replicate the optimized pattern sixfold around the stator.

## Design Domain Setup
- **Sector scope**: one 15° mechanical sector covering stator back-iron, tooth, and slot area that repeats 6× to form 90° for the 12-slot layout.
- **Symmetry**: build/retain the FEMM model as a quarter (or other) symmetric model, but keep the design variables confined to the 15° sector so that tessellation is straightforward.
- **Discretization options**:
  - **Pixel/voxel grid**: subdivide the 15° sector into a polar grid (r–θ) of cells; each cell is binary (material vs. void/air).
  - **Parameterized primitives**: define a handful of radial/arc control lines to shape tooth walls and slot opening, encoded as continuous parameters.
  - **Hybrid**: use a coarse binary grid for slot interior plus continuous parameters for tooth/bridge thickness to reduce dimensionality.
- **Boundary constraints**: lock the outer radius (stator yoke) and inner airgap circle; enforce minimum tooth/bridge thickness to keep mesh solvable and manufacturable.

## Encoding for CLONALG
- **Binary genes**: map each grid cell to one bit (1 = iron, 0 = air/copper/void). For a 20×8 (θ×r) grid, that is 160 bits. Extend the genome to include continuous parameters by dedicating bit fields per parameter and decoding them with the existing `decode` pattern.
- **Fitness function**: replace `myfunc` to accept the decoded material map/parameters and:
  1. Regenerate the FEMM geometry for the 15° sector by toggling cells/primitives (e.g., drawing blocks or deleting regions) and assigning materials.
  2. Run `scantorque` (or a fast surrogate) to obtain average torque and ripple.
  3. Apply constraints/penalties (e.g., minimum area, demag flux density, manufacturability) and compute scalar fitness such as `J = -T_avg + λ·T_ripple + μ·penalties`.
- **Population size**: keep total bit-length manageable (<256 bits) so hypermutation remains effective and runtime per evaluation stays feasible.

## FEMM Automation Notes
- Prepare a **base template** FEMM file with coils, magnets, boundary conditions, and periodic sector references already defined.
- Script the **geometry rebuild** within `femmfunc` (or a new helper) to:
  - Clear only the design-domain entities.
  - Draw or delete cells according to the bit grid (use `mi_drawpolygon`/`mi_addblocklabel` with material assignment).
  - Re-mesh the sector (`mi_createmesh`) and solve.
- After solving the single-sector model, rely on periodic boundaries; if torque is integrated over the sector only, scale by the number of sectors represented when reporting full-machine torque.

## Symmetry Replication
- For visualization or final validation, replicate the optimized 15° sector six times around the stator (and mirror as needed for full 360°) using rotation copy commands in FEMM.
- Keep the optimization loop on the single sector to minimize per-iteration solve time; only replicate for final checks.

## Runtime & Acceleration
- **Caching**: memoize evaluations for identical genomes to avoid repeated solves.
- **Surrogates**: if FEMM is too slow, train a regression model (e.g., RBF or small neural net) on sampled designs to approximate fitness, then refine with FEMM for elites.
- **Parallelization**: distribute fitness evaluations across CPU cores if FEMM batch mode allows.

## Next Steps
1. Choose a grid resolution and map it to bit indices in `clonalg` (adjust chromosome length accordingly).
2. Prototype a geometry-builder function that reads a bit grid and modifies the FEMM sector template.
3. Replace `myfunc` with a wrapper that decodes the genome, rebuilds the FEMM sector, runs torque/ripple evaluation, and returns the scalar fitness.
4. Validate on a few hand-crafted patterns to ensure the torque scaling and boundary conditions are correct before running full optimization.
