# blend.fillet_corner_p3: the closed Turn reports the blend corner where C++ builds the closed Nx1 corner (2026-11-05)

Units: `blend.fillet_corner_p1`, `blend.fillet_corner_p2`, `blend.fillet_corner_p3`,
`blend.fillet_executive_p3`. PR: #2349's closed Turn.

## Change

A closed blended edge meeting itself at an angle is the blend corner's Turn, Mitred: an
`NByTwo` corner with one solver edge use and `BlendCornerEvent::Turn`, as the closed Join is
for one meeting itself G1. Pinned C++ (v1.0.1) builds `SmFilletNx1ClosedCorner` there.
`FilletCornerKind::NByOneClosed` and `topology/closed_nx1.rs` are deleted. The construction
is the old corner's, step for step. The filleted #2349 teardrop is unchanged bit for bit
(vertices and edge samples at full precision, radii 0.5, 1, 1.5 and 3, default and
generic-blend kernels).

## Differences from C++, recorded as known mismatches

- The corner's type: C++ `SM_FCR_N_x_1_CLOSED` (7), Rust `SM_FCR_BLEND` (9), and no
  corner keeps the closed Nx1 corner's `m_bG1Case` (Rust -1, C++ 0):
  `blend.fillet_corner_p3` gen0011 and gen0332 (`SmFilletCorner::Create` on
  `Fillet/rollover_corner.smb`).
- The #2349 teardrop's executive checkpoints (corners, first topology, first vertex
  geometry, solver geometry and corner edges at r 0.5; the first three at r 3): the corner
  type only. Every other value of each checkpoint agrees.

The checkpoint runner now builds the first corner's topology by its class
(`FilletSession::make_corner_topology`, the virtual `MakeCornerTopology` the shim calls)
in place of the Nx1 corner's.

## Retired cases

Cases that build `SmFilletNx1ClosedCorner` directly (the shim's class index 2) have no
Rust class left to compare and are retired, as #2400 retired the 1x1 and 2x2 corners':
16 in `blend.fillet_corner_p1`, 44 in `blend.fillet_corner_p2` (which runs part 1's
cases) and 5 in `blend.fillet_corner_p3`. The 2 known mismatches #2400 recorded for the
closed Nx1 corner's `GetTangentFilletSolver` are among them.

## Unchanged on main

Before this change, on main 4b10d7d6e and on the rebased #2381 head 9aa43f084,
`blend.fillet_executive_p3` already fails preview-radius-tiny,
source-file-box_w_slot-all-edges, source-file-cyl_on_box_w_fillet-all-edges and
sweep4-sphere_diff_3cylinders-all-r0.1 (main's CI report, run 41980837473). This change
does not touch them.
