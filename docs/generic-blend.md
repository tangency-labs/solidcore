# The generic rolling-ball blend

This page describes the blend that is replacing SMLib's fillet code step by step: what it builds,
how, why it is built that way, what is merged, and what is open. It is written for a reader who
knows blending surfaces and is new to this code. Numbers in brackets are references (end of page).
The same page lives in the SMLib fork (tangency-labs/smlib, `docs/generic-blend.md`); this copy adds
how the work is ported to solidcore (§9).

Contents:
1. [Purpose and the one rule](#1-purpose-and-the-one-rule)
2. [The construction](#2-the-construction)
3. [Events of the trace](#3-events-of-the-trace)
4. [The exact classes](#4-the-exact-classes)
5. [Design decisions that shape behaviour](#5-design-decisions-that-shape-behaviour)
6. [Status](#6-status)
7. [Case gallery](#7-case-gallery)
8. [For contributors](#8-for-contributors)
9. [Porting to solidcore](#9-porting-to-solidcore)
10. [References](#10-references)

---

## 1. Purpose and the one rule

**One construction for every pair of faces.** A round, a curvature-continuous round and a chamfer
of any edge, between faces of any kind, are built by the same steps: a spine, two rails, a section,
and the events met along the way. Corners, ends and overflows are events of that one trace, not
separate algorithms chosen by the shape of the input.

**Failures are fixed in that construction.** A case that fails is a defect of the construction,
fixed where the construction is wrong. It is never fixed by a point solution for that case, a
fallback to the old code, or a chain of methods tried in turn until one returns something. The old
fillet code had all three: an analytic try-chain for the surface, fourteen vertex kinds and twelve
corner classes chosen by edge counts, and a march into tangent points. Each step of the rework
deletes one of them.

The other rules follow from that one:
- every result is exact or certified (§2.4); nothing is fitted and then hoped to be close enough;
- a blend that cannot be built is refused with the size that would build (§3.5); the user's
  geometry is never changed silently;
- the design is grounded in the published literature, cited here.

The order of the work follows the survey's own order of the problem [46]: spine and rails, then
section, then vertex blends.

## 2. The construction

![The rolling ball between two faces: the offsets meet in the spine; the ball's feet are the rails; the section joins them.](generic-blend/rolling-ball.svg)

*Figure 1. One cross-section of a round between faces F1 and F2: the offsets O1 and O2 meet in the
spine m; the ball of radius r centred on m touches each face at its rail; the section joins the rails.*

### 2.1 Spine and rails

Take a blend of edge e between faces F1 and F2, with outward normals n1, n2 and radius r, constant
or a law r(t).

- **Spine.** m = O1 ∩ O2, where Oi = Fi + r·ni is face i offset towards the blend [6, 37]. With a
  radius law the offsets move with r(t); the spine is defined implicitly and traced by marching with
  the law [7]. The solver's existing station solve (Newton on both offsets) is kept; what changes is
  what is made from its points.
- **Rails.** ρi(t) = m(t) − r(t)·ni at the foot point: the ball's feet. A spine point's parameters
  on Oi are its foot's parameters on Fi, so the rails cost no further solve. The rails are all the
  rolling ball fixes [6, 14].
- **Cut.** Each face is trimmed on its rail. The pcurve is the spine's parameters on the offset.
- **Chamfer.** The same construction with the contacts at the setback distances, joined by the chord.

Every geometric quantity downstream (rails, sections, events, corner boundaries) is built from the
same exact spine points.

### 2.2 The section

The section joins the two rails. It is one of three kinds, chosen by its order n:

| n | section | surface |
|---|---|---|
| 0 | the chord | the chamfer |
| 1 | the circle of the ball between the rails | the exact canal surface, the envelope of the balls [35] |
| 2 to 4 | the Gⁿ section | Hartmann's blend between the rails [14, 17] |

![One cross-section of a square edge with the chord, the circle and the order-2 section on the same rails.](generic-blend/gn-sections.svg)

*Figure 2. The section family on one cross-section, both rails set back s: chord (n = 0), circle
(n = 1), order-2 section (n = 2, drawn without thumbweight).*

**The circle (n = 1)** is the default round. It is stored as an exact rational circle: the rational
quadratic arc with middle weight cos(θ/2), θ the angle between the contact normals [6, 7]. It is
never an approximated arc.

**The Gⁿ section (n ≥ 2)** is Hartmann's construction [14, 17]: two curves, one on each face from
its rail towards the edge, blended by f(t) = (1 − t)ⁿ⁺¹ / ((1 − t)ⁿ⁺¹ + tⁿ⁺¹). The first n
derivatives of f vanish at both ends, so the section has Cⁿ contact with each face along its rail,
hence Gⁿ. Its n = 1 member is not the circle, so the circle stays the order-1 section and Hartmann's
family is used from n = 2 up. Two parameters shape it:
- **Thumbweight w** [17, 46], the fullness. Its default is **equal area**: at each station, the w at
  which the section removes the same cross-section area as the circle of that radius would there
  ((1 − π/4) r² at a square edge). Each station takes its own w from its own wedge angle, so the
  curvature round removes about the same material as the round of the same radius everywhere along
  the edge.
- **The faces' order-n jets.** The section is built from curves on each face whose first n − 1
  coefficients are the feet's and whose n-th is fixed by the face's n-th order data in the normal
  direction only (the tangential part is free under reparametrization [10]). So a Gⁿ round needs a
  Cⁿ face, not Cⁿ⁺¹: a G2 round is continuous on a cubic B-spline face, and the span certificate
  holds across its knot lines.

The section's contact order with its face is n + 1 in general, and n + 2 between planes. That order
gives the tangent-end laws of §3.1.

A **variable round** is the union of its balls (Figure 3). The ball touches the envelope along the
circle in the plane through the rails square to the centre line's tangent c′, not along the great
circle: its centre is shifted by −r r′/|c′| along the centre line and its radius is
r √(1 − r′²/|c′|²) [35, 7, 29]. The great-circle section is off the envelope whenever r′ ≠ 0 and can
cut through the balls [29]. The envelope is real only where |c′| > |r′| [35].

![Side view of balls along a centre line, their envelope, the touch circle and the great circle.](generic-blend/variable-round.svg)

*Figure 3. A variable round is the envelope of its balls. The section is the touch circle (accent),
shifted from the ball's centre; the great circle (dashed) is not on the envelope.*

### 2.3 The round stored as its definition

A blend face is stored as the system that defines it, following Hoffmann's dimensionality paradigm
[21]: more unknowns and simple equations, nothing eliminated, evaluated on demand [1, 41]. The class
is `SmProceduralBlend`, which derives from `SmBSplineSurface` the way the analytic surfaces carry
their identity on a net. It holds:

| field | content |
|---|---|
| supports S1, S2 | the support faces' surfaces, by value; a support may itself be a definition, nested to depth 3 |
| sides | ±1, the ball's side of each support |
| law | constant r, or a curve r(s) over the leader's parameter |
| order, thumbweight | n (1 circle, 2 to 4 Gⁿ) and w |
| leader | the centre line the stations are solved in; it fixes the u parameter |
| seeds | the stations' parameter curves, to start Newton |
| ε | the bound of the certified net |

Evaluation solves one station per point (Newton from the seeds, then derivatives order by order by
the implicit function theorem) and the section in closed form along it. The Jacobian is singular
exactly at the events of §3, so a definition never evaluates through an event silently.

The net B is a certified approximation of the definition S, kept beside it. B serves bounding boxes,
the surface cache, tessellation and display; final intersection points, projections, the model check
and continuity grading are solved on S. B is never an operand, so error does not stack under later
booleans: every edge on a blend face lies within its own tolerance of S, however many operations
follow [22, 20]. Only export flattens to B and carries ε, once.

Rounds of the exact classes (§4) are stored as their analytic or rational surface, with no
definition.

### 2.4 Certification

"Within tolerance at the span midpoints" is a sampled check, not a certificate [11]. The net is
certified per span instead: the error e = S − B is known to machine precision on the span's grid,
and between grid points it is bounded by the interpolation remainder with second derivatives
enclosed rigorously, S's from interval evaluation over the span box. Existence and uniqueness of the
station's root over the span box are proved by Moore's interval test on the Krawczyk operator
[32, 28, 33]: K(X) ⊆ X proves a root in X, and K(X) ⊆ X with ‖I − Y F′(X)‖ < 1 proves it unique.
Knots are inserted until the bound reaches a tenth of the zone tolerance; the achieved bound is
stored as ε. Where plain intervals over-estimate, the remedy is a tighter enclosure [3, 4, 31], never
a sampled check.

## 3. Events of the trace

The spine runs until it meets something known exactly. Each meeting is an event; each event is a
root of a scalar or a small system along the trace, solved where it is well conditioned.

| the spine meets | event | result |
|---|---|---|
| an end face's offset: the ball touches the end face E | End | the blend ends on E; its end curve is S ∩ E |
| an end face tangent to F along a seam T | tangent end | the end curve built outward from P (§3.1) |
| another spine | corner | Turn, Join, Ball or Setback (§3.2) |
| its face's boundary: the rail reaches a rim | overflow | the blend kept and trimmed at its crease (§3.4) |
| a face that cannot hold the ball | fold or room | refused, naming the size that fits (§3.3, §3.5) |

### 3.1 The tangent end

Where an end face E runs tangent into the rail's face F along a seam T, both normals agree at the
rail's end P, and marching into P fails (the intersection residual's Jacobian vanishes there). So P
is solved as an event: the rail crossing T in F's parameters, a transversal curve–curve intersection
solved by Newton from the bracketing stations.

In coordinates at P, a along the rail and s along T, E stands off F by ½ Δκ a² and the section by
c sᵏ with k = n + 1. The end curve therefore follows **a ∝ s^(k/2)**:

| order | end curve at P |
|---|---|
| chamfer, n = 0 | a conic tangent to the rail |
| round, n = 1 | a straight branch of the node at P, from the second-order tangent cone |
| curvature round, n = 2 | tangent to T, curvature growing as s^(−1/2) |

![Plan view of the tangent end: the rail, the seam and the three end laws from P.](generic-blend/tangent-end.svg)

*Figure 4. The end curve leaves P with a ∝ s^(1/2), s, s^(3/2) for the chamfer, the round and the
curvature round. Each is built from P outward, never marched into P.*

No regular parametrization has unbounded curvature, so the order-2 end is held in t = √s: the end
curve and the round's end patch are B-splines in t with zero speed at P. In the regularized chart
v = t², a = α tᵏ, the end condition becomes ĝ(α, t) = ½ Δκ α² − c + O(t), regular at t = 0 with the
simple root α* = ±√(2c/Δκ), so Newton converges quadratically down to P. The tangent-end vertex is
first class: curve-dropping and face orientation take the first non-zero derivative there, and the
model check accepts the s^(−1/2) curvature growth at that vertex only. Where E osculates F (Δκ = 0),
the end is refused. A round fading to nothing (its two contacts merging) is singular in the round's
own system and is refused until its own chart is derived.

### 3.2 Corners

One class, `SmBlendCorner`, classifies every blended vertex once and builds one of these events. It
never returns "no kind": a configuration the table does not name is a Setback.

| event | when | construction |
|---|---|---|
| **End** | one blend ends at the vertex | the end rib on the end face, S ∩ E; a tangent side is the tangent end of §3.1 |
| **Join** | two blends meeting G1 and continuing each other | one shared section |
| **Turn** | two blends sharing a face, turning at a sharp unrounded side edge (a boss base, a pocket's floor or rim) | **Mitre** (default) or **Swept** (option) |
| **Ball** | three equal rounds of one convexity at a vertex of three faces | the corner ball's sphere |
| **Setback** | everything else: unequal radii, four or more blends, mixed convexity, and any corner where G2 is selected | the setback hole and its certified convex fill |

![Plan views of the four corner kinds: mitre crease, swept torus, ball sphere patch, setback hole.](generic-blend/corner-events.svg)

*Figure 5. Corner events. Mitre: the rounds run on and are cut on their crossing, a crease on the
bisector. Swept: the section revolved about the side edge. Ball: one rational patch on the corner
ball's sphere. Setback: the rounds stop short and a 2n-sided hole is filled.*

- **Mitre.** Each round runs on to the other and both are trimmed on their crossing. Two equal
  cylinders meet in a planar conic on the bisector plane, built exactly. The seam is G0 by
  construction and graded G0.
- **Swept.** Each round ends on its own section in the plane through the side edge square to its
  edge; the section is revolved about the side edge through the corner's turn. Between planes at a
  straight side edge this is an exact torus segment (a horn torus touching the axis at the side
  vertex), G1 to both rounds. It is the ball kept touching the side edge: the offset of an edge is a
  cylinder about it, so the spine there is a circle and the envelope a torus [37, 6, 18, 8]. A ball
  cannot pass a sharp concave side edge, so Swept there is refused, naming the corner radius that
  would let it pass (R ≥ r).
- **Ball.** Three equal rounds at three faces end on the corner ball's sphere [37, 6]. The corner
  face is always one rational patch exactly on the sphere, of degree (4, 2) by stereographic
  construction [8], with no analytic sphere special case: one code path.
- **Setback** [44, 43]. Each blend stops at a setback distance from the vertex. The hole has 2n
  sides: n profile curves (the blends' end sections) alternating with n spring curves on the faces.
  Plain trimlines that simply meet give poor shapes for unequal radii, edges far from 90° and
  convex–concave mixes [44]. The fill:
  - follows the setback split: one setback patch per blend plus an n-sided interior region. The
    central split of the 2n-sided hole has a twist incompatibility for 2n = 6, 8, 10; the setback
    split does not [44];
  - takes its boundary rows exactly from the rounds' end jets; each setback patch meets the
    interior through a connecting map, which is Gᵏ exactly [13];
  - makes the interior one polynomial map of a regular n-gon, so its parts meet each other as one
    surface and there is no vertex-enclosure condition at the centre [47, 12];
  - is chosen as the least thin-plate energy among the fills that meet every condition exactly, and
    certified convex: at a convex corner the fill stays inside the rounds' prisms with K ≥ 0
    (mirrored at concave corners). Convexity is required because rounding a convex body with
    constant balls is an opening by a ball, which keeps it convex;
  - uses no least squares and no averaging.

  The default setback is the larger of 1.5 r and the contact-line minimum plus 0.25 r. Where that
  does not certify, the automatic setback grows to the smallest common scale that certifies (found by
  bisection; certification is monotone in scale), up to 2×. A setback the user gives explicitly that
  does not fit or does not certify is refused, naming the size that does. No uncertified fill is ever
  returned. At a convex edge meeting concave ones, any setback blend dips under the faces' planes
  regardless of representation [44]; there the spring curve between adjacent convex and concave
  blends is dropped and their trimlines meet.

**Corners in a radius law.** A law with a corner (r continuous, r′ not) is a corner of the union of
balls. At a dip the corner ball's sphere zone between the two sides' touch circles lies on the
round, G1 to both sides; it is built as an exact rational biquadratic patch. At a peak the two
sides' envelopes cross and each ends on the crossing curve, a crease.

### 3.3 Folds and room: the checks

The checks are Lukács's, stated exactly [29, 30]:
1. **Each face locally encloses the ball.** At every contact, h_i(r) = (E − rL)(G − rN) − (F − rM)² > 0,
   with (E, F, G) and (L, M, N) the face's fundamental forms. h_i = 0 means r is a principal radius
   of curvature there (a focal point) and the offset degenerates.
2. **No local self-intersection.** Every point of the section arc is progressive (f_tt > 0), Lukács's
   condition (4.28).
3. **Variable rounds use the envelope**, never the great-circle sweep (§2.2).

A transversality check (spine against the Voronoi surface) was built and then shown redundant with
the station solve: it can fire only where the two offsets are tangent, where the solve already finds
no section. It was deleted.

A blend past a check is refused, naming the radius where the check first fails. It is never built
folded and never trimmed back.

### 3.4 Overflow: trimmed at the crease, in place

Where a rail reaches the rim of its face (a cliff), the blend is kept and the sharp edge moves: the
blend is built on its faces continued past the rim and trimmed along its crease with the face beyond,
S ∩ X. For a plane X the crease is exact (a line on a cylinder, a conic on a cone); otherwise it is
certified. This replaces building a blend and subtracting the overflow afterwards.

![Cross-section of a round whose rail would land past a narrow top face, trimmed where it meets the cliff wall.](generic-blend/overflow-crease.svg)

*Figure 6. Keep the blend: the round keeps its radius and is trimmed in place at its crease with the
continued cliff face, so the rim moves down.*

Keeping the edge instead (the ball rolling on the rim, the spine continuing as the other offset
meeting the rim's pipe) is a possible future option, never the default.

### 3.5 Refusals name the size that fits

![Left: offsets of a concave face collapse at its curvature radius. Right: a ball in a slot fits up to half the width.](generic-blend/refusal-size.svg)

*Figure 7. Two refusal causes. Fold: the offset collapses at the face's curvature radius. Room: the
ball must fit between the faces.*

Every refusal (no room, fold, too large, overflow that cannot be trimmed, a setback that does not fit)
returns the largest size that builds, computed from the cause's own geometry: in closed form where
one exists, else by solving that cause's geometry, never by building again. Where the cause's bound is
open (a size at it would leave a face or an edge of no extent), it is the bound less a thousand of the
edge's tolerances. The C API returns it next to the reason:

```c
int32_t smlib_blend_reason(char *out, int64_t size, int64_t *out_length);
int32_t smlib_blend_fit_size(double *out_size, int32_t *out_what); /* SMLIB_FIT_SIZE, _SETBACK or _SCALE */
```

The application shows the refusal and offers that size; it does not search by retrying builds.

## 4. The exact classes

Exact rational fixed-radius rolling-ball blends exist for exactly the configurations Dahl and
Krasauskas classify [9]; everything else gets the certified spine. One
classification replaced the old try-and-fall-through chain.

| faces | spine | surface, minimal bidegree |
|---|---|---|
| plane/plane | a line | the cylinder |
| plane/cylinder, plane square to a coaxial cylinder's axis, coaxial pairs | a circle | the torus |
| plane/cone at any angle (plane/cylinder as a limit) | a conic or a line | (4, 2) [9, Thms 6, 8, 9] |
| two cones with two touching points | a pair of conics | (6, 2) [9, Thms 16–18] |
| two cones or cylinders with one touching point | a singular quartic | (12, 2) [9, Thms 19–21] |
| cone/sphere (sphere/cylinder) with one touching point | a singular quartic | (8, 2) [9, Thm 22] |

A rational spine with a rational radius gives a real rational canal surface [35, Thm 5.1], and a
rolling ball on a surface with rational offsets along a rational rail is rational [35, §4]; minimal
degrees with the two rails as boundary are in [27], and the variable-radius cylinder and cone cases
in [24]. A pair whose offsets happen to touch at one radius is not taken as exact unless proved for
that radius. Between planes, the Gⁿ section is rational too (the feet are linear), and a Gⁿ section
with rotational symmetry is that section revolved; both are stored as rational surfaces.

## 5. Design decisions that shape behaviour

Each line is a decision the code holds to. Where common CAD practice is documented, the comparison is
given.

- **Default round: the circular G1 round.** The curvature round (G2) and orders up to 4 are
  selectable, never the default. Matches common CAD practice.
- **Round sections are exact rational circles**, not approximated arcs; expectations moved with
  them are re-derived, not loosened.
- **G2 shape: equal-area thumbweight per station**, with no user shape factor. Common CAD practice
  either offers a shape factor or fixes an undocumented shape.
- **Turn corners: Mitre by default, Swept as an option**; Swept at a sharp concave side edge refused
  with the radius that would pass. A crease by default with a rounded turn as an option matches the
  documented practice of most systems.
- **Three or more rounds of unequal radii at a vertex: Setback with a certified convex fill by
  default.** Common practice offers setback corners as an option the user sets up.
- **Default setback: max(1.5 r, contact-line minimum + 0.25 r)**, grown automatically to the
  smallest scale that certifies; an explicit setback that does not fit is refused, naming the size.
  Common practice has the user enter it.
- **Mitre creases are graded G0.**
- **Ball corners are one rational patch on the sphere**, box corners included.
- **Overflow keeps the blend and moves the sharp edge**, matching the constant-radius overflow mode
  common CAD practice offers; keeping the edge is a future option.
- **A boss or cut standing inside a round's band is kept by default** (its walls carried down to the
  round); an option lets the round take it. Matches the "keep features" option of common practice.
- **A partly selected tangent chain** rounds exactly the selected edges and ends at the smooth vertex
  with a real End.
- **Folds are refused, never trimmed**, naming the radius.
- **The user's geometry is never changed without consent**: no silent shrink, no auto-fit.
- **The round is stored as its definition** with a certified net.
- **No Chebyshev or Müntz bases for blends**: B-splines, in t = √s at a tangent end.
- **Parallel structure, serial for now**: chains, spines, stations, spans and corners are pure
  functions over immutable inputs with results in a fixed order, so the output is bit-identical
  whatever the thread count once threading is switched on.

## 6. Status

The rework lands in seven steps, each behind tests committed failing first.

| step | content | state |
|---|---|---|
| 1 | the tangent end as an event; the end rib built outward, order 2 in t = √s; the march into P deleted | merged, [smlib #76](https://github.com/tangency-labs/smlib/pull/76) |
| 2 | one classification of the exact rounds, each built as its rational surface; the try-chain deleted | merged, [smlib #77](https://github.com/tangency-labs/smlib/pull/77) |
| 3 | one Gⁿ section generator (chord, circle, curvature) certified per span; variable rounds as the union of their balls | merged, [smlib #78](https://github.com/tangency-labs/smlib/pull/78) |
| 4 | corners: Turn (Mitre, Swept), Join, Ball, End merged; Setback with its convex fill, mixed convexity and the remaining old corner classes open | part merged, [smlib #79](https://github.com/tangency-labs/smlib/pull/79) |
| 5 | fold and room checks, exact; refusals sized through `smlib_blend_fit_size` | merged in #79; remaining classes open |
| 6 | overflow trimmed at the crease, in place | merged in #79 |
| 7 | the round stored as its definition, with a certified net | merged in #79; nesting and the fading end open |

The code is in `source/SMLib/src`: `SmBlendCorner*.cpp` (classification and the corner events),
`SmFilletSection.cpp` (`SmGnSectionFSG`, the section generator), `SmFilletSolver.cpp` (the trace
and the §3.3 checks), `SmFilletOverflow.cpp` (the crease trim), `SmFilletExact.cpp` (the exact
classes), `SmProceduralBlend.cpp` (the definition), `SmBallPatch.cpp` (the Ball's sphere patch) and
`SmBlendSetback.cpp` (the setback fill); the C API is in `capi/include/smlib.h`.

**Open work**, each an issue to claim:

| issue | topic |
|---|---|
| [#80](https://github.com/tangency-labs/smlib/issues/80) | Setback with a certified convex fill (G1) and setback growth |
| [#81](https://github.com/tangency-labs/smlib/issues/81) | a round ending at the smooth vertex of a partly selected tangent chain |
| [#82](https://github.com/tangency-labs/smlib/issues/82) | mixed convexity: convex and concave blends at one vertex (the test peg) |
| [#83](https://github.com/tangency-labs/smlib/issues/83) | delete the old corner hierarchy |
| [#84](https://github.com/tangency-labs/smlib/issues/84) | several concave rounds at one vertex, unequal or N ≥ 4 |
| [#85](https://github.com/tangency-labs/smlib/issues/85) | a ball or setback corner with a curved face |
| [#86](https://github.com/tangency-labs/smlib/issues/86) | four or more equal rounds at a vertex, and open setbacks |
| [#87](https://github.com/tangency-labs/smlib/issues/87) | the G2 setback fill held convex |
| [#88](https://github.com/tangency-labs/smlib/issues/88) | four blends over three faces |
| [#89](https://github.com/tangency-labs/smlib/issues/89) | chamfers meeting at a vertex |
| [#90](https://github.com/tangency-labs/smlib/issues/90) | a closed edge meeting itself at an angle |
| [#91](https://github.com/tangency-labs/smlib/issues/91) | size the remaining refusal classes; refuse the collapsed offset |
| [#92](https://github.com/tangency-labs/smlib/issues/92) | a round on a definition face; the fading end |
| [#93](https://github.com/tangency-labs/smlib/issues/93) | the pipe-on-plate curvature round grades a seam G3 where G2 is wanted |
| [#94](https://github.com/tangency-labs/smlib/issues/94) | a tangent chain over an earlier round's curved end curve fails |
| [#95](https://github.com/tangency-labs/smlib/issues/95) | building near a refusal's size bound fails |
| [#96](https://github.com/tangency-labs/smlib/issues/96) | performance of the order-2 section certificate and of a grown setback |
| [#97](https://github.com/tangency-labs/smlib/issues/97) | follow-ups: refitted frames, exact plane∩cylinder/torus edges, pcurves, untested paths |
| [#98](https://github.com/tangency-labs/smlib/issues/98) | a rail that rolls over a tangent edge before a smooth-vertex End |

## 7. Case gallery

One render per fillet and chamfer case of the test suites, with its status: **passes**, **refused**
(with the size SMLib names), or **known failing** (with its reason and issue). The cases are grouped
by the event they exercise:

- ends: plain ends and tangent ends;
- corners: Turn (Mitre, Swept), Join, Ball, Setback, law corners;
- overflow;
- refusals;
- variable rounds;
- Gⁿ sections.

The renders are produced with the SolidCraft command-line tool, one PNG per case, named
`<case>.png`.

<!-- The image list is filled from the render folder: one row per case,
     | render | case | status |, under the group headings below. -->

### Ends

### Corners: Turn

### Corners: Join

### Corners: Ball

### Corners: Setback

### Corners in a radius law

### Overflow

### Refusals

### Variable rounds

### Gⁿ sections

## 8. For contributors

- **Claim first.** Pick an open issue from §6, say on it that you are taking it, and keep one topic
  per pull request.
- **Tests first.** Commit the failing test before the fix: a C API case in
  `capi/tests/ops_kernel.cpp` (or the known-failing list it joins), gated on values a wrong result
  would change: exact volumes, distances to the exact surface, seam angles, counts of faces and
  spans. Timings are machine-dependent and are never gates.
- **Root causes only.** Find why the construction is wrong and fix that. No special case for the
  failing input, no fallback to old code, no try-chain, no loosened expectation. Several failing
  corner cases usually share one cause; fix the cause.
- **Literature first.** Before designing a construction or a fix, check the literature; most of
  these problems have a published analysis. This page's bibliography is where references live; code
  comments and commit messages describe what the code does, not which paper it came from.
- **Exact or certified.** Every surface and curve you add is exact (analytic or rational) or
  certified against its definition with a stored bound. A sampled check is not a certificate.
- **Refuse with a size.** A new refusal returns the largest size that builds through
  `smlib_blend_fit_size`, or says why there is none.
- **Every commit builds**, and the tests it touches pass or are listed known-failing with a reason.
- **Delete what the new construction replaces** in the same step, once nothing calls it.

## 9. Porting to solidcore

The blend is designed and landed in the SMLib fork first; solidcore ports each merged step. The port
is **bit-identical** to the fork: the same inputs give the same output, bit for bit, and the difftest
(`crates/difftest`, Rust and C++ run side by side on the same inputs) checks it unit by unit. A
difference is a port defect to find, never a tolerance to widen. Fixing a defect in the construction
happens in the fork first, then is ported.

| port issue | fork step |
|---|---|
| [solidcore #2243](https://github.com/tangency-labs/solidcore/issues/2243) | smlib #76, step 1: the tangent end as an event of the trace |
| [solidcore #2244](https://github.com/tangency-labs/solidcore/issues/2244) | smlib #77, step 2: one classification of the exact rounds |
| [solidcore #2250](https://github.com/tangency-labs/solidcore/issues/2250) | smlib #78, step 3: one Gⁿ section generator; variable rounds as the union of their balls |
| [solidcore #2273](https://github.com/tangency-labs/solidcore/issues/2273) | smlib #79, steps 4 (part) to 7: corners, fold and room checks with fit sizes, overflow in place, the round as its definition |

Port the steps in order: each fork step deletes code the next one assumes gone. Claim the port issue
before starting, and keep the fork's tests as the port's tests.

## 10. References

Titles, pages or DOIs are given only where the project's reading notes record them; an entry without one is incomplete, not unpublished.

1. R. E. Barnhill, G. Farin, Q. Chen. Constant-radius blending of parametric surfaces. *Computing Suppl.* 8, Springer (1993) 1–20. doi:10.1007/978-3-7091-6916-2_1
2. R. E. Barnhill, S. N. Kersey. A marching method for parametric surface/surface intersection. *Computer Aided Geometric Design* 7(1–4) (1990) 257–280. doi:10.1016/0167-8396(90)90035-P
3. Bühler. On validated enclosures of parametric objects (title, venue and year not recorded). doi:10.1145/584458.584479
4. Bühler, Barth. On validated enclosures of parametric objects (title, venue and year not recorded). doi:10.1007/978-1-4757-6484-0_15
5. B. K. Choi. *Surface Modeling for CAD/CAM*. Elsevier (1991).
6. B. K. Choi, S. Y. Ju. Constant-radius blending in surface modelling. *Computer-Aided Design* 21(4) (1989) 213–220. doi:10.1016/0010-4485(89)90046-8
7. J.-H. Chuang, C.-H. Lin, W.-C. Hwang. Variable-radius blending of parametric surfaces. *The Visual Computer* 11 (1995) 513–525. doi:10.1007/BF02434038
8. H. E. I. Dahl. Piecewise rational parametrizations of canal surfaces. *Mathematical Methods for Curves and Surfaces 2012*, LNCS 8177, Springer (2014) 88–111. doi:10.1007/978-3-642-54382-1_6
9. H. E. I. Dahl, R. Krasauskas. Rational fixed radius rolling ball blends between natural quadrics. *Computer Aided Geometric Design* 29 (2012) 691–706. doi:10.1016/j.cagd.2012.06.008
10. DeRose. Technical report UCB/CSD-86-255, University of California, Berkeley (1985).
11. Filkins, Tuohy, Patrikalakis. *Engineering with Computers* 9 (1993) 49–62. doi:10.1007/BF01198253
12. J. A. Gregory, J. M. Hahn. A C² polygonal surface patch. *Computer Aided Geometric Design* 6(1) (1989) 69–75. doi:10.1016/0167-8396(89)90007-1
13. Hahn. Geometric continuous patch complexes. Brunel University technical report TR/05/87 (1987).
14. E. Hartmann. Gⁿ-blending with rolling ball contact curves. *Proc. Geometric Modeling and Processing 2000*, IEEE (2000) 385–389. doi:10.1109/GMAP.2000.838269
15. E. Hartmann. Numerical parameterization of curves and surfaces. *Computer Aided Geometric Design* 17(3) (2000) 251–266. doi:10.1016/S0167-8396(99)00050-3
16. E. Hartmann. Implicit Gⁿ-blending of vertices. *Computer Aided Geometric Design* 18 (2001) 267–285.
17. E. Hartmann. Parametric Gⁿ blending of curves and surfaces. *The Visual Computer* 17(1) (2001) 1–13. doi:10.1007/PL00013398
18. E. Hartmann. Gⁿ-continuous connections between normal ringed surfaces. *Computer Aided Geometric Design* 18 (2001) 751–770.
19. G. J. Hettinga, J. Kosinka. Multisided B-spline patches over extraordinary regions. *STAG 2020*, Eurographics (2020). doi:10.2312/stag.20201235
20. Hoffmann. The problems of accuracy and robustness in geometric computation. *IEEE Computer* 22(3) (1989). doi:10.1109/2.16223
21. Hoffmann. A dimensionality paradigm for surface interrogations. *Computer Aided Geometric Design* 7 (1990) 517–532. doi:10.1016/0167-8396(90)90013-H
22. Jackson. Boundary representation modelling with local tolerances. *ACM Symposium on Solid Modeling* (1995). doi:10.1145/218013.218067
23. S. K. Jaiswal, A. Ghosal, B. Gurumoorthy. An efficient algorithm for circular blending of edges. *ASME DETC 2001*, DETC01/DAC-21065 (2001).
24. M. Kazakevičiūtė, R. Krasauskas. Blending cylinders and cones using canal surfaces. *Nonlinear Analysis: Modelling and Control* 5 (2000) 77–89.
25. K. Kim, G. Elber. A symbolic approach to freeform parametric surface blends. University of Utah technical report (1995).
26. Kos, Martin, Várady. *Computer Aided Geometric Design* 17 (2000). doi:10.1016/S0167-8396(99)00043-6
27. R. Krasauskas. Minimal rational parametrizations of canal surfaces. *Computing* 79 (2007) 281–290. doi:10.1007/s00607-006-0204-0
28. Krawczyk. *Computing* 4 (1969). doi:10.1007/BF02234767
29. G. Lukács. Differential geometry of G1 variable radius rolling ball blend surfaces. *Computer Aided Geometric Design* 15 (1998) 585–613. doi:10.1016/S0167-8396(98)00006-5
30. G. Lukács, T. Hermann, T. Várady. Continuity and self-intersections of variable radius rolling ball blend surfaces. *Product Modeling for Computer Integrated Design and Manufacture*, IFIP (1997) 369–370. doi:10.1007/978-0-387-35187-2_31
31. Makino, Berz. *International Journal of Pure and Applied Mathematics* 4 (2003).
32. Moore. *SIAM Journal on Numerical Analysis* 14 (1977) 611–615. doi:10.1137/0714040
33. Neumaier. *Interval Methods for Systems of Equations* (1990). doi:10.1017/CBO9780511526473
34. N. M. Patrikalakis, T. Maekawa. *Shape Interrogation for Computer Aided Design and Manufacturing*. Springer (2002).
35. M. Peternell, H. Pottmann. Computing rational parametrizations of canal surfaces. *Journal of Symbolic Computation* 23 (1997) 255–266. doi:10.1006/jsco.1996.0087
36. K. Qin, Y. Li, C. Deng. Blending Bézier patch for multi-sided surface modeling. *Computer Aided Geometric Design* 105 (2023) 102222. doi:10.1016/j.cagd.2023.102222
37. J. R. Rossignac, A. A. G. Requicha. Constant-radius blending in solid modeling. *Computers in Mechanical Engineering* 3 (1984) 65–73.
38. P. Salvi. A circular parameterization for multi-sided patches. arXiv:2305.07534 (2023).
39. P. Salvi, T. Várady. G2 surface interpolation over general topology curve networks. *Computer Graphics Forum* 33(7) (2014) 151–160. doi:10.1111/cgf.12483
40. M. A. Sanglikar, P. Koparkar, V. N. Joshi. Modelling rolling ball blends for computer aided geometric design. *Computer Aided Geometric Design* 7(5) (1990) 399–414.
41. Stroud. *Boundary Representation Modelling Techniques*. Springer (2006). doi:10.1007/978-1-84628-616-2
42. T. Várady, C. M. Hoffmann. Vertex blending: problems and solutions. *Mathematical Methods for Curves and Surfaces II*, Vanderbilt University Press (1998) 501–527.
43. T. Várady, A. P. Rockwood. Vertex blending based on the setback split. *Mathematical Methods for Curves and Surfaces*, Vanderbilt University Press (1995).
44. T. Várady, A. Rockwood. Geometric construction for setback vertex blending. *Computer-Aided Design* 29(6) (1997) 413–425. doi:10.1016/S0010-4485(96)00070-X
45. T. Várady, A. Rockwood, P. Salvi. Transfinite surface interpolation over irregular n-sided domains. *Computer-Aided Design* 43 (2011) 1330–1340. doi:10.1016/j.cad.2011.08.028
46. J. Vida, R. R. Martin, T. Várady. A survey of blending methods that use parametric surfaces. *Computer-Aided Design* 26(5) (1994) 341–365. doi:10.1016/0010-4485(94)90023-X
47. X. Ye, H. Nowacki. Ensuring compatibility of G²-continuous surface patches around a nodepoint. *Computer Aided Geometric Design* 13(9) (1996) 931–949.
