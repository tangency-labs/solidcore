# Port progress

Generated from `plan/units/*.toml` and `crates/goldens/expected/*.jsonl` by `harness/plan/progress.py`; do not edit by hand.

**281 / 282 units ported** (254 green, 27 partial), **405,737 / 408,742 upstream lines (99.3%)**. 0 claimed, 1 blocked, 0 todo.

Plan status is historical port progress, not a current CI result. Golden record counts describe stored expectations, not correctness or runtime passes. Known failures are listed separately in `crates/goldens/expected/known_failing.json`.

Legend: ✅ green (port-plan status) · 🟡 partial (port-plan status; some code deferred) · 🔨 claimed · ▶️ ready (all dependencies done) · ⏳ claimable, dependencies still missing · ⛔ blocked.

## Ready to claim now

_None right now; see the next section._

## Claimable now, dependencies still missing

Anything `todo` can be claimed. These need stubs for, or co-porting of, the listed dependencies (AGENTS.md section 4).


## In progress

_None._

## By module

### geom (29/29)

- [x] ✅ `geom.extent3d` (1122 lines) — green, 3086 golden records
- [x] ✅ `geom.vector2d` (333 lines) — green, 4060 golden records
- [x] ✅ `geom.attribute` (920 lines) — green, 2029 golden records
- [x] ✅ `geom.tol` (932 lines) — green, 3616 golden records
- [x] ✅ `geom.tree` (1264 lines) — green, 639 golden records
- [x] ✅ `geom.assert_array` (373 lines) — green, 3735 golden records
- [x] ✅ `geom.context` (355 lines) — green, 2013 golden records
- [x] ✅ `geom.extent1d` (974 lines) — green, 3091 golden records
- [x] ✅ `geom.matrix` (816 lines) — green, 4658 golden records
- [x] ✅ `geom.periodic_extent1d` (406 lines) — green, 4054 golden records
- [x] ✅ `geom.vector3d` (1145 lines) — green, 6114 golden records
- [x] ✅ `geom.axis2_placement` (479 lines) — green, 4051 golden records
- [x] ✅ `geom.extent2d` (766 lines) — green, 3489 golden records
- [x] ✅ `geom.geom_utility_p2` (1617 lines) — green, 3282 golden records
- [x] ✅ `geom.hash_table` (593 lines) — green, 621 golden records
- [x] ✅ `geom.local_solve1d` (1294 lines) — green, 3318 golden records
- [x] ✅ `geom.local_solve_nd` (860 lines) — green, 3043 golden records
- [x] ✅ `geom.object` (679 lines) — green, 624 golden records
- [x] ✅ `geom.polynomial` (386 lines) — green, 6038 golden records
- [x] ✅ `geom.geom_utility_p1` (2178 lines) — green, 4185 golden records
- [x] ✅ `geom.grid` (876 lines) — green, 621 golden records
- [x] ✅ `geom.polar_box` (979 lines) — green, 3033 golden records
- [x] ✅ `geom.solution_array` (1319 lines) — green, 810 golden records
- [x] ✅ `geom.global_solver` (1775 lines) — green, 628 golden records
- [x] ✅ `geom.pseudo_box` (1508 lines) — green, 3012 golden records
- [x] ✅ `geom.tools` (326 lines) — green, 915 golden records
- [x] ✅ `geom.transform` (681 lines) — green, 3042 golden records
- [x] ✅ `geom.a_object` (951 lines) — green, 2040 golden records
- [x] ✅ `geom.point_set` (694 lines) — green, 2032 golden records

### nurbs (49/49)

- [x] ✅ `nlib.knots` (768 lines) — green, 2959 golden records
- [x] ✅ `nlib.math` (2971 lines) — green, 3086 golden records
- [x] ✅ `nlib.geom` (1244 lines) — green, 3176 golden records
- [x] ✅ `nlib.basis_basic` (786 lines) — green, 2605 golden records
- [x] ✅ `nlib.srf_basic` (1088 lines) — green, 3220 golden records
- [x] ✅ `nlib.basis_adv` (361 lines) — green, 1751 golden records
- [x] ✅ `nlib.bezier` (1154 lines) — green, 2223 golden records
- [x] ✅ `nlib.funcs_basic` (538 lines) — green, 3130 golden records
- [x] ✅ `nlib.crv_basic` (1542 lines) — green, 3417 golden records
- [x] ✅ `nlib.crv_fit_p1` (1707 lines) — green, 2089 golden records
- [x] ✅ `nlib.crv_tool_p1` (1380 lines) — green, 2143 golden records
- [x] ✅ `nlib.iges` (1055 lines) — green, 2029 golden records
- [x] ✅ `nlib.srf_tool_p1` (1678 lines) — green, 3133 golden records
- [x] ✅ `nlib.crv_fit_p2` (1586 lines) — green, 2109 golden records
- [x] ✅ `nlib.crv_nurb` (1854 lines) — green, 2886 golden records
- [x] ✅ `nlib.funcs_adv` (1676 lines) — green, 2080 golden records
- [x] ✅ `nlib.srf_fit` (1359 lines) — green, 1556 golden records
- [x] ✅ `nlib.crv_fit_p3` (1234 lines) — green, 2304 golden records
- [x] ✅ `nlib.crv_geom` (729 lines) — green, 1725 golden records
- [x] ✅ `nlib.crv_tool_p2` (1748 lines) — green, 2126 golden records
- [x] ✅ `nurbs.nurbs_crv` (2200 lines) — green, 4113 golden records
- [x] ✅ `nlib.crv_adv_p1` (2031 lines) — green, 2670 golden records
- [x] ✅ `nlib.crv_common` (803 lines) — green, 2140 golden records
- [x] ✅ `nlib.spiral` (392 lines) — green, 2032 golden records
- [x] ✅ `nlib.crv_adv_p2` (1722 lines) — green, 2820 golden records
- [x] ✅ `nlib.srf_common` (661 lines) — green, 1550 golden records
- [x] ✅ `nurbs.b_spline_curve_p1` (1894 lines) — green, 2216 golden records
- [x] ✅ `nlib.crv_approx` (908 lines) — green, 1646 golden records
- [x] ✅ `nlib.crv_offset` (346 lines) — green, 1555 golden records
- [x] ✅ `nlib.srf_nurb` (2518 lines) — green, 2711 golden records
- [x] ✅ `nlib.srf_tool_p2` (1455 lines) — green, 3102 golden records
- [x] ✅ `nurbs.b_spline_curve_n_lib` (2794 lines) — green, 2436 golden records
- [x] ✅ `nurbs.nurbs_srf` (1188 lines) — green, 3092 golden records
- [x] ✅ `nlib.srf_adv_p2` (2097 lines) — green, 1520 golden records
- [x] ✅ `nlib.srf_adv_p3` (2233 lines) — green, 3275 golden records
- [x] ✅ `nlib.srf_tool_p3` (1380 lines) — green, 2524 golden records
- [x] ✅ `nurbs.b_spline_curve_p2` (1930 lines) — green, 2335 golden records
- [x] ✅ `nurbs.b_spline_surface_p1` (1848 lines) — green, 2126 golden records
- [x] ✅ `nlib.srf_tool_p4` (1005 lines) — green, 2681 golden records
- [x] ✅ `nurbs.b_spline_curve_p3` (1256 lines) — green, 2460 golden records
- [x] ✅ `nlib.srf_adv_p4` (2228 lines) — green, 3164 golden records
- [x] ✅ `nlib.srf_adv_p1` (1188 lines) — green, 2758 golden records
- [x] ✅ `nlib.srf_adv_p5` (1970 lines) — green, 3086 golden records
- [x] ✅ `nlib.srf_adv_p6` (2214 lines) — green, 1020 golden records
- [x] ✅ `nlib.srf_offset` (613 lines) — green, 438 golden records
- [x] ✅ `nurbs.b_spline_surface_n_lib_p1` (2033 lines) — green, 701 golden records
- [x] ✅ `nurbs.b_spline_surface_n_lib_p2` (1933 lines) — green, 774 golden records
- [x] ✅ `nurbs.b_spline_surface_p2` (1343 lines) — green, 815 golden records
- [x] ✅ `nlib.crv_fit_p4` (520 lines) — green, 932 golden records

### curve (25/25)

- [x] ✅ `curve.hermite_curve` (383 lines) — green, 2731 golden records
- [x] ✅ `curve.curve_p2` (2302 lines) — green, 3029 golden records
- [x] ✅ `curve.conic` (859 lines) — green, 3240 golden records
- [x] ✅ `curve.crv_in_volume` (949 lines) — green, 2777 golden records
- [x] ✅ `curve.curve_cache` (1484 lines) — green, 642 golden records
- [x] ✅ `curve.hyperbola` (468 lines) — green, 2062 golden records
- [x] ✅ `curve.iso_curve` (310 lines) — green, 2874 golden records
- [x] ✅ `curve.line` (1884 lines) — green, 2105 golden records
- [x] ✅ `curve.offset_curve` (492 lines) — green, 3123 golden records
- [x] ✅ `curve.parabola` (441 lines) — green, 3044 golden records
- [x] ✅ `curve.projected_curve` (569 lines) — green, 3047 golden records
- [x] ✅ `curve.tangent_field` (319 lines) — green, 1355 golden records
- [x] ✅ `curve.ellipse_p2` (1916 lines) — green, 2107 golden records
- [x] ✅ `curve.curve_p4` (2381 lines) — green, 3260 golden records
- [ ] 🟡 `curve.curve_p3` (2072 lines) — partial, 3085 golden records
- [x] ✅ `curve.curve_p1` (2097 lines) — green, 3096 golden records
- [x] ✅ `curve.curve_p6` (2376 lines) — green, 3138 golden records
- [x] ✅ `curve.crv_on_surf` (1728 lines) — green, 2407 golden records
- [x] ✅ `curve.curve_p7` (1666 lines) — green, 3067 golden records
- [x] ✅ `curve.ellipse_p1` (1932 lines) — green, 2527 golden records
- [x] ✅ `curve.curve_p8` (1670 lines) — green, 3069 golden records
- [x] ✅ `curve.composite_curve_p1` (1700 lines) — green, 3096 golden records
- [x] ✅ `curve.composite_curve_p2` (1345 lines) — green, 676 golden records
- [ ] 🟡 `curve.curve_p5` (2299 lines) — partial, 3140 golden records
- [x] ✅ `curve.curve_p9` (1783 lines) — green, 3126 golden records

### surface (42/42)

- [x] ✅ `surface.volume_p1` (2014 lines) — green, 2886 golden records
- [x] ✅ `surface.b_spline_volume` (1226 lines) — green, 2333 golden records
- [x] ✅ `surface.twist_volume` (649 lines) — green, 1572 golden records
- [ ] 🟡 `surface.surface_p2` (2398 lines) — partial, 17593 golden records
- [x] ✅ `surface.curve_bounded_surface` (324 lines) — green, 537 golden records
- [x] ✅ `surface.surface_cache_p1` (1944 lines) — green, 872 golden records
- [ ] 🟡 `surface.plane_p1` (1934 lines) — partial, 2698 golden records
- [x] ✅ `surface.surf_of_revolution_p1` (1522 lines) — green, 895 golden records
- [x] ✅ `surface.surface_p9` (1252 lines) — green, 833 golden records
- [x] ✅ `surface.trim_srf_cache` (966 lines) — green, 779 golden records
- [x] ✅ `surface.b_spline_volume_n_lib` (310 lines) — green, 2420 golden records
- [x] ✅ `surface.bend_volume` (945 lines) — green, 2347 golden records
- [x] ✅ `surface.cone_p1` (1835 lines) — green, 20736 golden records
- [x] ✅ `surface.surf_of_extrusion_p1` (1865 lines) — green, 2071 golden records
- [x] ✅ `surface.surf_of_revolution_p3` (1620 lines) — green, 26731 golden records
- [x] ✅ `surface.surface_cache_p2` (1252 lines) — green, 913 golden records
- [ ] 🟡 `surface.torus` (1180 lines) — partial, 758 golden records
- [x] ✅ `surface.unbend_volume` (920 lines) — green, 1874 golden records
- [x] ✅ `surface.cone_p2` (1687 lines) — green, 1713 golden records
- [x] ✅ `surface.surf_of_revolution_p2` (1369 lines) — green, 1167 golden records
- [x] ✅ `surface.surface_cache_p3` (1703 lines) — green, 600 golden records
- [x] ✅ `surface.volume_p2` (1490 lines) — green, 11220 golden records
- [x] ✅ `surface.surface_p7` (2407 lines) — green, 3741 golden records
- [x] ✅ `surface.plane_p2` (1627 lines) — green, 2228 golden records
- [x] ✅ `surface.surface_p4` (2270 lines) — green, 3377 golden records
- [x] ✅ `surface.sphere` (2861 lines) — green, 816 golden records
- [x] ✅ `surface.surf_of_extrusion_p2` (1420 lines) — green, 1528 golden records
- [x] ✅ `surface.surface_p5` (2321 lines) — green, 3118 golden records
- [x] ✅ `surface.cone_p3` (1940 lines) — green, 1490 golden records
- [x] ✅ `surface.surface_p1` (2458 lines) — green, 1000 golden records
- [x] ✅ `surface.surface_p6` (2150 lines) — green, 2104 golden records
- [x] ✅ `surface.srf_in_volume` (779 lines) — green, 4567 golden records
- [x] ✅ `surface.surface_p10` (2172 lines) — green, 867 golden records
- [x] ✅ `surface.surface_p3` (1794 lines) — green, 2507 golden records
- [x] ✅ `surface.offset_surface` (1787 lines) — green, 3079 golden records
- [x] ✅ `surface.surface_assert_valid` (1100 lines) — green, 2902 golden records
- [x] ✅ `surface.surface_p11` (1789 lines) — green, 814 golden records
- [x] ✅ `surface.surface_p8` (2374 lines) — green, 2108 golden records
- [x] ✅ `surface.surface_approximate` (640 lines) — green, 437 golden records
- [x] ✅ `surface.surface_point_solve_p1` (1771 lines) — green, 1035 golden records
- [x] ✅ `surface.surface_point_solve_p2` (1882 lines) — green, 1366 golden records
- [x] ✅ `surface.surface_point_solve_p3` (1282 lines) — green, 1365 golden records

### topology (40/40)

- [x] ✅ `topology.database_io` (366 lines) — green, 2066 golden records
- [x] ✅ `topology.topology` (393 lines) — green, 658 golden records
- [x] ✅ `topology.brep_data` (1466 lines) — green, 2379 golden records
- [x] ✅ `topology.proto_topology_p1` (2773 lines) — green, 1922 golden records
- [x] ✅ `topology.edgeuse_p1` (1158 lines) — green, 5586 golden records
- [x] ✅ `topology.brep_validation` (1500 lines) — green, 112 golden records
- [ ] 🟡 `topology.face_p2` (1884 lines) — partial, 2692 golden records
- [ ] 🟡 `topology.assert_valid` (2419 lines) — partial, 965 golden records
- [x] ✅ `topology.edgeuse_p2` (1013 lines) — green, 4795 golden records
- [x] ✅ `topology.shell` (1085 lines) — green, 944 golden records
- [x] ✅ `topology.vertexuse` (1475 lines) — green, 2686 golden records
- [x] ✅ `topology.vertex` (1791 lines) — green, 5954 golden records
- [x] ✅ `topology.loopuse` (958 lines) — green, 9013 golden records
- [x] ✅ `topology.brep_p8` (1768 lines) — green, 355 golden records
- [x] ✅ `topology.edge` (2306 lines) — green, 2949 golden records
- [x] ✅ `topology.edgeuse_p3` (2273 lines) — green, 5611 golden records
- [x] ✅ `topology.trimming_tools` (2885 lines) — green, 3305 golden records
- [x] ✅ `topology.face_p4` (2273 lines) — green, 5930 golden records
- [x] ✅ `topology.brep_p10` (1890 lines) — green, 1436 golden records
- [x] ✅ `topology.face_p3` (1878 lines) — green, 2407 golden records
- [x] ✅ `topology.topology_traverser` (784 lines) — green, 909 golden records
- [x] ✅ `topology.brep_p7` (2503 lines) — green, 462 golden records
- [x] ✅ `topology.face_p1` (2230 lines) — green, 2569 golden records
- [x] ✅ `topology.brep_p1` (2511 lines) — green, 2636 golden records
- [x] ✅ `topology.brep_p6` (1479 lines) — green, 595 golden records
- [x] ✅ `topology.brep_cache` (389 lines) — green, 663 golden records
- [x] ✅ `topology.brep_p4` (2113 lines) — green, 327 golden records
- [x] ✅ `topology.brep_p9` (2484 lines) — green, 710 golden records
- [x] ✅ `topology.brep_p2` (2578 lines) — green, 1720 golden records
- [x] ✅ `topology.brep_p5` (2193 lines) — green, 797 golden records
- [x] ✅ `topology.brep_mirror` (330 lines) — green, 341 golden records
- [ ] 🟡 `topology.brep_p3` (2527 lines) — partial, 1358 golden records
- [ ] 🟡 `topology.brep_planar_section` (520 lines) — partial, 167 golden records
- [x] ✅ `topology.proto_topology_p2` (2115 lines) — green, 1767 golden records
- [ ] 🟡 `topology.api_files` (132 lines) — partial, 339 golden records
- [ ] 🟡 `topology.brep_cutting` (960 lines) — partial, 108 golden records
- [ ] 🟡 `topology.brep_local_operation` (330 lines) — partial, 161 golden records
- [x] ✅ `topology.proto_topology_p3` (1581 lines) — green, 3013 golden records
- [x] ✅ `topology.brep_constructor` (850 lines) — green, 185 golden records
- [ ] 🟡 `topology.api_projection` (330 lines) — partial, 48 golden records

### intersect (26/26)

- [x] ✅ `intersect.surface_tracer` (1900 lines) — green, 427 golden records
- [x] ✅ `intersect.point_class` (909 lines) — green, 2520 golden records
- [x] ✅ `intersect.curve_class_p1` (2209 lines) — green, 2131 golden records
- [ ] 🟡 `intersect.surface_drop_lib_p1` (1478 lines) — partial, 2926 golden records
- [x] ✅ `intersect.surface_silhouette` (979 lines) — green, 709 golden records
- [x] ✅ `intersect.curve_class_p7` (2248 lines) — green, 1673 golden records
- [x] ✅ `intersect.ray_tracer` (1246 lines) — green, 433 golden records
- [x] ✅ `intersect.surface_drop_curve` (2014 lines) — green, 547 golden records
- [x] ✅ `intersect.surface_drop_lib_p2` (1080 lines) — green, 1531 golden records
- [x] ✅ `intersect.surface_intersector_p1` (1816 lines) — green, 1510 golden records
- [x] ✅ `intersect.topology_solver_p1` (1838 lines) — green, 3592 golden records
- [x] ✅ `intersect.adv_surface_intersector` (1646 lines) — green, 1949 golden records
- [x] ✅ `intersect.curve_class_p2` (2057 lines) — green, 1671 golden records
- [x] ✅ `intersect.surface_drop_lib_p3` (2070 lines) — green, 1532 golden records
- [ ] 🟡 `intersect.surface_intersector_p2` (1677 lines) — partial, 1686 golden records
- [x] ✅ `intersect.curve_class_p3` (1986 lines) — green, 3708 golden records
- [x] ✅ `intersect.curve_class_p4` (336 lines) — green, 5410 golden records
- [x] ✅ `intersect.surface_intersector_p3` (1932 lines) — green, 1527 golden records
- [x] ✅ `intersect.curve_class_p6` (2213 lines) — green, 4470 golden records
- [x] ✅ `intersect.surface_intersector_p4` (1570 lines) — green, 1542 golden records
- [x] ✅ `intersect.topology_solver_p2` (1391 lines) — green, 2575 golden records
- [x] ✅ `intersect.curve_class_p5` (3546 lines) — green, 1411 golden records
- [x] ✅ `intersect.surface_intersector_p5` (622 lines) — green, 1731 golden records
- [x] ✅ `intersect.curve_class_p8` (1900 lines) — green, 1127 golden records
- [x] ✅ `intersect.curve_class_p9` (1220 lines) — green, 457 golden records
- [ ] 🟡 `intersect.api_intersectors` (200 lines) — partial, 477 golden records

### primitive (8/8)

- [x] ✅ `primitive.api_curves` (649 lines) — green, 2046 golden records
- [x] ✅ `primitive.api_curve_ops` (190 lines) — green, 638 golden records
- [x] ✅ `primitive.primitive_creation_p1` (2367 lines) — green, 773 golden records
- [x] ✅ `primitive.api_surface_ops` (200 lines) — green, 570 golden records
- [ ] 🟡 `primitive.api_trim_ops` (200 lines) — partial, 211 golden records
- [x] ✅ `primitive.primitive_creation_p3` (1833 lines) — green, 330 golden records
- [x] ✅ `primitive.primitive_creation_p2` (2258 lines) — green, 374 golden records
- [x] ✅ `primitive.api_primitives` (910 lines) — green, 2245 golden records

### query (12/12)

- [x] ✅ `query.edge` (1025 lines) — green, 1368 golden records
- [x] ✅ `query.loop` (830 lines) — green, 3120 golden records
- [x] ✅ `query.face_p2` (1161 lines) — green, 1401 golden records
- [x] ✅ `query.face_p3` (1625 lines) — green, 1199 golden records
- [x] ✅ `query.poly` (2119 lines) — green, 966 golden records
- [x] ✅ `query.face_p1` (465 lines) — green, 2266 golden records
- [x] ✅ `query.face_p4` (1506 lines) — green, 503 golden records
- [x] ✅ `query.poly_copy` (700 lines) — green, 524 golden records
- [x] ✅ `query.poly_properties` (360 lines) — green, 641 golden records
- [x] ✅ `query.tess` (838 lines) — green, 676 golden records
- [x] ✅ `query.brep` (517 lines) — green, 1172 golden records
- [x] ✅ `query.api_queries` (843 lines) — green, 1399 golden records

### mass (1/1)

- [ ] 🟡 `mass.integrator` (1066 lines) — partial, 1628 golden records

### sweep (5/5)

- [x] ✅ `sweep.topology_sweep` (2456 lines) — green, 703 golden records
- [x] ✅ `sweep.rotational_sweep_geometry` (1542 lines) — green, 823 golden records
- [x] ✅ `sweep.trans_sweep_geometry` (817 lines) — green, 788 golden records
- [ ] 🟡 `sweep.api_brep` (529 lines) — partial, 696 golden records
- [x] ✅ `sweep.api_surfaces` (390 lines) — green, 885 golden records

### heal (7/7)

- [x] ✅ `heal.gap` (1475 lines) — green, 749 golden records
- [x] ✅ `heal.face_props` (1482 lines) — green, 1873 golden records
- [x] ✅ `heal.heal_data_p2` (2164 lines) — green, 741 golden records
- [x] ✅ `heal.heal_data_p1` (922 lines) — green, 495 golden records
- [x] ✅ `heal.stitch_p1` (2037 lines) — green, 674 golden records
- [x] ✅ `heal.heal_data_p3` (1939 lines) — green, 702 golden records
- [x] ✅ `heal.stitch_p2` (1906 lines) — green, 841 golden records

### offset (3/3)

- [ ] 🟡 `offset.offset_geometry_creation` (2661 lines) — partial, 926 golden records
- [ ] 🟡 `offset.offset_executive_p1` (2218 lines) — partial, 965 golden records
- [ ] 🟡 `offset.offset_executive_p2` (1847 lines) — partial, 568 golden records

### boolean (9/9)

- [x] ✅ `boolean.topology_intersector_p1` (2491 lines) — green, 845 golden records
- [x] ✅ `boolean.topology_intersector_p2` (2195 lines) — green, 913 golden records
- [x] ✅ `boolean.topology_intersector_p3` (1849 lines) — green, 271 golden records
- [x] ✅ `boolean.merge_p1` (993 lines) — green, 251 golden records
- [x] ✅ `boolean.merge_p2` (1967 lines) — green, 45 golden records
- [x] ✅ `boolean.merge_p5` (475 lines) — green, 140 golden records
- [x] ✅ `boolean.merge_p4` (2596 lines) — green, 311 golden records
- [x] ✅ `boolean.merge_p3` (1152 lines) — green, 438 golden records
- [ ] 🟡 `boolean.api_brep` (410 lines) — partial, 1126 golden records

### blend (24/25)

- [x] ✅ `blend.fillet_custom_solver` (1236 lines) — green, 451 golden records
- [x] ✅ `blend.identify_fillet` (1008 lines) — green, 510 golden records
- [x] ✅ `blend.fillet_geom_p1` (2113 lines) — green, 2199 golden records
- [x] ✅ `blend.fillet_solver_p1` (1861 lines) — green, 351 golden records
- [x] ✅ `blend.hcr` (1875 lines) — green, 481 golden records
- [x] ✅ `blend.fillet_intersector` (1708 lines) — green, 1292 golden records
- [x] ✅ `blend.fillet_solver_p5` (848 lines) — green, 100 golden records
- [x] ✅ `blend.fillet_standard_solver_p1` (2336 lines) — green, 3619 golden records
- [x] ✅ `blend.fillet_solver_p4` (1404 lines) — green, 651 golden records
- [ ] 🟡 `blend.fillet_executive_p3` (1608 lines) — partial, 793 golden records
- [x] ✅ `blend.fillet_executive_p1` (2113 lines) — green, 155 golden records
- [x] ✅ `blend.fillet_corner_p1` (1971 lines) — green, 263 golden records
- [x] ✅ `blend.fillet_solver_p2` (1939 lines) — green, 311 golden records
- [x] ✅ `blend.fillet_solver_p3` (1112 lines) — green, 400 golden records
- [x] ✅ `blend.fillet_geom_p3` (1032 lines) — green, 325 golden records
- [ ] 🟡 `blend.fillet_geom_p4` (2074 lines) — partial, 400 golden records
- [x] ✅ `blend.fillet_geom_p6` (1813 lines) — green, 400 golden records
- [x] ✅ `blend.fillet_standard_solver_p2` (1989 lines) — green, 525 golden records
- [ ] 🟡 `blend.fillet_corner_p3` (1697 lines) — partial, 900 golden records
- [x] ✅ `blend.fillet_executive_p2` (1155 lines) — green, 278 golden records
- [ ] ⛔ `blend.feature_executive` (3005 lines) — blocked (see notes / plan/findings.md)
- [x] ✅ `blend.fillet_geom_p2` (2077 lines) — green, 1151 golden records
- [x] ✅ `blend.fillet_removal` (1757 lines) — green, 170 golden records
- [x] ✅ `blend.fillet_geom_p5` (1952 lines) — green, 300 golden records
- [ ] 🟡 `blend.fillet_corner_p2` (1892 lines) — partial, 547 golden records

### (crate-wide) (1/1)

- [x] ✅ `audit.recursion` (0 lines) — green
