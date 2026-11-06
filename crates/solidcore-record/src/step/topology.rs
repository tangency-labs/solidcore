//! STEP topology (ISO 10303-42 topology and geometric model entities) as plain data: each
//! solid item of a file with its shells, faces, bounds, oriented edges, edges and vertices as the
//! file states them, linked by indices, with their geometry read through [`super::read`] in the
//! units of the item's representation context ([`super::units`]).
//!
//! Each `MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS` and `SHELL_BASED_SURFACE_MODEL` is one
//! [`StepSolid`] (a surface model with closed and open shells is one per closed shell, then one
//! for its open shells: [`items`]), and each face no shell lists is one of its own (a loose face,
//! [`loose_faces`]). A solid that cannot be read is listed with its reason and the rest of the
//! file is still read.
//!
//! The data is the file's: every orientation flag is kept where the file writes it (a shell's,
//! an oriented face's, a face's `same_sense`, a bound's, an oriented edge's, an edge's
//! `same_sense`), a chain of `ORIENTED_FACE`s or `ORIENTED_EDGE`s is folded into one flag (the
//! flags agree or not: no arithmetic), bounds and edges are in the file's order, an edge's
//! curve is the 3D curve under its `SURFACE_CURVE`s and `TRIMMED_CURVE`s as the file has it, and
//! a vertex is its point. Nothing is evaluated, projected or measured: an edge's interval on its
//! curve, a face's parameter window, the seams a closed surface needs, the outer loop of a face
//! that marks none, tolerances raised to measured gaps and the side a void faces are the
//! kernel's when ops builds a body from this data (#2578, slice 4: #2609, #2622, #2635, #2623,
//! and the kernel's build, sew and heal). It is not the record's `Brep`.
//!
//! Input is hostile: every reference is checked, nesting is capped ([`MAX_NESTING`]), counts are
//! capped ([`MAX_FACES`], [`MAX_LOOP`], [`MAX_BOUNDS`], [`MAX_SHELLS`], [`MAX_SOLIDS`]) and the
//! work of a whole file is bounded ([`MAX_FILE_FACES`], a shell counted again for each solid
//! that uses it). Every error is a message, never a panic.
//!
//! Origin: SolidCraft's `crates/step/src/topology.rs` at `a5b71ca2` (MIT OR Apache-2.0), its
//! data mapping only; loose faces as solidcore's `step-import` maps them (`map_loose_faces`,
//! `crates/step-import/src/map.rs`).
//!
//! | Here                                  | SolidCraft (`topology.rs`)                    |
//! |---------------------------------------|-----------------------------------------------|
//! | [`read_solids`], [`read_solid`]       | `read_solids`, `read_solids_in`, `read_solid` |
//! | [`items`], `model_shells`, `item_shells`, `shell_faces` | the same               |
//! | `Topology::{vertex, edge, edge_geometry, oriented_edge, bound, face, shell}` | `Builder::…`, without the geometry |
//! | [`StepSolid::ends`], [`StepShell::closed`], [`StepFace::loop_order`] | `Builder::{ends, closed, loop_order}` |
//! | [`StepLoop::seam_only`], [`reverse`]  | `seam_only`, `reverse`                        |
//! | [`face_entity`]                       | `Builder::face_entity`                        |
//! | [`StepSolid::edge_uses`]              | the uses `finish` lists per edge (`EdgeUseRef`) |
//! | [`loose_faces`]                       | step-import's `map_loose_faces` candidates    |
//!
//! Deviations from SolidCraft, each without a computation in its place:
//! - The output is this data, not a schema-4 `Brep`: no pcurves, intervals, seams, frames,
//!   domains, NURBS forms or canonical order, and no record validation (all of `finish`).
//! - Loops keep the file's direction and orientation flags (was: reversed for a bound's
//!   orientation `.F.` and a face's `same_sense`, read relative to the surface's natural normal);
//!   [`StepBound::face_loop`] gives a bound's loop as the face uses it.
//! - An edge keeps its vertices and `same_sense` as written (was: the schema edge run with or
//!   against its curve by the vertices' parameters, `map_edge`), so a loop's connectivity is
//!   checked on the file's vertices.
//! - A shell's closedness ([`StepShell::closed`], each edge used at least twice among its faces)
//!   is data; a solid whose outer shell or a void is not closed is not refused here (was: "the
//!   outer shell is not closed"): whether sewing closes it is the build's.
//! - The context's distance uncertainty is kept as the file states it ([`StepSolid::tolerance`],
//!   `None` when it states none); the default when there is none and the clamp to SolidCraft's
//!   least and loosest tolerances are the build's choice (was `tolerance`).
//! - A face bounded by vertex loops alone (a whole sphere bounded by one vertex loop) is data
//!   (was refused: "not supported yet: a face without an edge loop"); a face with no bound at
//!   all is refused ("has no bounds"), as ISO 10303-42 gives a face one at least.
//! - Loose faces are new (SolidCraft read no face outside a shell).
//! - Maps and sets are `BTreeMap`/`BTreeSet` (were `HashMap`/`HashSet`); nothing is pre-reserved
//!   from a count in the file. Recursion: `ORIENTED_EDGE`, `ORIENTED_FACE` and edge geometry
//!   chains each carry a depth, refused past [`MAX_NESTING`] (as in SolidCraft).

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use super::geom::{Curve, P3, Surface};
use super::p21::{Entity, Exchange, Param};
use super::read::{R, Reader};
use super::units::{self, Units};

/// Most faces in one solid.
pub const MAX_FACES: usize = 500_000;
/// Most edges in one loop, and in all of one face's loops.
pub const MAX_LOOP: usize = 100_000;
/// Most bounds of one face.
pub const MAX_BOUNDS: usize = 10_000;
/// Most faces read in one file, counting a shell again for each solid that uses it.
pub const MAX_FILE_FACES: usize = 5_000_000;
/// Most solid items (loose faces included) read from one file.
pub const MAX_SOLIDS: usize = 100_000;
/// Most shells of one solid item (the record's own bound, [`crate::brep::MAX_SHELLS`]).
pub const MAX_SHELLS: usize = crate::brep::MAX_SHELLS;
/// Deepest nesting of `ORIENTED_EDGE` in `ORIENTED_EDGE`, `ORIENTED_FACE` in `ORIENTED_FACE`,
/// and of an edge's curve in `SURFACE_CURVE`s and `TRIMMED_CURVE`s.
pub const MAX_NESTING: usize = 4;

/// The STEP entities a solid item is (`FACETED_BREP` is listed so that it is refused by name).
const SOLID_ITEMS: [&str; 4] = [
    "MANIFOLD_SOLID_BREP",
    "BREP_WITH_VOIDS",
    "SHELL_BASED_SURFACE_MODEL",
    "FACETED_BREP",
];

/// What a [`StepSolid`] is in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolidKind {
    /// A `MANIFOLD_SOLID_BREP`: one closed shell.
    ManifoldSolid,
    /// A `BREP_WITH_VOIDS`: an outer shell and its voids.
    BrepWithVoids,
    /// A `SHELL_BASED_SURFACE_MODEL` (only its open shells when it has closed ones too).
    SurfaceModel,
    /// A closed shell of a surface model that has open shells too, read as a solid of its own.
    ModelClosedShell,
    /// A face no shell lists, read as a sheet of its own.
    LooseFace,
}

/// A solid item of the file, with its topology as plain data. Indices refer to this solid's
/// own `faces`, `edges` and `vertices`.
#[derive(Clone, Debug, PartialEq)]
pub struct StepSolid {
    /// The STEP item: the solid entity, the closed shell (for [`SolidKind::ModelClosedShell`]),
    /// or the face (for [`SolidKind::LooseFace`]).
    pub item: u64,
    /// What the item is.
    pub kind: SolidKind,
    /// The item's own name, else its product's.
    pub name: String,
    /// The name of the product whose shape holds the item, when the file says.
    pub product: Option<String>,
    /// The units its geometry was read in (millimetres and radians when its context names none).
    pub units: Units,
    /// The distance uncertainty its context states, in millimetres, if any.
    pub tolerance: Option<f64>,
    /// The outer shell, or the surface model's shells, in the file's order.
    pub shells: Vec<StepShell>,
    /// A `BREP_WITH_VOIDS`' voids, in the file's order.
    pub voids: Vec<StepShell>,
    /// Every face, in the order the shells first use them.
    pub faces: Vec<StepFace>,
    /// Every edge, in the order the faces' loops first use them.
    pub edges: Vec<StepEdge>,
    /// Every vertex, in the order edges and vertex loops first use them.
    pub vertices: Vec<StepVertex>,
}

/// A shell as the solid uses it.
#[derive(Clone, Debug, PartialEq)]
pub struct StepShell {
    /// The shell entity the solid lists (the `ORIENTED_…_SHELL` when it is one); `None` for the
    /// sheet of a loose face.
    pub entity: Option<u64>,
    /// The file calls it closed (`CLOSED_SHELL`, `ORIENTED_CLOSED_SHELL`).
    pub declared_closed: bool,
    /// An `ORIENTED_…_SHELL`'s orientation (`.F.`: the shell uses its faces turned over);
    /// `true` for a plain shell.
    pub orientation: bool,
    /// Its faces, in the file's order.
    pub faces: Vec<StepOrientedFace>,
    /// Every edge of its faces' edge loops is used at least twice among them.
    pub closed: bool,
}

/// A face as a shell lists it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StepOrientedFace {
    /// The entity the shell lists (an `ORIENTED_FACE` or the face itself).
    pub entity: u64,
    /// Index into [`StepSolid::faces`].
    pub face: u32,
    /// The orientations of the `ORIENTED_FACE`s on the way to the face, folded (`true` when they
    /// turn it over an even number of times, or there are none).
    pub orientation: bool,
}

/// What a [`StepFace`] is in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FaceKind {
    /// `ADVANCED_FACE`.
    Advanced,
    /// `FACE_SURFACE`.
    FaceSurface,
}

/// An `ADVANCED_FACE` or `FACE_SURFACE`.
#[derive(Clone, Debug, PartialEq)]
pub struct StepFace {
    /// The face entity.
    pub entity: u64,
    /// Which entity it is.
    pub kind: FaceKind,
    /// The `face_geometry` entity.
    pub surface_entity: u64,
    /// Its surface.
    pub surface: Arc<Surface>,
    /// `same_sense`: the face's normal is the surface's.
    pub same_sense: bool,
    /// Its bounds, in the file's order.
    pub bounds: Vec<StepBound>,
}

/// A `FACE_BOUND` or `FACE_OUTER_BOUND`.
#[derive(Clone, Debug, PartialEq)]
pub struct StepBound {
    /// The bound entity.
    pub entity: u64,
    /// A `FACE_OUTER_BOUND`.
    pub outer: bool,
    /// `orientation`: the face uses the loop as written (`.F.`: reversed).
    pub orientation: bool,
    /// Its loop.
    pub r#loop: StepLoop,
}

/// A bound's loop.
#[derive(Clone, Debug, PartialEq)]
pub enum StepLoop {
    /// An `EDGE_LOOP`: its oriented edges in the file's order, connected and closed.
    Edges {
        /// The loop entity.
        entity: u64,
        /// Its oriented edges.
        edges: Vec<StepOrientedEdge>,
    },
    /// A `VERTEX_LOOP`.
    Vertex {
        /// The loop entity.
        entity: u64,
        /// Index into [`StepSolid::vertices`].
        vertex: u32,
    },
}

/// An `ORIENTED_EDGE` (or a bare `EDGE_CURVE`) as a loop lists it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StepOrientedEdge {
    /// The entity the loop lists.
    pub entity: u64,
    /// Index into [`StepSolid::edges`].
    pub edge: u32,
    /// The loop runs from the edge's start to its end (the orientations of nested
    /// `ORIENTED_EDGE`s folded; `true` for a bare edge).
    pub orientation: bool,
}

/// An `EDGE_CURVE`.
#[derive(Clone, Debug, PartialEq)]
pub struct StepEdge {
    /// The edge entity.
    pub entity: u64,
    /// `edge_start`: index into [`StepSolid::vertices`].
    pub start: u32,
    /// `edge_end`: index into [`StepSolid::vertices`].
    pub end: u32,
    /// The entity its curve was read from (under `edge_geometry`'s surface and trimmed curves).
    pub curve_entity: u64,
    /// The 3D curve: `edge_geometry` without its `SURFACE_CURVE`s and `TRIMMED_CURVE`s (the
    /// vertices bound it).
    pub curve: Arc<Curve>,
    /// `same_sense`: the edge runs along its curve's parameter.
    pub same_sense: bool,
}

/// A `VERTEX_POINT`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepVertex {
    /// The vertex entity.
    pub entity: u64,
    /// Its point, in millimetres.
    pub point: P3,
}

/// Where a loop uses an edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StepEdgeUse {
    /// Index into [`StepSolid::faces`].
    pub face: u32,
    /// Index into the face's bounds.
    pub bound: u32,
    /// Index into the bound's loop, as written.
    pub position: u32,
    /// The oriented edge's orientation.
    pub orientation: bool,
}

/// A solid item that could not be read, and why.
#[derive(Clone, Debug, PartialEq)]
pub struct StepFailure {
    /// The STEP item.
    pub item: u64,
    /// Why.
    pub message: String,
}

/// Every solid item of a file: those read and those refused.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StepSolids {
    /// The solids read, solid items in id order, then loose faces in id order.
    pub solids: Vec<StepSolid>,
    /// The items refused, in the same order.
    pub failures: Vec<StepFailure>,
}

impl StepSolid {
    /// The vertices an oriented edge runs from and to.
    pub fn ends(&self, oe: StepOrientedEdge) -> R<[u32; 2]> {
        ends(&self.edges, oe)
    }

    /// Each edge's uses (by index into [`Self::edges`]), walking faces, their bounds and their
    /// loops in order.
    pub fn edge_uses(&self) -> R<Vec<Vec<StepEdgeUse>>> {
        let mut uses: Vec<Vec<StepEdgeUse>> = Vec::new();
        uses.resize_with(self.edges.len(), Vec::new);
        for (f, face) in self.faces.iter().enumerate() {
            for (b, bound) in face.bounds.iter().enumerate() {
                let StepLoop::Edges { edges, .. } = &bound.r#loop else {
                    continue;
                };
                for (p, oe) in edges.iter().enumerate() {
                    uses.get_mut(oe.edge as usize)
                        .ok_or("edge reference")?
                        .push(StepEdgeUse {
                            face: index(f)?,
                            bound: index(b)?,
                            position: index(p)?,
                            orientation: oe.orientation,
                        });
                }
            }
        }
        Ok(uses)
    }
}

impl StepFace {
    /// The face's bounds, outer first: the first `FACE_OUTER_BOUND`, else the first bound, then
    /// the rest in the file's order. Which loop is outer where the file marks none is the
    /// kernel's to find.
    pub fn loop_order(&self) -> Vec<usize> {
        let outer = self.bounds.iter().position(|b| b.outer).unwrap_or(0);
        std::iter::once(outer)
            .chain((0..self.bounds.len()).filter(|l| *l != outer))
            .filter(|l| *l < self.bounds.len())
            .collect()
    }
}

impl StepBound {
    /// The bound's oriented edges as the face uses them: as written, or reversed (each edge
    /// the other way) when the bound's orientation is `.F.`. A vertex loop has none.
    pub fn face_loop(&self) -> Vec<StepOrientedEdge> {
        let StepLoop::Edges { edges, .. } = &self.r#loop else {
            return Vec::new();
        };
        let mut uses = edges.clone();
        if !self.orientation {
            reverse(&mut uses);
        }
        uses
    }
}

impl StepLoop {
    /// Every edge of the loop is used in it exactly twice, once each way: a loop of seam uses
    /// alone (a ball bounded by its meridian). `false` for a vertex loop.
    pub fn seam_only(&self) -> bool {
        let StepLoop::Edges { edges, .. } = self else {
            return false;
        };
        let mut ways: BTreeMap<u32, [u32; 2]> = BTreeMap::new();
        for oe in edges {
            if let Some(n) = ways
                .entry(oe.edge)
                .or_default()
                .get_mut(usize::from(oe.orientation))
            {
                *n = n.saturating_add(1);
            }
        }
        !edges.is_empty() && ways.values().all(|n| *n == [1, 1])
    }
}

/// The loop run the other way: the edges in reverse order, each the other way.
pub fn reverse(uses: &mut [StepOrientedEdge]) {
    uses.reverse();
    for u in uses {
        u.orientation = !u.orientation;
    }
}

fn ends(edges: &[StepEdge], oe: StepOrientedEdge) -> R<[u32; 2]> {
    let e = edges.get(oe.edge as usize).ok_or("edge reference")?;
    Ok(if oe.orientation {
        [e.start, e.end]
    } else {
        [e.end, e.start]
    })
}

fn entity(ex: &Exchange, id: u64) -> R<&Entity> {
    ex.get(id).ok_or_else(|| format!("missing entity #{id}"))
}

fn rid(p: Option<&Param>, what: &str) -> R<u64> {
    p.and_then(Param::as_ref_id)
        .ok_or_else(|| format!("{what}: expected a reference"))
}

fn refs(p: Option<&Param>, what: &str) -> R<Vec<u64>> {
    let list = p
        .and_then(Param::as_list)
        .ok_or_else(|| format!("{what}: expected a list"))?;
    list.iter()
        .map(|x| {
            x.as_ref_id()
                .ok_or_else(|| format!("{what}: expected references"))
        })
        .collect()
}

fn flag(p: Option<&Param>, what: &str) -> R<bool> {
    p.and_then(Param::as_bool)
        .ok_or_else(|| format!("{what}: expected a boolean"))
}

fn unsupported(what: impl std::fmt::Display) -> String {
    format!("not supported yet: {what}")
}

fn index(i: usize) -> R<u32> {
    u32::try_from(i).map_err(|_| "too many items".to_string())
}

/// Items and context of a representation entity (plain or complex).
fn representation(e: &Entity) -> Option<(Vec<u64>, Option<u64>)> {
    let p = if e.is_complex() {
        e.record("REPRESENTATION").or_else(|| {
            e.records
                .iter()
                .find(|r| r.name.ends_with("REPRESENTATION") && r.params.len() >= 3)
                .map(|r| r.params.as_slice())
        })?
    } else if e.name().ends_with("REPRESENTATION") {
        e.params()
    } else {
        return None;
    };
    let items = p
        .get(1)?
        .as_list()?
        .iter()
        .filter_map(Param::as_ref_id)
        .collect();
    Some((items, p.get(2).and_then(Param::as_ref_id)))
}

/// Where each item sits: its representation's context, and the product that representation
/// (or one related to it) is the shape of.
struct Index {
    context: BTreeMap<u64, u64>,
    product: BTreeMap<u64, String>,
}

impl Index {
    fn new(ex: &Exchange) -> Index {
        let mut context = BTreeMap::new();
        let mut rep_of: BTreeMap<u64, u64> = BTreeMap::new();
        // representation → PRODUCT_DEFINITION_SHAPE
        let mut shape_of: BTreeMap<u64, u64> = BTreeMap::new();
        let mut links: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (&id, e) in &ex.entities {
            if let Some((items, ctx)) = representation(e) {
                for item in items {
                    if let Some(c) = ctx {
                        context.entry(item).or_insert(c);
                    }
                    rep_of.entry(item).or_insert(id);
                }
            }
            let p = e.params();
            let rf = |i: usize| p.get(i).and_then(Param::as_ref_id);
            match e.name() {
                "SHAPE_DEFINITION_REPRESENTATION" => {
                    if let (Some(d), Some(r)) = (rf(0), rf(1)) {
                        shape_of.entry(r).or_insert(d);
                    }
                }
                "SHAPE_REPRESENTATION_RELATIONSHIP" | "REPRESENTATION_RELATIONSHIP"
                    if !e.is_complex() =>
                {
                    if let (Some(a), Some(b)) = (rf(2), rf(3)) {
                        links.entry(a).or_default().push(b);
                        links.entry(b).or_default().push(a);
                    }
                }
                _ => {}
            }
        }
        let name_of = |pds: u64| -> Option<String> {
            let pd = ex.get(pds)?.params().get(2)?.as_ref_id()?;
            let pdf = ex.get(pd)?.params().get(2)?.as_ref_id()?;
            let product = ex.get(pdf)?.params().get(2)?.as_ref_id()?;
            let p = ex.get(product)?;
            (p.name() == "PRODUCT")
                .then(|| p.params().get(1).or_else(|| p.params().first()))
                .flatten()?
                .as_str()
                .map(str::to_string)
        };
        let mut product = BTreeMap::new();
        for (&item, &rep) in &rep_of {
            if !ex
                .get(item)
                .is_some_and(|e| SOLID_ITEMS.contains(&e.name()))
            {
                continue;
            }
            let direct = shape_of.get(&rep).copied();
            let linked = || {
                links
                    .get(&rep)
                    .into_iter()
                    .flatten()
                    .find_map(|r| shape_of.get(r).copied())
            };
            if let Some(name) = direct.or_else(linked).and_then(name_of) {
                product.insert(item, name);
            }
        }
        Index { context, product }
    }
}

/// Read every solid item of the file, then every loose face ([`loose_faces`]). An item that
/// cannot be read is listed in [`StepSolids::failures`] with its reason; the others are still
/// read.
pub fn read_solids(ex: &Exchange) -> StepSolids {
    read_solids_in(ex, MAX_FILE_FACES, MAX_SOLIDS)
}

fn read_solids_in(ex: &Exchange, mut budget: usize, max_solids: usize) -> StepSolids {
    let index = Index::new(ex);
    let mut context = index.context.clone();
    let mut product = index.product.clone();
    // (item, loose face); a surface model's closed shells read on their own sit where the
    // model does.
    let mut todo: Vec<(u64, bool)> = Vec::new();
    for (&solid, e) in &ex.entities {
        if !SOLID_ITEMS.contains(&e.name()) {
            continue;
        }
        for piece in items(ex, solid) {
            if piece != solid {
                if let Some(c) = index.context.get(&solid).copied() {
                    context.insert(piece, c);
                }
                if let Some(p) = index.product.get(&solid).cloned() {
                    product.insert(piece, p);
                }
            }
            todo.push((piece, false));
        }
    }
    todo.extend(loose_faces(ex).into_iter().map(|f| (f, true)));
    let mut out = StepSolids::default();
    let mut readers: BTreeMap<u64, (Reader<'_>, units::Context)> = BTreeMap::new();
    for (n, (item, loose)) in todo.into_iter().enumerate() {
        if n >= max_solids {
            out.failures.push(StepFailure {
                item,
                message: format!("more than {max_solids} solids in one file"),
            });
            continue;
        }
        let ctx = context.get(&item).copied();
        let (reader, c) = readers.entry(ctx.unwrap_or(0)).or_insert_with(|| {
            let c = ctx.map(|c| units::context(ex, c)).unwrap_or_default();
            (Reader::new(ex, c.units), c)
        });
        let c = *c;
        let read = if loose {
            read_loose(ex, reader, item, c, &mut budget)
        } else {
            read_with(ex, reader, item, c, &mut budget)
        };
        match read {
            Ok(mut s) => {
                s.product = product.get(&item).cloned();
                if s.name.is_empty()
                    && let Some(p) = &s.product
                {
                    s.name.clone_from(p);
                }
                out.solids.push(s);
            }
            Err(message) => out.failures.push(StepFailure { item, message }),
        }
    }
    out
}

/// Read one solid item (a solid entity, or a closed shell of a surface model) in the given
/// context.
pub fn read_solid(ex: &Exchange, item: u64, context: units::Context) -> R<StepSolid> {
    let mut reader = Reader::new(ex, context.units);
    let mut budget = MAX_FILE_FACES;
    read_with(ex, &mut reader, item, context, &mut budget)
}

fn read_with<'a>(
    ex: &'a Exchange,
    reader: &mut Reader<'a>,
    item: u64,
    context: units::Context,
    budget: &mut usize,
) -> R<StepSolid> {
    let e = entity(ex, item)?;
    let name = e
        .params()
        .first()
        .and_then(Param::as_str)
        .unwrap_or("")
        .to_string();
    let (ids, kind) = item_shells(ex, item)?;
    let mut t = Topology::new(ex, reader, budget);
    let mut shells = Vec::new();
    for s in ids {
        shells.push(t.shell(s)?);
    }
    let voids = if kind == SolidKind::BrepWithVoids {
        shells.split_off(shells.len().min(1))
    } else {
        Vec::new()
    };
    Ok(t.finish(item, kind, name, context, shells, voids))
}

fn read_loose<'a>(
    ex: &'a Exchange,
    reader: &mut Reader<'a>,
    item: u64,
    context: units::Context,
    budget: &mut usize,
) -> R<StepSolid> {
    let e = entity(ex, item)?;
    let name = e
        .params()
        .first()
        .and_then(Param::as_str)
        .unwrap_or("")
        .to_string();
    let mut t = Topology::new(ex, reader, budget);
    let (face, orientation) = t.face(item, 0)?;
    let mut shell = StepShell {
        entity: None,
        declared_closed: false,
        orientation: true,
        faces: vec![StepOrientedFace {
            entity: item,
            face,
            orientation,
        }],
        closed: false,
    };
    shell.closed = closed(&t.faces, &shell.faces);
    Ok(t.finish(
        item,
        SolidKind::LooseFace,
        name,
        context,
        vec![shell],
        Vec::new(),
    ))
}

/// A surface model's shells the file calls closed and the others, in its order.
fn model_shells(ex: &Exchange, shells: &[u64]) -> (Vec<u64>, Vec<u64>) {
    shells.iter().partition(|s| {
        ex.get(**s)
            .is_some_and(|e| matches!(e.name(), "CLOSED_SHELL" | "ORIENTED_CLOSED_SHELL"))
    })
}

/// The items a solid item is read as: the item itself, but for a surface model with both
/// closed and open shells (a part with a sheet beside it) each closed shell on its own (a
/// solid), then the model (its open shells, a sheet).
pub fn items(ex: &Exchange, item: u64) -> Vec<u64> {
    let shells: Vec<u64> = match ex.get(item) {
        Some(e) if e.name() == "SHELL_BASED_SURFACE_MODEL" => e
            .params()
            .get(1)
            .and_then(Param::as_list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Param::as_ref_id)
            .collect(),
        _ => Vec::new(),
    };
    let (closed, open) = model_shells(ex, &shells);
    if closed.is_empty() || open.is_empty() {
        return vec![item];
    }
    closed.into_iter().chain([item]).collect()
}

/// The shells of a solid item (the outer shell first, then its voids), and what kind it is.
fn item_shells(ex: &Exchange, item: u64) -> R<(Vec<u64>, SolidKind)> {
    let e = entity(ex, item)?;
    let p = e.params();
    match e.name() {
        "MANIFOLD_SOLID_BREP" => Ok((
            vec![rid(p.get(1), "outer shell")?],
            SolidKind::ManifoldSolid,
        )),
        "BREP_WITH_VOIDS" => {
            let mut shells = vec![rid(p.get(1), "outer shell")?];
            let voids = refs(p.get(2), "voids")?;
            if voids.len() >= MAX_SHELLS {
                return Err(format!("{} voids", voids.len()));
            }
            shells.extend(voids);
            Ok((shells, SolidKind::BrepWithVoids))
        }
        "SHELL_BASED_SURFACE_MODEL" => {
            let ids = refs(p.get(1), "surface model shells")?;
            if ids.is_empty() || ids.len() > MAX_SHELLS {
                return Err(format!("a surface model with {} shells", ids.len()));
            }
            // Its closed shells, where it has open ones too, are items of their own ([`items`]).
            let (closed, open) = model_shells(ex, &ids);
            Ok((
                if closed.is_empty() || open.is_empty() {
                    ids
                } else {
                    open
                },
                SolidKind::SurfaceModel,
            ))
        }
        // A closed shell of a surface model that has open shells too ([`items`]).
        "CLOSED_SHELL" | "ORIENTED_CLOSED_SHELL" => Ok((vec![item], SolidKind::ModelClosedShell)),
        n => Err(unsupported(n)),
    }
}

/// The faces a shell lists, whether the file calls it closed, and its orientation (`.F.`: an
/// oriented shell used against its faces).
fn shell_faces(ex: &Exchange, id: u64) -> R<(Vec<u64>, bool, bool)> {
    let e = entity(ex, id)?;
    let p = e.params();
    let (faces, declared_closed, orientation) = match e.name() {
        n @ ("CLOSED_SHELL" | "OPEN_SHELL") => {
            (refs(p.get(1), "shell faces")?, n == "CLOSED_SHELL", true)
        }
        n @ ("ORIENTED_CLOSED_SHELL" | "ORIENTED_OPEN_SHELL") => {
            let inner = rid(p.get(2), "oriented shell element")?;
            let inner_e = entity(ex, inner)?;
            if !matches!(inner_e.name(), "CLOSED_SHELL" | "OPEN_SHELL") {
                return Err(format!(
                    "#{inner}: {}",
                    unsupported(format!("shell {}", inner_e.name()))
                ));
            }
            (
                refs(inner_e.params().get(1), "shell faces")?,
                n == "ORIENTED_CLOSED_SHELL",
                flag(p.get(3), "oriented shell orientation")?,
            )
        }
        n => return Err(format!("#{id}: {}", unsupported(format!("shell {n}")))),
    };
    if faces.is_empty() {
        return Err(format!("shell #{id} has no faces"));
    }
    if faces.len() > MAX_FACES {
        return Err(format!("shell #{id} with {} faces", faces.len()));
    }
    Ok((faces, declared_closed, orientation))
}

/// The `ADVANCED_FACE` or `FACE_SURFACE` (or whatever else) a face entity `id` is under its
/// `ORIENTED_FACE`s.
pub fn face_entity(ex: &Exchange, id: u64) -> R<u64> {
    face_entity_at(ex, id, 0)
}

fn face_entity_at(ex: &Exchange, id: u64, depth: usize) -> R<u64> {
    if depth > MAX_NESTING {
        return Err("oriented faces nested too deeply".into());
    }
    let e = entity(ex, id)?;
    if e.name() == "ORIENTED_FACE" {
        return face_entity_at(
            ex,
            rid(e.params().get(2), "oriented face element")?,
            depth + 1,
        );
    }
    Ok(id)
}

/// The faces no shell lists, in id order, as solidcore's step-import finds them: an
/// `ADVANCED_FACE`, `FACE_SURFACE` or `ORIENTED_FACE` that no `CLOSED_SHELL` or `OPEN_SHELL`
/// lists (itself or under its `ORIENTED_FACE`s), that is not the face of an `ORIENTED_FACE`,
/// and, for an `ORIENTED_FACE`, whose face no shell lists either. Each is a sheet of its own.
pub fn loose_faces(ex: &Exchange) -> Vec<u64> {
    let mut listed = BTreeSet::new();
    let mut wrapped = BTreeSet::new();
    for e in ex.entities.values() {
        if e.is_complex() {
            continue;
        }
        match e.name() {
            "CLOSED_SHELL" | "OPEN_SHELL" => {
                let faces = e.params().get(1).and_then(Param::as_list).unwrap_or(&[]);
                for f in faces.iter().filter_map(Param::as_ref_id) {
                    listed.insert(f);
                    if let Ok(base) = face_entity(ex, f) {
                        listed.insert(base);
                    }
                }
            }
            "ORIENTED_FACE" => {
                if let Some(inner) = e.params().get(2).and_then(Param::as_ref_id) {
                    wrapped.insert(inner);
                }
            }
            _ => {}
        }
    }
    let mut out = Vec::new();
    for (&id, e) in &ex.entities {
        if e.is_complex()
            || !matches!(e.name(), "ADVANCED_FACE" | "FACE_SURFACE" | "ORIENTED_FACE")
            || listed.contains(&id)
            || wrapped.contains(&id)
        {
            continue;
        }
        if e.name() == "ORIENTED_FACE" && face_entity(ex, id).is_ok_and(|b| listed.contains(&b)) {
            continue;
        }
        out.push(id);
    }
    out
}

/// Every edge of the shell's faces' edge loops used at least twice among them.
fn closed(faces: &[StepFace], shell: &[StepOrientedFace]) -> bool {
    let mut count: BTreeMap<u32, u32> = BTreeMap::new();
    for of in shell {
        for b in faces.get(of.face as usize).map_or(&[][..], |f| &f.bounds) {
            if let StepLoop::Edges { edges, .. } = &b.r#loop {
                for oe in edges {
                    let n = count.entry(oe.edge).or_default();
                    *n = n.saturating_add(1);
                }
            }
        }
    }
    count.values().all(|n| *n >= 2)
}

/// One solid's topology as it is read: each entity read once, indexed in the order first used.
struct Topology<'b, 'a> {
    ex: &'a Exchange,
    reader: &'b mut Reader<'a>,
    budget: &'b mut usize,
    vertex_ids: BTreeMap<u64, u32>,
    vertices: Vec<StepVertex>,
    edge_ids: BTreeMap<u64, u32>,
    edges: Vec<StepEdge>,
    face_ids: BTreeSet<u64>,
    faces: Vec<StepFace>,
}

impl<'b, 'a> Topology<'b, 'a> {
    fn new(ex: &'a Exchange, reader: &'b mut Reader<'a>, budget: &'b mut usize) -> Self {
        Topology {
            ex,
            reader,
            budget,
            vertex_ids: BTreeMap::new(),
            vertices: Vec::new(),
            edge_ids: BTreeMap::new(),
            edges: Vec::new(),
            face_ids: BTreeSet::new(),
            faces: Vec::new(),
        }
    }

    fn finish(
        self,
        item: u64,
        kind: SolidKind,
        name: String,
        context: units::Context,
        shells: Vec<StepShell>,
        voids: Vec<StepShell>,
    ) -> StepSolid {
        StepSolid {
            item,
            kind,
            name,
            product: None,
            units: context.units,
            tolerance: context.uncertainty,
            shells,
            voids,
            faces: self.faces,
            edges: self.edges,
            vertices: self.vertices,
        }
    }

    fn vertex(&mut self, id: u64) -> R<u32> {
        if let Some(v) = self.vertex_ids.get(&id) {
            return Ok(*v);
        }
        let e = entity(self.ex, id)?;
        if e.name() != "VERTEX_POINT" {
            return Err(format!(
                "#{id}: {}",
                unsupported(format!("vertex {}", e.name()))
            ));
        }
        let point = self
            .reader
            .point(rid(e.params().get(1), "vertex point")?)
            .map_err(|m| format!("vertex #{id}: {m}"))?;
        let i = index(self.vertices.len())?;
        self.vertices.push(StepVertex { entity: id, point });
        self.vertex_ids.insert(id, i);
        Ok(i)
    }

    /// An `EDGE_CURVE`'s index.
    fn edge(&mut self, id: u64) -> R<u32> {
        if let Some(e) = self.edge_ids.get(&id) {
            return Ok(*e);
        }
        let e = entity(self.ex, id)?;
        if e.name() != "EDGE_CURVE" {
            return Err(format!(
                "#{id}: {}",
                unsupported(format!("edge {}", e.name()))
            ));
        }
        let p = e.params();
        let (vs, ve) = (rid(p.get(1), "edge start")?, rid(p.get(2), "edge end")?);
        let (start, end) = (self.vertex(vs)?, self.vertex(ve)?);
        let same_sense = flag(p.get(4), "edge same_sense")?;
        let (curve_entity, curve) = self
            .edge_geometry(rid(p.get(3), "edge geometry")?, 0)
            .map_err(|m| format!("edge #{id}: {m}"))?;
        let i = index(self.edges.len())?;
        self.edges.push(StepEdge {
            entity: id,
            start,
            end,
            curve_entity,
            curve,
            same_sense,
        });
        self.edge_ids.insert(id, i);
        Ok(i)
    }

    /// An edge's curve: `SURFACE_CURVE` and its subtypes by their 3D curve, `TRIMMED_CURVE` by
    /// its basis (the edge's vertices bound it). A simple surface curve is followed by
    /// reference, so its pcurves are not read.
    fn edge_geometry(&mut self, id: u64, depth: usize) -> R<(u64, Arc<Curve>)> {
        if depth > MAX_NESTING {
            return Err("edge geometry nested too deeply".into());
        }
        let e = entity(self.ex, id)?;
        if !e.is_complex() {
            match e.name() {
                "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "BOUNDED_SURFACE_CURVE" => {
                    return self.edge_geometry(
                        rid(e.params().get(1), "surface curve 3D curve")?,
                        depth + 1,
                    );
                }
                "TRIMMED_CURVE" => {
                    return self
                        .edge_geometry(rid(e.params().get(1), "trimmed curve basis")?, depth + 1);
                }
                _ => {}
            }
        }
        let mut c = self.reader.curve(id)?;
        // A complex entity read whole: unwrap it the same way.
        for _ in 0..=MAX_NESTING {
            let inner = match &*c {
                Curve::SurfaceCurve { curve_3d, .. } => Arc::clone(curve_3d),
                Curve::Trimmed { basis, .. } => Arc::clone(basis),
                _ => return Ok((id, c)),
            };
            c = inner;
        }
        Err("edge geometry nested too deeply".into())
    }

    /// An `ORIENTED_EDGE` (or a bare `EDGE_CURVE`) as a loop uses it: its edge and whether the
    /// loop runs from the edge's start to its end.
    fn oriented_edge(&mut self, id: u64, depth: usize) -> R<(u32, bool)> {
        if depth > MAX_NESTING {
            return Err("oriented edges nested too deeply".into());
        }
        let e = entity(self.ex, id)?;
        if e.name() == "ORIENTED_EDGE" {
            let p = e.params();
            let (edge, along) =
                self.oriented_edge(rid(p.get(3), "oriented edge element")?, depth + 1)?;
            let orientation = flag(p.get(4), "oriented edge orientation")?;
            return Ok((edge, along == orientation));
        }
        Ok((self.edge(id)?, true))
    }

    /// A face bound, its loop checked connected and closed through its vertices as written.
    /// `edges_so_far` counts the face's edges over its bounds.
    fn bound(&mut self, id: u64, edges_so_far: &mut usize) -> R<StepBound> {
        let b = entity(self.ex, id)?;
        let outer = match b.name() {
            "FACE_OUTER_BOUND" => true,
            "FACE_BOUND" => false,
            n => return Err(format!("#{id}: {}", unsupported(format!("face bound {n}")))),
        };
        let p = b.params();
        let lp_id = rid(p.get(1), "face bound loop")?;
        let orientation = flag(p.get(2), "face bound orientation")?;
        let lp = entity(self.ex, lp_id)?;
        let r#loop = match lp.name() {
            "EDGE_LOOP" => {
                let ids = refs(lp.params().get(1), "edge loop")?;
                *edges_so_far = edges_so_far.saturating_add(ids.len());
                if ids.is_empty() || *edges_so_far > MAX_LOOP {
                    return Err(format!("edge loop #{lp_id} with {} edges", ids.len()));
                }
                let mut edges = Vec::new();
                for oe in ids {
                    let (edge, orientation) = self.oriented_edge(oe, 0)?;
                    edges.push(StepOrientedEdge {
                        entity: oe,
                        edge,
                        orientation,
                    });
                }
                // Connected and closed, as the loop is written.
                let mut last = None;
                for u in &edges {
                    let [s, e] = ends(&self.edges, *u)?;
                    if last.is_some_and(|l| l != s) {
                        return Err(format!("edge loop #{lp_id} is not connected"));
                    }
                    last = Some(e);
                }
                if let (Some(first), Some(last)) = (edges.first(), last)
                    && ends(&self.edges, *first)?[0] != last
                {
                    return Err(format!("edge loop #{lp_id} is not closed"));
                }
                StepLoop::Edges {
                    entity: lp_id,
                    edges,
                }
            }
            "VERTEX_LOOP" => StepLoop::Vertex {
                entity: lp_id,
                vertex: self.vertex(rid(lp.params().get(1), "vertex loop")?)?,
            },
            n => return Err(format!("#{lp_id}: {}", unsupported(format!("loop {n}")))),
        };
        Ok(StepBound {
            entity: id,
            outer,
            orientation,
            r#loop,
        })
    }

    /// A face as a shell uses it: its index and the orientations of its `ORIENTED_FACE`s,
    /// folded.
    fn face(&mut self, id: u64, depth: usize) -> R<(u32, bool)> {
        if depth > MAX_NESTING {
            return Err("oriented faces nested too deeply".into());
        }
        let e = entity(self.ex, id)?;
        let p = e.params();
        let kind = match e.name() {
            "ORIENTED_FACE" => {
                let (f, along) = self.face(rid(p.get(2), "oriented face element")?, depth + 1)?;
                return Ok((f, along == flag(p.get(3), "oriented face orientation")?));
            }
            "ADVANCED_FACE" => FaceKind::Advanced,
            "FACE_SURFACE" => FaceKind::FaceSurface,
            n => return Err(format!("#{id}: {}", unsupported(format!("face {n}")))),
        };
        if !self.face_ids.insert(id) {
            return Err(format!("face #{id} is used twice"));
        }
        *self.budget = self
            .budget
            .checked_sub(1)
            .ok_or_else(|| format!("more than {MAX_FILE_FACES} faces in one file"))?;
        if self.faces.len() >= MAX_FACES {
            return Err(format!("more than {MAX_FACES} faces"));
        }
        let surface_entity = rid(p.get(2), "face geometry")?;
        let surface = self
            .reader
            .surface(surface_entity)
            .map_err(|m| format!("face #{id}: {m}"))?;
        let same_sense = flag(p.get(3), "face same_sense")?;
        let ids = refs(p.get(1), "face bounds")?;
        // ISO 10303-42 gives a face one bound at least (`SET [1:?]`).
        if ids.is_empty() {
            return Err(format!("face #{id} has no bounds"));
        }
        if ids.len() > MAX_BOUNDS {
            return Err(format!("face #{id} with {} bounds", ids.len()));
        }
        let mut bounds = Vec::new();
        let mut edges = 0;
        for b in ids {
            bounds.push(
                self.bound(b, &mut edges)
                    .map_err(|m| format!("face #{id}: {m}"))?,
            );
        }
        let i = index(self.faces.len())?;
        self.faces.push(StepFace {
            entity: id,
            kind,
            surface_entity,
            surface,
            same_sense,
            bounds,
        });
        Ok((i, true))
    }

    /// Shell `id` as the solid lists it.
    fn shell(&mut self, id: u64) -> R<StepShell> {
        let (ids, declared_closed, orientation) = shell_faces(self.ex, id)?;
        let mut faces = Vec::new();
        for f in ids {
            let (face, o) = self.face(f, 0)?;
            faces.push(StepOrientedFace {
                entity: f,
                face,
                orientation: o,
            });
        }
        let closed = closed(&self.faces, &faces);
        Ok(StepShell {
            entity: Some(id),
            declared_closed,
            orientation,
            faces,
            closed,
        })
    }
}

#[cfg(test)]
mod tests;
