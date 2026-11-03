// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `blend.fillet_standard_solver_p1`.
// Wraps source/SMLib/src/SmFilletStandardSolver.cpp: the SmCircularCrossSectionFSG and
// SmBlendCurveCrossSectionFSG CreateSurfaceFromCurves, SmLinearCrossSectionFSG::CreateSurface,
// SmCurveBasedFS::EvaluateBaseCurve, SmSurfaceSurfaceFS::SetupOffsetExtension, the
// SmConstantRadiusFS surface constructor, SetupOffsetValues and LoadJacobian, the
// SmConstantDistanceFS SetupOffsetValues, LoadJacobian and LoadInitialValues and the
// SmVariableRadiusFS radius queries and SetupOffsetValues. The SMLib fork the difftest is
// pinned to has none of the three cross-section surface generators: their ops, and the
// variable-radius ones that need a generator, return SC_ERR_UNLINKED (the v1.0.1 shim, with its
// dry runs of the generators' out-of-bounds reads, is in this file's history).
//
// One entry point, sc_blend_fillet_standard_solver_p1_eval(op, in, n_in, out, cap, out_len).
// Layouts are documented in crates/solidcore-ffi/src/units/blend_fillet_standard_solver_p1.rs,
// which validates every input so upstream only touches valid memory. `out` receives the op
// status (0 = SM_SUCCESS, else the SmStatus value) and the op's outputs.
//
// The SmConstantDistanceFS constructor takes an edge use: the solver is built on the first edge
// of an SmApiCreateBox box, then a derived class replaces its offset surfaces, orientations,
// distance and extension factor and clears its edge uses.
#include "sc_shim.h"

#include "StdAfx.h"
#include <SmApiPrimitives.h>
#include <SmAxis2Placement.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmContext.h>
#include <SmEdge.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmExtentNd.h>
#include <SmFilletExecutive.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmHermiteCurve.h>
#include <SmLine.h>
#include <SmMatrix.h>
#include <SmNurbsCrv.h>
#include <SmNurbsSrf.h>
#include <SmOffsetSurface.h>
#include <SmPlane.h>
#include <SmSurfaceIntersector.h>
#include <SmTArray.h>
#include <nurbs.h>

#include <vector>

namespace
{
// Ops (keep in sync with crates/solidcore-ffi/src/units/blend_fillet_standard_solver_p1.rs).
enum : int32_t
{
    OP_CIRC = 0,
    OP_BLEND = 1,
    OP_LINEAR = 2,
    OP_BASE = 3,
    OP_SS_EXT = 4,
    OP_CR = 5,
    OP_CR_JAC = 6,
    OP_CD_SETUP = 7,
    OP_CD_JAC = 8,
    OP_CD_INIT = 9,
    OP_VR = 10,
    OP_VR_SETUP = 11,
};

struct Rd
{
    const double* in;
    int64_t at;
    double num() { return in[at++]; }
    int64_t idx() { return (int64_t)in[at++]; }
    bool flag() { return in[at++] != 0.0; }
    SmPoint3d pt3()
    {
        double x = num();
        double y = num();
        double z = num();
        return SmPoint3d(x, y, z);
    }
};

double st(SmStatus s)
{
    return s == SM_SUCCESS ? 0.0 : (double)s;
}

double flg(bool b)
{
    return b ? 1.0 : 0.0;
}

void put3(std::vector<double>& o, const SmVector3d& v)
{
    o.push_back(v.x);
    o.push_back(v.y);
    o.push_back(v.z);
}

// Reads `n p m cpts knots` into a gw_CURVE.
gw_CURVE* read_nurb_curve(Rd& rd)
{
    gw_INDEX n = (gw_INDEX)rd.idx();
    gw_DEGREE p = (gw_DEGREE)rd.idx();
    gw_INDEX m = (gw_INDEX)rd.idx();
    gw_CURVE* cur = sm_AllocateNurbCurve(n, p, m);
    for (gw_INDEX i = 0; i <= n; i++)
    {
        cur->pol->Pw[i].x = rd.num();
        cur->pol->Pw[i].y = rd.num();
        cur->pol->Pw[i].z = rd.num();
        cur->pol->Pw[i].w = rd.num();
    }
    for (gw_INDEX i = 0; i <= m; i++)
        cur->knt->U[i] = rd.num();
    return cur;
}

// A B-spline curve `n p m cpts knots` (SmBSplineCurve(gw_CURVE*, 3, FALSE)).
SmBSplineCurve* read_bcurve(Rd& rd, const SmContext& ctx)
{
    return new (ctx) SmBSplineCurve(read_nurb_curve(rd), 3, FALSE);
}

// A curve `kind ...`: 0 B-spline `oob n p m cpts knots`; 1 line `start end`
// (SmLine(start, end, 3, ctx)); 2 Hermite `p1 d1 p2 d2`.
SmCurve* read_curve(Rd& rd, const SmContext& ctx)
{
    int64_t kind = rd.idx();
    if (kind == 0)
    {
        bool oob = rd.flag();
        SmBSplineCurve* c = read_bcurve(rd, ctx);
        c->SetOutOfBoundsEnabled(oob ? TRUE : FALSE);
        return c;
    }
    if (kind == 1)
    {
        SmPoint3d a = rd.pt3();
        SmPoint3d b = rd.pt3();
        return new (ctx) SmLine(a, b, 3, &ctx);
    }
    SmPoint3d p1 = rd.pt3();
    SmVector3d d1 = rd.pt3();
    SmPoint3d p2 = rd.pt3();
    SmVector3d d2 = rd.pt3();
    return new (ctx) SmHermiteCurve(p1, d1, p2, d2, 3);
}

// A surface `kind ...`: 0 B-spline `n m p q r s cpts u v`; 1 plane `o x y sx sy u0 v0 u1 v1`
// (SmPlane(o, x, y, (sx, sy), domain, ctx)).
SmSurface* read_surface(Rd& rd, const SmContext& ctx)
{
    int64_t kind = rd.idx();
    if (kind == 1)
    {
        SmPoint3d o = rd.pt3();
        SmVector3d x = rd.pt3();
        SmVector3d y = rd.pt3();
        double sx = rd.num();
        double sy = rd.num();
        double u0 = rd.num();
        double v0 = rd.num();
        double u1 = rd.num();
        double v1 = rd.num();
        return new (ctx) SmPlane(o, x, y, SmVector2d(sx, sy), SmExtent2d(u0, v0, u1, v1), &ctx);
    }
    gw_INDEX n = (gw_INDEX)rd.idx();
    gw_INDEX m = (gw_INDEX)rd.idx();
    gw_DEGREE p = (gw_DEGREE)rd.idx();
    gw_DEGREE q = (gw_DEGREE)rd.idx();
    gw_INDEX r = (gw_INDEX)rd.idx();
    gw_INDEX s = (gw_INDEX)rd.idx();
    gw_SURFACE* sur = sm_AllocateNurbSurface(n, m, p, q, r, s);
    for (gw_INDEX i = 0; i <= n; i++)
        for (gw_INDEX j = 0; j <= m; j++)
        {
            sur->net->Pw[i][j].x = rd.num();
            sur->net->Pw[i][j].y = rd.num();
            sur->net->Pw[i][j].z = rd.num();
            sur->net->Pw[i][j].w = rd.num();
        }
    for (gw_INDEX i = 0; i <= r; i++)
        sur->knu->U[i] = rd.num();
    for (gw_INDEX i = 0; i <= s; i++)
        sur->knv->U[i] = rd.num();
    return new (ctx) SmBSplineSurface(sur, FALSE, NULL);
}

// `has [n p m cpts knots]` of a B-spline curve's NURB.
void put_curve(std::vector<double>& o, const SmBSplineCurve* c)
{
    const gw_CURVE* cur = c ? c->GetGwNurbPointer() : NULL;
    o.push_back(flg(cur != NULL));
    if (!cur)
        return;
    gw_INDEX n = cur->pol->n;
    gw_INDEX m = cur->knt->m;
    o.push_back((double)n);
    o.push_back((double)cur->p);
    o.push_back((double)m);
    for (gw_INDEX i = 0; i <= n; i++)
    {
        o.push_back(cur->pol->Pw[i].x);
        o.push_back(cur->pol->Pw[i].y);
        o.push_back(cur->pol->Pw[i].z);
        o.push_back(cur->pol->Pw[i].w);
    }
    for (gw_INDEX i = 0; i <= m; i++)
        o.push_back(cur->knt->U[i]);
}

// `has [n m p q r s cpts u v]` of a B-spline surface's NURB.
void put_surface(std::vector<double>& o, const SmBSplineSurface* s)
{
    const gw_SURFACE* sur = s ? s->GetGwNurbPointer() : NULL;
    o.push_back(flg(sur != NULL));
    if (!sur)
        return;
    gw_INDEX n = sur->net->n;
    gw_INDEX m = sur->net->m;
    gw_INDEX r = sur->knu->m;
    gw_INDEX ss = sur->knv->m;
    o.push_back((double)n);
    o.push_back((double)m);
    o.push_back((double)sur->p);
    o.push_back((double)sur->q);
    o.push_back((double)r);
    o.push_back((double)ss);
    for (gw_INDEX i = 0; i <= n; i++)
        for (gw_INDEX j = 0; j <= m; j++)
        {
            o.push_back(sur->net->Pw[i][j].x);
            o.push_back(sur->net->Pw[i][j].y);
            o.push_back(sur->net->Pw[i][j].z);
            o.push_back(sur->net->Pw[i][j].w);
        }
    for (gw_INDEX i = 0; i <= r; i++)
        o.push_back(sur->knu->U[i]);
    for (gw_INDEX i = 0; i <= ss; i++)
        o.push_back(sur->knv->U[i]);
}

// The type queries of a solver: GetType, `len chars` of GetTypeString, IsKindOf of
// `nk` codes.
void put_types(std::vector<double>& o, const SmObject* obj, Rd& rd)
{
    o.push_back((double)obj->GetType());
    const TCHAR* s = obj->GetTypeString();
    std::vector<double> chars;
    for (const TCHAR* c = s; *c; c++)
        chars.push_back((double)*c);
    o.push_back((double)chars.size());
    o.insert(o.end(), chars.begin(), chars.end());
    int64_t nk = rd.idx();
    for (int64_t i = 0; i < nk; i++)
        o.push_back(flg(obj->IsKindOf((SM_TYPE)rd.idx()) != FALSE));
}

// ---- SmCurveBasedFS -----------------------------------------------------------------------

struct CurveBased : SmCurveBasedFS
{
    CurveBased(const SmContext& ctx) : SmCurveBasedFS(ctx, 1.0e-6, 30.0 * SM_PI / 180.0, 2.0 * SM_PI / 180.0)
    {
    }
    SmStatus eval(double t, ULONG n, SmVector3d* out) { return EvaluateBaseCurve(t, n, out); }
    const SmBSplineCurve* extended() const { return m_pExtended; }
};

// `curve ncalls` then `ncalls` x `t nderiv`, `nk kinds[nk]`.
int32_t base(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    SmCurve* c = read_curve(rd, ctx);
    CurveBased* fs = new (ctx) CurveBased(ctx);
    fs->SetOriginalCurve(c);
    int64_t calls = rd.idx();
    o.push_back(0.0);
    for (int64_t i = 0; i < calls; i++)
    {
        double t = rd.num();
        ULONG n = (ULONG)rd.idx();
        std::vector<SmVector3d> pv(n + 1);
        SmStatus s = fs->eval(t, n, pv.data());
        o.push_back(st(s));
        if (s == SM_SUCCESS)
            for (ULONG j = 0; j <= n; j++)
                put3(o, pv[j]);
    }
    put_curve(o, fs->extended());
    o.push_back((double)fs->GetSolverType());
    put_types(o, fs, rd);
    fs->SetOriginalCurve(NULL);
    delete fs;
    delete c;
    return SC_OK;
}

// ---- SmSurfaceSurfaceFS -------------------------------------------------------------------

// `ext surface u0 v0 u1 v1 nk kinds[nk]`.
int32_t ss_ext(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    double ext = rd.num();
    SmSurface* s = read_surface(rd, ctx);
    double u0 = rd.num();
    double v0 = rd.num();
    double u1 = rd.num();
    double v1 = rd.num();
    SmBSplineSurface* sa = NULL;
    SmBSplineSurface* sb = NULL;
    SmBSplineSurface::CreateBilinearSurface(ctx, SmPoint3d(0, 0, 0), SmPoint3d(1, 0, 0),
                                            SmPoint3d(0, 1, 0), SmPoint3d(1, 1, 0), sa);
    SmBSplineSurface::CreateBilinearSurface(ctx, SmPoint3d(0, 0, 0), SmPoint3d(0, 0, 1),
                                            SmPoint3d(1, 0, 0), SmPoint3d(1, 0, 1), sb);
    if (sa == NULL || sb == NULL)
    {
        delete sa;
        delete sb;
        delete s;
        return SC_ERR_ARG;
    }
    SmSurfaceSurfaceFS* fs = new (ctx) SmSurfaceSurfaceFS(
        ctx, 1.0e-6, 30.0 * SM_PI / 180.0, 2.0 * SM_PI / 180.0, 1.0, 1.0, *sa, *sb, TRUE, TRUE,
        SmPoint2d(0.5, 0.5), SmPoint2d(0.5, 0.5));
    fs->SetExtensionFactor(ext);
    SmExtent2d d(u0, v0, u1, v1);
    fs->SetupOffsetExtension(s, d);
    o.push_back(0.0);
    o.push_back(d.GetMin().x);
    o.push_back(d.GetMin().y);
    o.push_back(d.GetMax().x);
    o.push_back(d.GetMax().y);
    o.push_back((double)fs->GetSolverType());
    put_types(o, fs, rd);
    delete fs;
    delete sa;
    delete sb;
    delete s;
    return SC_OK;
}

// ---- SmConstantRadiusFS -------------------------------------------------------------------

struct CrBuilt
{
    SmSurface* a = NULL;
    SmSurface* b = NULL;
    SmConstantRadiusFS* fs = NULL;
    SmSurface* copy = NULL;
    SmOffsetSurface* extra_a = NULL;
    SmOffsetSurface* extra_copy = NULL;
    ~CrBuilt()
    {
        delete extra_a;
        delete extra_copy;
        delete fs;
        delete copy;
        delete a;
        if (b != a)
            delete b;
    }
    // `same surfaceA [surfaceB] r0 r1 o0 o1`.
    bool build(Rd& rd, const SmContext& ctx)
    {
        bool same = rd.flag();
        a = read_surface(rd, ctx);
        b = same ? a : read_surface(rd, ctx);
        double r0 = rd.num();
        double r1 = rd.num();
        bool o0 = rd.flag();
        bool o1 = rd.flag();
        fs = new (ctx) SmConstantRadiusFS(ctx, 1.0e-6, 30.0 * SM_PI / 180.0, 2.0 * SM_PI / 180.0, r0, r1,
                                          *a, *b, o0 ? TRUE : FALSE, o1 ? TRUE : FALSE);
        if (a->Copy(ctx, copy) != SM_SUCCESS || copy == NULL)
            return false;
        extra_a = new (ctx) SmOffsetSurface(r0, *a);
        extra_copy = new (ctx) SmOffsetSurface(r0, *copy);
        return true;
    }
    // Selector 0 NULL, 1 GetSurface(0), 2 GetSurface(1), 3 offset of base A, 4 offset of a copy.
    SmSurface* select(int64_t sel)
    {
        switch (sel)
        {
        case 1: return fs->GetSurface(0);
        case 2: return fs->GetSurface(1);
        case 3: return extra_a;
        case 4: return extra_copy;
        default: return NULL;
        }
    }
};

// `solver set_offsets d0 d1 do_setup user nk kinds[nk]`.
int32_t cr(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    CrBuilt b;
    if (!b.build(rd, ctx))
        return SC_ERR_ARG;
    SmConstantRadiusFS* fs = b.fs;
    bool set = rd.flag();
    double d0 = rd.num();
    double d1 = rd.num();
    bool setup = rd.flag();
    double user = rd.num();
    if (set)
    {
        fs->GetSurface(0)->SetOffsetDistance(d0);
        fs->GetSurface(1)->SetOffsetDistance(d1);
    }
    double dist = -999.0;
    SmStatus s = SM_SUCCESS;
    if (setup)
    {
        SmTsectPnt tp;
        tp.m_adUserDoubles[0] = user;
        s = fs->SetupOffsetValues(tp, &dist);
    }
    o.push_back(st(s));
    o.push_back(dist);
    o.push_back(fs->GetFilletRadius(0));
    o.push_back(fs->GetFilletRadius(1));
    o.push_back(fs->GetOrientation(0));
    o.push_back(fs->GetOrientation(1));
    o.push_back(fs->GetSurface(0)->GetOffsetDistance());
    o.push_back(fs->GetSurface(1)->GetOffsetDistance());
    o.push_back(fs->GetConversionTol());
    o.push_back(fs->GetExtensionFactor());
    o.push_back((double)fs->GetJacobianSize());
    o.push_back(flg(fs->OffsetRadiiCanChange() != FALSE));
    o.push_back((double)fs->GetSolverType());
    put_types(o, fs, rd);
    return SC_OK;
}

// `solver rail ne np sel1 off1 sel2 off2 nx x[nx] nf f[nf] has_j rows cols`.
int32_t cr_jac(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    CrBuilt b;
    if (!b.build(rd, ctx))
        return SC_ERR_ARG;
    ULONG rail = (ULONG)rd.idx();
    ULONG ne = (ULONG)rd.idx();
    ULONG np = (ULONG)rd.idx();
    SmSurface* s1 = b.select(rd.idx());
    ULONG off1 = (ULONG)rd.idx();
    SmSurface* s2 = b.select(rd.idx());
    ULONG off2 = (ULONG)rd.idx();
    int64_t nx = rd.idx();
    SmTArray<double> x;
    for (int64_t i = 0; i < nx; i++)
        x.Add(rd.num());
    int64_t nf = rd.idx();
    SmTArray<double> f;
    for (int64_t i = 0; i < nf; i++)
        f.Add(rd.num());
    bool has_j = rd.flag();
    ULONG rows = (ULONG)rd.idx();
    ULONG cols = (ULONG)rd.idx();
    SmMatrix j(has_j ? rows : 1, has_j ? cols : 1);
    if (has_j)
        for (ULONG r = 0; r < rows; r++)
            for (ULONG c = 0; c < cols; c++)
                j[r][c] = 0.0;
    SmBoolean found = TRUE;
    SmStatus stt = b.fs->LoadJacobian(ne, np, rail, s1, off1, s2, off2, x, f, has_j ? &j : NULL, found);
    o.push_back(st(stt));
    o.push_back((double)ne);
    o.push_back((double)np);
    o.push_back((double)off1);
    o.push_back((double)off2);
    o.push_back(flg(found != FALSE));
    for (ULONG i = 0; i < f.GetSize(); i++)
        o.push_back(f[i]);
    if (has_j)
        for (ULONG r = 0; r < rows; r++)
            for (ULONG c = 0; c < cols; c++)
                o.push_back(j[r][c]);
    o.push_back(b.fs->GetSurface(0)->GetOffsetDistance());
    o.push_back(b.fs->GetSurface(1)->GetOffsetDistance());
    return SC_OK;
}

// ---- SmConstantDistanceFS -----------------------------------------------------------------

struct Cd : SmConstantDistanceFS
{
    Cd(const SmContext& ctx, double dist, SmEdgeuse* eu)
        : SmConstantDistanceFS(ctx, 1.0e-6, 30.0 * SM_PI / 180.0, 2.0 * SM_PI / 180.0, dist, eu)
    {
    }
    void setup(SmOffsetSurface* s0, SmOffsetSurface* s1, double o0, double o1, double dist, double ext)
    {
        delete this->m_pSurfaces[0];
        delete this->m_pSurfaces[1];
        this->m_pSurfaces[0] = s0;
        this->m_pSurfaces[1] = s1;
        this->m_dOrientations[0] = o0;
        this->m_dOrientations[1] = o1;
        this->m_dDistance = dist;
        this->m_pEdgeuses[0] = NULL;
        this->m_pEdgeuses[1] = NULL;
        this->m_dSurfaceExtensionFactor = ext;
    }
};

// The solver block `same surfaceA [surfaceB] d0 d1 o0 o1 dist ext`: rail i offsets its base by
// d_i with orientation o_i.
struct CdBuilt
{
    SmBrep* brep = NULL;
    SmSurface* a = NULL;
    SmSurface* b = NULL;
    SmSurface* copy = NULL;
    Cd* fs = NULL;
    SmOffsetSurface* extra_a = NULL;
    SmOffsetSurface* extra_copy = NULL;
    ~CdBuilt()
    {
        delete extra_a;
        delete extra_copy;
        delete fs;
        delete copy;
        delete a;
        if (b != a)
            delete b;
        delete brep;
    }
    bool build(Rd& rd, const SmContext& ctx)
    {
        bool same = rd.flag();
        a = read_surface(rd, ctx);
        b = same ? a : read_surface(rd, ctx);
        double d0 = rd.num();
        double d1 = rd.num();
        double o0 = rd.num();
        double o1 = rd.num();
        double dist = rd.num();
        double ext = rd.num();
        if (SmApiCreateBox(SmVector3d(0.0, 0.0, 0.0), 1.0, 1.0, 1.0, brep) != SM_SUCCESS || brep == NULL)
            return false;
        SmTArray<SmEdge*> edges;
        brep->GetEdges(edges);
        if (edges.GetSize() == 0)
            return false;
        SmEdgeuse* eu = edges[0]->GetPrimaryEdgeuse();
        if (eu == NULL)
            return false;
        fs = new (ctx) Cd(ctx, dist, eu);
        SmOffsetSurface* s0 = new (ctx) SmOffsetSurface(d0, *a);
        SmOffsetSurface* s1 = new (ctx) SmOffsetSurface(d1, *b);
        s0->SetAllowExtension(TRUE);
        s1->SetAllowExtension(TRUE);
        fs->setup(s0, s1, o0, o1, dist, ext);
        if (a->Copy(ctx, copy) != SM_SUCCESS || copy == NULL)
            return false;
        extra_a = new (ctx) SmOffsetSurface(d0, *a);
        extra_copy = new (ctx) SmOffsetSurface(d0, *copy);
        return true;
    }
    SmSurface* select(int64_t sel)
    {
        switch (sel)
        {
        case 1: return fs->GetSurface(0);
        case 2: return fs->GetSurface(1);
        case 3: return extra_a;
        case 4: return extra_copy;
        default: return NULL;
        }
    }
};

// `solver user nk kinds[nk]`.
int32_t cd_setup(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    CdBuilt b;
    if (!b.build(rd, ctx))
        return SC_ERR_ARG;
    double user = rd.num();
    SmTsectPnt tp;
    tp.m_adUserDoubles[0] = user;
    double dist = -999.0;
    SmStatus s = b.fs->SetupOffsetValues(tp, &dist);
    o.push_back(st(s));
    o.push_back(dist);
    o.push_back(b.fs->GetSurface(0)->GetOffsetDistance());
    o.push_back(b.fs->GetSurface(1)->GetOffsetDistance());
    o.push_back(b.fs->GetDistance());
    o.push_back(flg(b.fs->CanMakeAnalytic() != FALSE));
    o.push_back((double)b.fs->GetJacobianSize());
    o.push_back((double)b.fs->GetSolverType());
    put_types(o, b.fs, rd);
    return SC_OK;
}

// `solver rail ne np sel1 off1 sel2 off2 nx x[nx] nf f[nf] has_j rows cols`.
int32_t cd_jac(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    CdBuilt b;
    if (!b.build(rd, ctx))
        return SC_ERR_ARG;
    ULONG rail = (ULONG)rd.idx();
    ULONG ne = (ULONG)rd.idx();
    ULONG np = (ULONG)rd.idx();
    SmSurface* s1 = b.select(rd.idx());
    ULONG off1 = (ULONG)rd.idx();
    SmSurface* s2 = b.select(rd.idx());
    ULONG off2 = (ULONG)rd.idx();
    int64_t nx = rd.idx();
    SmTArray<double> x;
    for (int64_t i = 0; i < nx; i++)
        x.Add(rd.num());
    int64_t nf = rd.idx();
    SmTArray<double> f;
    for (int64_t i = 0; i < nf; i++)
        f.Add(rd.num());
    bool has_j = rd.flag();
    ULONG rows = (ULONG)rd.idx();
    ULONG cols = (ULONG)rd.idx();
    SmMatrix j(has_j ? rows : 1, has_j ? cols : 1);
    if (has_j)
        for (ULONG r = 0; r < rows; r++)
            for (ULONG c = 0; c < cols; c++)
                j[r][c] = 0.0;
    SmBoolean found = TRUE;
    SmStatus stt = b.fs->LoadJacobian(ne, np, rail, s1, off1, s2, off2, x, f, has_j ? &j : NULL, found);
    o.push_back(st(stt));
    o.push_back((double)ne);
    o.push_back((double)np);
    o.push_back(flg(found != FALSE));
    for (ULONG i = 0; i < f.GetSize(); i++)
        o.push_back(f[i]);
    if (has_j)
        for (ULONG r = 0; r < rows; r++)
            for (ULONG c = 0; c < cols; c++)
                o.push_back(j[r][c]);
    o.push_back(b.fs->GetSurface(0)->GetOffsetDistance());
    o.push_back(b.fs->GetSurface(1)->GetOffsetDistance());
    return SC_OK;
}

// `solver do2 rail sel1 s1_index s2_index dim` `dim` x `min max`, `np per[np] ng guess[ng]`.
int32_t cd_init(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    CdBuilt b;
    if (!b.build(rd, ctx))
        return SC_ERR_ARG;
    ULONG do2 = rd.flag() ? 1 : 0;
    ULONG rail = (ULONG)rd.idx();
    SmSurface* s1 = b.select(rd.idx());
    ULONG s1_index = (ULONG)rd.idx();
    ULONG s2_index = (ULONG)rd.idx();
    ULONG dim = (ULONG)rd.idx();
    SmExtentNd ivls(dim);
    for (ULONG i = 0; i < dim; i++)
    {
        double lo = rd.num();
        double hi = rd.num();
        ivls[i] = SmExtent1d(lo, hi);
    }
    int64_t np = rd.idx();
    SmTArray<SmBoolean> per;
    for (int64_t i = 0; i < np; i++)
        per.Add(rd.flag() ? TRUE : FALSE);
    int64_t ng = rd.idx();
    SmTArray<double> g;
    for (int64_t i = 0; i < ng; i++)
        g.Add(rd.num());
    if (s1 == NULL)
        return SC_ERR_ARG;
    SmSurface* s2 = NULL;
    SmStatus stt = b.fs->LoadInitialValues(do2, rail, *s1, s1_index, s2, s2_index, ivls, per, g);
    o.push_back(st(stt));
    double s2_rail = -1.0;
    if (s2 != NULL)
        s2_rail = s2 == b.fs->GetSurface(0) ? 0.0 : s2 == b.fs->GetSurface(1) ? 1.0 : -2.0;
    o.push_back(s2_rail);
    o.push_back((double)s2_index);
    for (ULONG i = 0; i < dim; i++)
    {
        o.push_back(ivls[i].GetMin());
        o.push_back(ivls[i].GetMax());
    }
    o.push_back((double)per.GetSize());
    for (ULONG i = 0; i < per.GetSize(); i++)
        o.push_back(flg(per[i] != FALSE));
    o.push_back((double)g.GetSize());
    for (ULONG i = 0; i < g.GetSize(); i++)
        o.push_back(g[i]);
    o.push_back(b.fs->GetSurface(0)->GetOffsetDistance());
    o.push_back(b.fs->GetSurface(1)->GetOffsetDistance());
    return SC_OK;
}

} // namespace

SC_EXPORT int32_t sc_blend_fillet_standard_solver_p1_eval(int32_t op, const double* in, int64_t n_in, double* out,
                                                          int64_t cap, int64_t* out_len)
{
    (void)n_in;
    SC_TRY
    SmContext ctx;
    std::vector<double> o;
    Rd rd{in, 0};
    int32_t rc = SC_OK;
    switch (op)
    {
    case OP_BASE: rc = base(rd, ctx, o); break;
    case OP_SS_EXT: rc = ss_ext(rd, ctx, o); break;
    case OP_CR: rc = cr(rd, ctx, o); break;
    case OP_CR_JAC: rc = cr_jac(rd, ctx, o); break;
    case OP_CD_SETUP: rc = cd_setup(rd, ctx, o); break;
    case OP_CD_JAC: rc = cd_jac(rd, ctx, o); break;
    case OP_CD_INIT: rc = cd_init(rd, ctx, o); break;
    // Not linked: the SMLib fork has no SmLinearCrossSectionFSG, SmCircularCrossSectionFSG
    // (tangency-labs/smlib 3aae1c4) or SmBlendCurveCrossSectionFSG (0bd3f4d), which these ops
    // run or give the solver (the variable-radius ops read its type in OffsetRadiiCanChange).
    case OP_CIRC:
    case OP_BLEND:
    case OP_LINEAR:
    case OP_VR:
    case OP_VR_SETUP: return SC_ERR_UNLINKED;
    default: return SC_ERR_ARG;
    }
    if (rc != SC_OK)
        return rc;
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
