// SPDX-License-Identifier: Apache-2.0
// topology.brep_p5. Native operations use the real SmBrep, not copied algorithms.
// Shared fixture/export helpers read only SMLib. No shared harness changes.
#include "sc_topology_fixture.h"
#include <memory>

namespace {
// Fixture block: kind, two xyz points, two raw tolerances, two incident-edge counts.
// kind 0: two separate shell vertices with 0..4 outgoing wire edges apiece;
// kind 1: upstream prog_test my_test_topology's wire network; kind 2: planar face.
SmBrep* p5_build(const double* f, SmStatus& status) {
    const SmContext& ctx = the_context();
    // A scaling operation may update its context's model-size tolerance. Start
    // each independent fixture with the default estimate, including buffer retries.
    SmTol::SetModelSizeEstimate(const_cast<SmContext&>(ctx), 50.0);
    if (f[0] == 1) {
        auto* b = new (ctx) SmBrep();
        status = build_wire_network(ctx, b);
        return b;
    }
    if (f[0] == 2) {
        // Same real SmPlane/CreateFaceFromSurface constructor as the shared fixture.
        const double block[] = {2, 1, 0,0,0, 0,0,1, -2,-3,4,5, 0,0};
        Rd rd{block, 0};
        return build(rd, "", status);
    }
    auto* b = new (ctx) SmBrep();
    b->m_bEditingEnabled = TRUE;
    status = SM_SUCCESS;
    for (int i=0; i<2; ++i) {
        const int p = 1 + 3*i;
        SmVertex* v = nullptr;
        SmShell* shell = nullptr;
        status = b->MakeShellVertex(b->GetInfiniteRegion(), SmPoint3d(f[p],f[p+1],f[p+2]), shell, v);
        if (status != SM_SUCCESS) return b;
        v->SetTolerance(f[7+i], FALSE, FALSE);
        for (int j=0; j<(int)f[9+i]; ++j) {
            SmPoint3d end = v->GetPoint() + SmVector3d(2+j,3+j*j,5+j);
            SmBSplineCurve* line = nullptr;
            status = SmBSplineCurve::CreateLineSegment(ctx, 3, v->GetPoint(), end, line);
            if (status != SM_SUCCESS) return b;
            SmVertex* tail = nullptr;
            SmEdge* edge = nullptr;
            status = b->MakeWireEdgeVertex(b->GetInfiniteRegion(), v, line,
                       line->GetNaturalInterval(), SM_OT_SAME, end, edge, tail);
            if (status != SM_SUCCESS) return b;
        }
    }
    return b;
}

template<class T> T* select(const std::vector<T*>& v, double x) {
    if (v.empty() || x < 0) return nullptr;
    return v[(size_t)((int64_t)x % (int64_t)v.size())];
}


}

// Args: editing, first selector, second selector, aux selector, mode, raw tolerance,
// placement O/X/Y (9), have_scale, scale xyz (3), scale_tol, flag, spare.
SC_EXPORT int32_t sc_topology_brep_p5_eval(int32_t op, const double* fixture, int64_t nf,
        const double* a, int64_t na, double* out, int64_t cap, int64_t* out_len) {
    if (nf!=11 || na!=22 || !fixture || !a || !out_len || op<0 || op>13) return SC_ERR_ARG;
    SC_TRY
    SmStatus status=SM_SUCCESS;
    std::unique_ptr<SmBrep> owner(p5_build(fixture,status));
    if (!owner || status!=SM_SUCCESS) {
        std::vector<double> o={-999.0,st(status)};
        return sc_copy_out(o.data(),o.size(),out,cap,out_len);
    }
    SmBrep* b=owner.get();
    b->m_bEditingEnabled=a[0]!=0;
    Graph g; walk(b,g,true);
    std::vector<double> o;
    if (op==0) {
        o.push_back(0.0);
        export_graph(b,g,o);
        return sc_copy_out(o.data(),o.size(),out,cap,out_len);
    }
    SmVertex* v1=select(g.vertices,a[1]);
    SmVertex* v2=select(g.vertices,a[2]);
    SmEdge* edge=select(g.edges,a[3]);
    SmFace* face=select(g.faces,a[1]);
    if (op==1 || op==2) {
        // Equal vertices cause native use-after-free, and nulls are dereferenced.
        if (!v1 || !v2 || v1==v2) return SC_ERR_ARG;
        if (op==2) {
            SmTArray<SmEdge*> e1,e2; v1->GetEdges(e1); v2->GetEdges(e2);
            if (e1.GetSize()!=0 && e2.GetSize()!=0) return SC_ERR_ARG;
        }
        SmEdge* deleted=edge;
        status=b->CombineCoincidentVertices(v1,v2,deleted);
        o.push_back(st(status));
        if (op==2) {
            // Numerical stage only: compare the survivor target against the real
            // mutation's final state. Topology mutation is deliberately not certified.
            put_vec(o,v1->GetPoint()); o.push_back(v1->GetTolerance());
            return sc_copy_out(o.data(),o.size(),out,cap,out_len);
        }
        o.push_back(deleted ? 1.0 : 0.0);
    } else if (op==3) {
        status=b->DeleteTopologicalVertex(v1,edge); o.push_back(st(status));
    } else if (op==4) {
        status=b->FindAndSetInfiniteRegion(nullptr,(SmSetRegionVoidType)(int)a[4],nullptr);
        o.push_back(st(status));
    } else if (op==5) {
        if (!v1 || !v2 || v1==v2) return SC_ERR_ARG;
        status=b->GlueVertices(v1,v2); o.push_back(st(status));
    } else if (op==6) {
        // The disabled-edit and null-face paths are defined before geometry is read.
        if (a[0]!=0 && face) return SC_ERR_ARG;
        SmEdge* ne=edge; SmLoop* nl=select(g.loops,a[2]); SmFace* nf=face;
        status=b->MakeEdgeInFace(face,v1,v2,nullptr,nullptr,SmExtent1d(0,1),SM_OT_SAME,0,ne,nl,nf);
        o.insert(o.end(),{st(status),ne==edge?1.0:0.0,nl==select(g.loops,a[2])?1.0:0.0,nf==face?1.0:0.0});
    } else if (op==7) {
        SmTArray<SmRegion*> regions;
        if (a[4]==2) b->GetRegions(regions);
        status=b->MakeManifold(a[4]==0 ? nullptr : &regions,a[20]!=0);
        o.push_back(st(status));
    } else if (op==8) {
        SmTArray<SmEdge*> ke; SmTArray<SmVertex*> kv;
        if (a[20]!=0) { b->GetEdges(ke); b->GetVertices(kv); }
        status=b->RemoveTopologicalEdgesAndVertices(a[4]==0?nullptr:&ke,a[4]==0?nullptr:&kv);
        o.push_back(st(status));
    } else if (op==9) {
        SmTArray<SmFace*> fs;
        if (a[4]!=0 && face) fs.Add(face);
        const SmContext& ctx=*b->GetContext();
        std::unique_ptr<SmNewMarkAndLock> held;
        if (a[20]!=0) held=std::make_unique<SmNewMarkAndLock>(&ctx,SM_MT_ALLMARKS);
        const SmMarkType marks[]={SM_MT_MARK,SM_MT_MARK2,SM_MT_MARK3,SM_MT_MARKIO,SM_MT_MARKASSERT};
        ULONG before[5];
        for (int j=0;j<5;++j) before[j]=ctx.GetCurrentMark(marks[j]);
        status=b->RemoveFaces(fs); o.push_back(st(status));
        for (int j=0;j<5;++j) {
            o.push_back((double)(ctx.GetCurrentMark(marks[j])-before[j]));
            o.push_back(flg(ctx.IsMarkLocked(marks[j])));
        }
    } else if (op==10) {
        SmBoolean just=a[20]!=0;
        status=b->StitchAndOrient(a[5],a[21]!=0,nullptr,a[4]==0?nullptr:&just);
        o.push_back(st(status));
    } else if (op==11) {
        if (a[0]!=0 && !face) return SC_ERR_ARG;
        SmTArray<SmVertex*> dv; SmTArray<SmEdge*> ne,de;
        if (v1) dv.Add(v1); if (edge) {ne.Add(edge);de.Add(edge);}
        SmXSectTol3d tol=a[5];
        status=b->SqueezeFace(face,a[4]==0?nullptr:&tol,&dv,&ne,&de);
        o.insert(o.end(),{st(status),(double)dv.GetSize(),(double)ne.GetSize(),(double)de.GetSize()});
    } else if (op==12) {
        SmAxis2Placement pose(SmPoint3d(a[6],a[7],a[8]),SmVector3d(a[9],a[10],a[11]),SmVector3d(a[12],a[13],a[14]));
        SmVector3d scale(a[16],a[17],a[18]);
        status=b->Transform(pose,a[15]!=0?&scale:nullptr,a[19]!=0); o.push_back(st(status));
    } else if (op==13) {
        // A real SmLine; the tested interval is deliberately non-normalized.
        SmLine* line=new(the_context()) SmLine(SmPoint3d(0,0,0),SmVector3d(1,0,0),SmExtent1d(-2,3));
        SmEdge* new_edge=nullptr;
        status=b->CreateWireEdgeFromCurve(line,SmExtent1d(-2,3),new_edge);
        o.push_back(st(status));
    }
    std::vector<double> result={(double)o.size()};
    result.insert(result.end(),o.begin(),o.end());
    result.push_back(0.0);
    Graph after; walk(b,after,true);
    export_graph(b,after,result);
    return sc_copy_out(result.data(),result.size(),out,cap,out_len);
    SC_CATCH
}
