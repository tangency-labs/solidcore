// SPDX-License-Identifier: Apache-2.0
// Real public fillet API, preserving the selected pre-operation edge and the
// freshly traversed post-operation solid. No executive fixture owns fake objects.
#include "support/primitive_record.h"
#include <SmApiQueries.h>
#include <SmLoop.h>
#include <cmath>
#include <limits>

namespace fillet_real_solid {
inline double status(SmStatus value) { return value == SM_SUCCESS ? 0.0 : double(value); }
inline void section(std::vector<double>& out, const std::vector<double>& data) {
    out.push_back(double(data.size())); out.insert(out.end(), data.begin(), data.end());
}
template<class T> double label(const SmTArray<T*>& list, T* value) {
    if (!value) return -1.0;
    for (ULONG i=0; i<list.GetSize(); ++i) if (list[i] == value) return double(i);
    return -2.0;
}
inline void topology(std::vector<double>& out, SmBrep* b) {
    SmTArray<SmVertex*> vertices; b->GetVertices(vertices);
    SmTArray<SmEdge*> edges; b->GetEdges(edges);
    SmTArray<SmFace*> faces; b->GetFaces(faces);
    out.push_back(double(edges.GetSize()));
    for (ULONG i=0; i<edges.GetSize(); ++i) {
        auto* edge=edges[i];
        out.push_back(edge->GetCurve() ? double(edge->GetCurve()->GetType()) : 0.0);
        SmTArray<SmVertex*> ends; edge->GetVertices(ends);
        out.push_back(double(ends.GetSize()));
        for (ULONG j=0; j<ends.GetSize(); ++j) out.push_back(label(vertices,ends[j]));
        SmTArray<SmEdgeuse*> uses; edge->GetEdgeuses(uses);
        out.push_back(double(uses.GetSize()));
        for (ULONG j=0; j<uses.GetSize(); ++j)
            out.insert(out.end(), {label(faces,uses[j]->GetFace()),double(uses[j]->GetOrientation())});
    }
    out.push_back(double(faces.GetSize()));
    for (ULONG i=0; i<faces.GetSize(); ++i) {
        auto* face=faces[i];
        out.push_back(face->GetSurface() ? double(face->GetSurface()->GetType()) : 0.0);
        SmTArray<SmLoop*> loops; face->GetLoops(loops);
        out.push_back(double(loops.GetSize()));
        for (ULONG j=0; j<loops.GetSize(); ++j) {
            SmTArray<SmEdgeuse*> uses; loops[j]->GetEdgeuses(uses);
            out.push_back(double(uses.GetSize()));
            for (ULONG q=0; q<uses.GetSize(); ++q)
                out.insert(out.end(), {label(edges,uses[q]->GetEdge()),double(uses[q]->GetOrientation())});
        }
    }
}
inline void counts(std::vector<double>& out, SmBrep* b) {
    SmTArray<SmVertex*> vertices; b->GetVertices(vertices);
    SmTArray<SmEdge*> edges; b->GetEdges(edges);
    SmTArray<SmFace*> faces; b->GetFaces(faces);
    out.insert(out.end(), {double(vertices.GetSize()), double(edges.GetSize()), double(faces.GetSize())});
}
}

// Wire: dimensions XYZ, radius, xsect, continuity, thumbweight, edge count,
// then owning-list edge indices. This initial fixture constructor is a box.
SC_EXPORT int32_t sc_blend_fillet_executive_p3_real_solid(
    const double* input, int64_t n, double* out, int64_t cap, int64_t* out_len) {
    if (!input || n < 9) return SC_ERR_ARG;
    for (int64_t i=0; i<n; ++i) if (!std::isfinite(input[i])) return SC_ERR_ARG;
    for (int i=0; i<3; ++i) if (input[i] <= SM_EFF_ZERO) return SC_ERR_ARG;
    // The exclusive power-of-two bound is exact even when ULONG_MAX or
    // long double would round to 2^64. Check before the integer conversion.
    const double ulong_limit=std::ldexp(1.0,std::numeric_limits<ULONG>::digits);
    for (int i : {4,5}) if (input[i] < 0 || std::floor(input[i]) != input[i]
        || input[i] >= ulong_limit) return SC_ERR_ARG;
    if (input[7] < 1 || std::floor(input[7]) != input[7]
        || input[7] != double(n-8)) return SC_ERR_ARG;
    SC_TRY
    SmBrep* raw=nullptr;
    const auto setup=SmApiCreateBox(SmVector3d(0,0,0),input[0],input[1],input[2],raw);
    if (setup != SM_SUCCESS || !raw) return SC_ERR_ARG;
    std::unique_ptr<SmBrep> body(raw);
    SmTArray<SmEdge*> edges, selected; body->GetEdges(edges);
    std::vector<double> before;
    for (int64_t i=8; i<n; ++i) {
        if (input[i] < 0 || input[i] >= edges.GetSize() || std::floor(input[i]) != input[i]) return SC_ERR_ARG;
        auto* edge=edges[ULONG(input[i])]; selected.Add(edge);
        const auto interval=edge->GetInterval();
        before.insert(before.end(), {input[i],interval.GetMin(),interval.GetMax(),double(edge->GetCurve()!=nullptr)});
        if (auto* curve=edge->GetCurve()) {
            before.push_back(double(curve->GetType()));
            for (double fraction : {0.17,0.37,0.83}) {
                SmPoint3d p; const auto s=curve->EvaluatePoint(interval.Evaluate(fraction),p);
                primitive_record::sample(before,s,p);
            }
        }
    }
    const auto result=SmApiFilletEdges(body.get(),selected,input[3],ULONG(input[4]),ULONG(input[5]),input[6]);
    // Selected pointers may no longer be live; all result traversals are fresh.
    std::vector<double> values; fillet_real_solid::section(values,before);
    values.push_back(fillet_real_solid::status(result));
    fillet_real_solid::counts(values,body.get());
    std::vector<double> topology; fillet_real_solid::topology(topology,body.get());
    fillet_real_solid::section(values,topology);
    std::vector<double> geometry; primitive_record::append(geometry,body.get());
    fillet_real_solid::section(values,geometry);
    double area,volume,mass; SmPoint3d centroid; SmVector3d inertia,products;
    const auto mass_status=SmApiBrepComputeMassProperties(body.get(),1e-4,1.0,SmPoint3d(0,0,0),area,volume,mass,centroid,inertia,products);
    values.push_back(fillet_real_solid::status(mass_status));
    if (mass_status == SM_SUCCESS) {
        values.insert(values.end(), {area,volume,mass});
        primitive_record::point(values,centroid);
        primitive_record::point(values,inertia);
        primitive_record::point(values,products);
    }
    return sc_copy_out(values.data(),int64_t(values.size()),out,cap,out_len);
    SC_CATCH
}

#include <SmFilletCorner.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmVertexuse.h>

namespace fillet_real_solid {
// The corner-vertex component of CreateFilletBrep, cpp1940–2059. Its
// smooth-sequence prelude and later adjustment loop are omitted on BOTH sides.
SmStatus corner_vertices_component(SmFilletExecutive& exec) {
    for(ULONG i=0;i<exec.GetFilletCorners().GetSize();++i) {
        auto* corner=exec.GetFilletCorners()[i];
        const auto fail=[&](const TCHAR* comment) {
            TCHAR name[SM_TBLOCK_SIZE];corner->GetName(name,SM_TBLOCK_SIZE);
            exec.NoteFilletError(SM_FILERR_VERTEX_PROBLEM,corner->GetCornerType(),SM_FIL_FAILURE,
                name,corner->GetFilletedVertex(),nullptr,nullptr,nullptr,nullptr,comment);
            return SM_ERR;
        };
        if(corner->MakeCornerTopology()!=SM_SUCCESS)
            return fail(_T("Fillet Corner Error: Unable to Make Topology"));
        SmBoolean constant=corner->IsBlendingCorner();
        if(constant && corner->IsTangentCorner()) {
            SmTArray<SmFilletSolver*> solvers;corner->GetFilletSolvers(solvers);
            if(solvers.GetSize()==2) {
                auto* curve=solvers[0]->GetEdgeuse(0)->GetEdge()->GetCurve();
                SmPoint3d start,end;curve->GetEnds(start,end);
                const auto point=corner->GetFilletedVertex()->GetPoint();
                const auto tolerance=2.0*corner->GetThisApproxTol3d();
                if(start.DistanceBetween(point)>tolerance && end.DistanceBetween(point)>tolerance)
                    constant=FALSE;
            }
        }
        corner->SetCalcAsConstRad(constant);
        if(corner->CalcCornerVertGeom()!=SM_SUCCESS)
            return fail(_T("Fillet Corner Error: Unable to Compute Geometry"));
    }
    return SM_SUCCESS;
}
std::vector<double> snapshot(SmFilletExecutive& exec,SmTArray<SmVertex*>& target_vertices,SmTArray<SmEdge*>& target_edges,SmTArray<SmFace*>& target_faces, std::vector<double>& parameters) {
    SmTArray<SmEdgeuse*> target_uses;
    for(ULONG i=0;i<target_edges.GetSize();++i) {SmTArray<SmEdgeuse*> uses;target_edges[i]->GetEdgeuses(uses);target_uses.Append(uses);}
    auto& corners=exec.GetFilletCorners();auto& solvers=exec.GetFilletSolvers();
    auto* pseudo=exec.GetPseudoBrep();SmTArray<SmEdge*> edges;pseudo->GetEdges(edges);SmTArray<SmVertex*> vertices;pseudo->GetVertices(vertices);
    SmTArray<SmVertex*> corner_vertices;
    for(ULONG i=0;i<corners.GetSize();++i) {SmTArray<SmFilletVertex*> fv;corners[i]->GetFilletVertices(fv);for(ULONG j=0;j<fv.GetSize();++j)corner_vertices.Add(fv[j]);}
    SmTArray<SmFilletGeom*> all_geometries;
    for(ULONG i=0;i<solvers.GetSize();++i) {SmTArray<SmFilletGeom*> geoms;solvers[i]->GetFilletGeoms(geoms);all_geometries.Append(geoms);}
    std::vector<double> state={double(exec.GetDoGlobalMerge()),double(solvers.GetSize()),double(corners.GetSize()),double(vertices.GetSize()),double(edges.GetSize())};
    for(ULONG i=0;i<corners.GetSize();++i) {
        auto* c=corners[i];
        parameters.insert(parameters.end(),{c->GetThisApproxTol3d(),c->GetTangencyTolRadians()});
        double vertex_label=-2;
        for(ULONG j=0;j<target_vertices.GetSize();++j) if(target_vertices[j]==c->GetFilletedVertex())vertex_label=double(j);
        SmTArray<SmFilletVertex*> fv;c->GetFilletVertices(fv);SmTArray<SmFilletEdge*> fe;c->GetFilletEdges(fe);
        state.insert(state.end(),{double(c->GetCornerType()),vertex_label,double(fv.GetSize()),double(fe.GetSize())});
        for(ULONG j=0;j<fv.GetSize();++j) {
            auto* v=fv[j];auto* object=v->GetPointClassObject();double object_label=object ? -2 : -1;
            for(ULONG u=0;u<target_uses.GetSize();++u) if(static_cast<SmObject*>(target_uses[u])==object)object_label=double(u);
            state.insert(state.end(),{fillet_real_solid::label(corner_vertices,static_cast<SmVertex*>(v)),fillet_real_solid::label(vertices,static_cast<SmVertex*>(v)),double(v->GetFilletVertexType()),fillet_real_solid::label(corners,v->GetFilletCorner()),double(v->GetPointClass()),object?double(object->GetType()):0,object_label});
        }
        for(ULONG j=0;j<fe.GetSize();++j) state.push_back(fillet_real_solid::label(edges,static_cast<SmEdge*>(fe[j])));
    }
    for(ULONG i=0;i<solvers.GetSize();++i) {
        auto* solver=solvers[i];
        auto* constant=dynamic_cast<SmConstantRadiusFS*>(solver);
        state.insert(state.end(),{double(solver->GetSolverType()),constant?2.0:0.0});
        if(constant)parameters.insert(parameters.end(),{constant->GetFilletRadius(0),constant->GetFilletRadius(1)});
        parameters.insert(parameters.end(),{solver->GetThisApproxTol3d(),solver->GetThisAngTolRad(),solver->GetTangencyTolerance()});
        SmTArray<SmFilletGeom*> geometries;solver->GetFilletGeoms(geometries);state.push_back(double(geometries.GetSize()));
        for(ULONG j=0;j<geometries.GetSize();++j) for(ULONG rail=0;rail<2;++rail)state.push_back(fillet_real_solid::label(edges,static_cast<SmEdge*>(geometries[j]->GetRail(rail))));
    }
    for(ULONG i=0;i<edges.GetSize();++i) {
        auto* e=dynamic_cast<SmFilletEdge*>(edges[i]);
        state.push_back(e?double(e->GetFilletEdgeType()):-1);
        state.push_back(e?fillet_real_solid::label(target_faces,e->GetOriginalFace()):-1);
        state.push_back(e?fillet_real_solid::label(corners,e->GetFilletCorner()):-2);
        SmTArray<SmEdgeuse*> uses;edges[i]->GetEdgeuses(uses);state.push_back(double(uses.GetSize()));
        for(ULONG j=0;j<uses.GetSize();++j) {
            auto* vu=uses[j]->GetVertexuse();
            auto* fillet_use=dynamic_cast<SmFilletEdgeuse*>(uses[j]);
            state.insert(state.end(),{double(uses[j]->GetOrientation()),fillet_real_solid::label(corner_vertices,vu?vu->GetVertex():nullptr),fillet_real_solid::label(vertices,vu?vu->GetVertex():nullptr),fillet_use?fillet_real_solid::label(all_geometries,fillet_use->GetFilletGeom()):-2});
        }
    }
    return state;
}
}

SC_EXPORT int32_t sc_blend_fillet_executive_p3_real_prefix(
    const double* input,int64_t n,int32_t checkpoint,double* out,int64_t cap,int64_t* out_len) {
    if(checkpoint<1 || checkpoint>5 || !input || n<9) return SC_ERR_ARG;
    for(int64_t i=0;i<n;++i) if(!std::isfinite(input[i])) return SC_ERR_ARG;
    for(int i=0;i<3;++i) if(input[i]<=SM_EFF_ZERO) return SC_ERR_ARG;
    const double ulong_limit=std::ldexp(1.0,std::numeric_limits<ULONG>::digits);
    for(int i:{4,5}) if(input[i]<0 || std::floor(input[i])!=input[i] || input[i]>=ulong_limit) return SC_ERR_ARG;
    if(input[7]<1 || std::floor(input[7])!=input[7] || input[7]!=double(n-8))return SC_ERR_ARG;
    SC_TRY
    SmBrep* raw=nullptr;
    if(SmApiCreateBox(SmVector3d(0,0,0),input[0],input[1],input[2],raw)!=SM_SUCCESS || !raw)return SC_ERR_ARG;
    std::unique_ptr<SmBrep> body(raw);
    SmTArray<SmVertex*> vertices;body->GetVertices(vertices);
    SmTArray<SmEdge*> edges,selected;body->GetEdges(edges);
    SmTArray<SmFace*> faces;body->GetFaces(faces);
    std::vector<double> before;
    for(int64_t i=8;i<n;++i) {
        if(input[i]<0 || input[i]>=edges.GetSize() || std::floor(input[i])!=input[i])return SC_ERR_ARG;
        auto* edge=edges[ULONG(input[i])];selected.Add(edge);const auto interval=edge->GetInterval();
        before.insert(before.end(),{input[i],interval.GetMin(),interval.GetMax(),double(edge->GetCurve()!=nullptr)});
        if(auto* curve=edge->GetCurve()) {
            before.push_back(double(curve->GetType()));
            for(double fraction:{0.17,0.37,0.83}) {SmPoint3d p;auto s=curve->EvaluatePoint(interval.Evaluate(fraction),p);primitive_record::sample(before,s,p);}
        }
    }
    SmFilletExecutive exec(*body->GetContext(),body.get());exec.SetDoGlobalMerge(TRUE);
    const auto xsect=input[4]==0?SM_FSG_LINEAR:input[4]==2?SM_FSG_BLEND_CURVE:SM_FSG_CIRCULAR;
    const auto prepare=exec.SetFilletParameters(body.get(),selected,xsect,SM_FS_CONST_RADIUS,input[3],nullptr,-1,-1,-1,-1,ULONG(input[5]),input[6]);
    double corners=-1,topology=-1,vertex_geometry=-1;
    std::vector<double> prior,parameters;
    if(prepare==SM_SUCCESS && checkpoint>=2) {
        const auto c=exec.CreateFilletCorners();corners=fillet_real_solid::status(c);
        if(c==SM_SUCCESS && (checkpoint==3 || checkpoint==4)) {
            // Direct first-corner operation after setup/corner creation. This
            // is not the complete CreateFilletBrep prelude: smooth-sequence
            // processing is omitted on BOTH sides, never reported as success.
            std::vector<double> unused;
            prior=fillet_real_solid::snapshot(exec,vertices,edges,faces,unused);
            if(exec.GetFilletCorners().GetSize()) {
                auto* corner=exec.GetFilletCorners()[0];
                const auto t=corner->MakeCornerTopology();topology=fillet_real_solid::status(t);
                if(t==SM_SUCCESS && checkpoint==4)vertex_geometry=fillet_real_solid::status(corner->CalcCornerVertGeom());
            }
        }
        if(c==SM_SUCCESS && checkpoint==5) {
            std::vector<double> unused;
            prior=fillet_real_solid::snapshot(exec,vertices,edges,faces,unused);
            topology=fillet_real_solid::status(fillet_real_solid::corner_vertices_component(exec));
        }
    }
    auto state=fillet_real_solid::snapshot(exec,vertices,edges,faces,parameters);
    std::vector<double> values;fillet_real_solid::section(values,before);
    values.insert(values.end(),{fillet_real_solid::status(prepare),corners,topology});
    if(checkpoint==4)values.push_back(vertex_geometry);
    fillet_real_solid::section(values,prior);fillet_real_solid::section(values,state);fillet_real_solid::section(values,parameters);
    if(checkpoint==4) {
        std::vector<double> vertex_state,vertex_values;
        if(exec.GetFilletCorners().GetSize()) {
            SmTArray<SmFilletVertex*> vertices;exec.GetFilletCorners()[0]->GetFilletVertices(vertices);
            for(ULONG i=0;i<vertices.GetSize();++i) {
                auto* vertex=vertices[i];
                vertex_state.insert(vertex_state.end(),{double(i),double(vertex->GetStatus()),double(vertex->IsProcessed())});
                if(vertex->IsProcessed()) {
                    primitive_record::point(vertex_values,vertex->GetPoint());
                    vertex_values.insert(vertex_values.end(),{vertex->GetTolerance(),vertex->GetOriginalTParam()});
                }
            }
        }
        fillet_real_solid::section(values,vertex_state);fillet_real_solid::section(values,vertex_values);
    }
    if(checkpoint==5) {
        std::vector<double> flags,vertex_state,geometry;
        for(ULONG c=0;c<exec.GetFilletCorners().GetSize();++c) {
            auto* corner=exec.GetFilletCorners()[c];
            flags.insert(flags.end(),{double(c),double(corner->IsBlendingCorner()),
                double(corner->IsTangentCorner()),double(corner->IsCalcAsConstRad())});
            SmTArray<SmFilletVertex*> vertices;corner->GetFilletVertices(vertices);
            for(ULONG v=0;v<vertices.GetSize();++v) {
                auto* vertex=vertices[v];
                vertex_state.insert(vertex_state.end(),{double(c),double(v),double(vertex->GetStatus()),double(vertex->IsProcessed())});
                if(vertex->IsProcessed()) {
                    primitive_record::point(geometry,vertex->GetPoint());
                    geometry.insert(geometry.end(),{vertex->GetTolerance(),vertex->GetOriginalTParam()});
                }
            }
        }
        fillet_real_solid::section(values,flags);fillet_real_solid::section(values,vertex_state);fillet_real_solid::section(values,geometry);
    }
    return sc_copy_out(values.data(),int64_t(values.size()),out,cap,out_len);
    SC_CATCH
}
