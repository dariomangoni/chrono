// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Andrea Favali, Alessandro Tasora
// =============================================================================

#include "chrono/collision/ChCModelBullet.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChElementShellANCF.h"
#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChMesh.h"

using namespace std;

namespace chrono {
namespace fea {

//////////////////////////////////////////////////////////////////////////////
////  ChContactNodeXYZ

void ChContactNodeXYZ::ContactForceLoadResidual_F(const ChVector<>& F,
                                                  const ChVector<>& abs_point,
                                                  ChVectorDynamic<>& R) {
    R.PasteSumVector(F, this->mnode->NodeGetOffset_w() + 0, 0);
}

void ChContactNodeXYZ::ComputeJacobianForContactPart(const ChVector<>& abs_point,
                                                     ChMatrix33<>& contact_plane,
                                                     type_constraint_tuple& jacobian_tuple_N,
                                                     type_constraint_tuple& jacobian_tuple_U,
                                                     type_constraint_tuple& jacobian_tuple_V,
                                                     bool second) {
    ChMatrix33<> Jx1;

    Jx1.CopyFromMatrixT(contact_plane);
    if (!second)
        Jx1.MatrNeg();

    jacobian_tuple_N.Get_Cq()->PasteClippedMatrix(Jx1, 0, 0, 1, 3, 0, 0);
    jacobian_tuple_U.Get_Cq()->PasteClippedMatrix(Jx1, 1, 0, 1, 3, 0, 0);
    jacobian_tuple_V.Get_Cq()->PasteClippedMatrix(Jx1, 2, 0, 1, 3, 0, 0);
}

std::shared_ptr<ChMaterialSurface>& ChContactNodeXYZ::GetMaterialSurface() {
    return container->GetMaterialSurface();
}

ChPhysicsItem* ChContactNodeXYZ::GetPhysicsItem() {
    return (ChPhysicsItem*)container->GetMesh();
}

//////////////////////////////////////////////////////////////////////////////
////  ChContactNodeXYZsphere

ChContactNodeXYZsphere::ChContactNodeXYZsphere(ChNodeFEAxyz* anode, ChContactSurface* acontainer)
    : ChContactNodeXYZ(anode, acontainer) {
    this->collision_model = new collision::ChModelBullet;
    this->collision_model->SetContactable(this);
}

//////////////////////////////////////////////////////////////////////////////
////  ChContactNodeXYZROT

void ChContactNodeXYZROT::ContactForceLoadResidual_F(const ChVector<>& F,
                                                     const ChVector<>& abs_point,
                                                     ChVectorDynamic<>& R) {
    R.PasteSumVector(F, this->mnode->NodeGetOffset_w() + 0, 0);
}

void ChContactNodeXYZROT::ComputeJacobianForContactPart(const ChVector<>& abs_point,
                                                        ChMatrix33<>& contact_plane,
                                                        type_constraint_tuple& jacobian_tuple_N,
                                                        type_constraint_tuple& jacobian_tuple_U,
                                                        type_constraint_tuple& jacobian_tuple_V,
                                                        bool second) {
    ChMatrix33<> Jx1;

    Jx1.CopyFromMatrixT(contact_plane);
    if (!second)
        Jx1.MatrNeg();

    jacobian_tuple_N.Get_Cq()->PasteClippedMatrix(Jx1, 0, 0, 1, 3, 0, 0);
    jacobian_tuple_U.Get_Cq()->PasteClippedMatrix(Jx1, 1, 0, 1, 3, 0, 0);
    jacobian_tuple_V.Get_Cq()->PasteClippedMatrix(Jx1, 2, 0, 1, 3, 0, 0);
}

std::shared_ptr<ChMaterialSurface>& ChContactNodeXYZROT::GetMaterialSurface() {
    return container->GetMaterialSurface();
}

ChPhysicsItem* ChContactNodeXYZROT::GetPhysicsItem() {
    return (ChPhysicsItem*)container->GetMesh();
}

//////////////////////////////////////////////////////////////////////////////
////  ChContactNodeXYZROTsphere

ChContactNodeXYZROTsphere::ChContactNodeXYZROTsphere(ChNodeFEAxyzrot* anode, ChContactSurface* acontainer)
    : ChContactNodeXYZROT(anode, acontainer) {
    this->collision_model = new collision::ChModelBullet;
    this->collision_model->SetContactable(this);
}

//////////////////////////////////////////////////////////////////////////////
////  ChContactSurfaceNodeCloud

void ChContactSurfaceNodeCloud::AddNode(std::shared_ptr<ChNodeFEAxyz> mnode, const double point_radius, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZsphere>(mnode.get(), this);

    newp->GetCollisionModel()->AddPoint(point_radius);
    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes.push_back(newp);
}

void ChContactSurfaceNodeCloud::AddNode(std::shared_ptr<ChNodeFEAxyzrot> mnode, const double point_radius, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZROTsphere>(mnode.get(), this);


    newp->GetCollisionModel()->AddPoint(point_radius);
    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes_rot.push_back(newp);
}


/// Add all nodes of the mesh to this collision cloud
void ChContactSurfaceNodeCloud::AddAllNodes(const double point_radius, int collision_family) {
    if (!this->GetMesh())
        return;
    for (unsigned int i = 0; i < this->GetMesh()->GetNnodes(); ++i)
        if (auto mnodeFEA = std::dynamic_pointer_cast<ChNodeFEAxyz>(this->GetMesh()->GetNode(i)))
            this->AddNode(mnodeFEA, point_radius, collision_family);
        else if (auto mnodeFEArot = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(this->GetMesh()->GetNode(i)))
            this->AddNode(mnodeFEArot, point_radius, collision_family);
}


void ChContactSurfaceNodeCloud::AddNode_Sphere(std::shared_ptr<ChNodeFEAxyz> mnode, const double sphere_radius, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZsphere>(mnode.get(), this);

    newp->GetCollisionModel()->AddSphere(sphere_radius);
    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes.push_back(newp);
}

void ChContactSurfaceNodeCloud::AddNode_Sphere(std::shared_ptr<ChNodeFEAxyzrot> mnode, const double sphere_radius, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZROTsphere>(mnode.get(), this);

    newp->GetCollisionModel()->AddSphere(sphere_radius);
    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes_rot.push_back(newp);
}

/// Add all nodes of the mesh to this collision cloud, but with Spheres as collision shapes
void ChContactSurfaceNodeCloud::AddAllNodes_Spheres(const double sphere_radius, int collision_family) {
    if (!this->GetMesh())
        return;
    for (unsigned int i = 0; i < this->GetMesh()->GetNnodes(); ++i)
        if (auto mnodeFEA = std::dynamic_pointer_cast<ChNodeFEAxyz>(this->GetMesh()->GetNode(i)))
            this->AddNode_Sphere(mnodeFEA, sphere_radius, collision_family);
        else if (auto mnodeFEArot = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(this->GetMesh()->GetNode(i)))
            this->AddNode_Sphere(mnodeFEArot, sphere_radius, collision_family);
}


void ChContactSurfaceNodeCloud::AddNode_Cylinder(std::shared_ptr<ChNodeFEAxyz> mnode, const double cyl_radius, const double cyl_length, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZsphere>(mnode.get(), this);

    ChMatrix33<> rotMat;
    rotMat.Set_A_AngAxis(CH_C_PI_2, VECT_Z);
    newp->GetCollisionModel()->AddCylinder(cyl_radius, cyl_radius, cyl_length, ChVector<>(), rotMat);

    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes.push_back(newp);
}

void ChContactSurfaceNodeCloud::AddNode_Cylinder(std::shared_ptr<ChNodeFEAxyzrot> mnode, const double cyl_radius, const double cyl_length, int collision_family) {
    if (!mnode)
        return;

    auto newp = std::make_shared<ChContactNodeXYZROTsphere>(mnode.get(), this);
    ChMatrix33<> rotMat;
    rotMat.Set_A_AngAxis(CH_C_PI_2, VECT_Z);

    newp->GetCollisionModel()->AddCylinder(cyl_radius, cyl_radius, cyl_length, ChVector<>(), rotMat);
    newp->GetCollisionModel()->BuildModel();  // will also add to system, if collision is on.

    newp->GetCollisionModel()->SetFamily(collision_family);
    newp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(collision_family);

    this->vnodes_rot.push_back(newp);
}

/// Add all nodes of the mesh to this collision cloud, but with Spheres as collision shapes
void ChContactSurfaceNodeCloud::AddAllNodes_Cylinders(const double cyl_radius, const double cyl_length, int collision_family) {
    if (!this->GetMesh())
        return;
    for (unsigned int i = 0; i < this->GetMesh()->GetNnodes(); ++i)
        if (auto mnodeFEA = std::dynamic_pointer_cast<ChNodeFEAxyz>(this->GetMesh()->GetNode(i)))
            this->AddNode_Cylinder(mnodeFEA, cyl_radius, cyl_length, collision_family);
        else if (auto mnodeFEArot = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(this->GetMesh()->GetNode(i)))
            this->AddNode_Cylinder(mnodeFEArot, cyl_radius, cyl_length, collision_family);
}


/// Add nodes of the mesh, belonging to the node_set, to this collision cloud
void ChContactSurfaceNodeCloud::AddFacesFromNodeSet(std::vector<std::shared_ptr<ChNodeFEAbase> >& node_set,
                                                    const double point_radius, int collision_family) {
    if (!this->GetMesh())
        return;
    for (unsigned int i = 0; i < node_set.size(); ++i)
        if (auto mnodeFEA = std::dynamic_pointer_cast<ChNodeFEAxyz>(node_set[i]))
            this->AddNode(mnodeFEA, point_radius, collision_family);
        else if (auto mnodeFEArot = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_set[i]))
            this->AddNode(mnodeFEArot, point_radius, collision_family);
}

void ChContactSurfaceNodeCloud::SurfaceSyncCollisionModels() {
    for (unsigned int j = 0; j < vnodes.size(); j++) {
        this->vnodes[j]->GetCollisionModel()->SyncPosition();
    }
    for (unsigned int j = 0; j < vnodes_rot.size(); j++) {
        this->vnodes_rot[j]->GetCollisionModel()->SyncPosition();
    }
}

void ChContactSurfaceNodeCloud::SurfaceAddCollisionModelsToSystem(ChSystem* msys) {
    assert(msys);
    SurfaceSyncCollisionModels();
    for (unsigned int j = 0; j < vnodes.size(); j++) {
        msys->GetCollisionSystem()->Add(this->vnodes[j]->GetCollisionModel());
    }
    for (unsigned int j = 0; j < vnodes_rot.size(); j++) {
        msys->GetCollisionSystem()->Add(this->vnodes_rot[j]->GetCollisionModel());
    }
}

void ChContactSurfaceNodeCloud::SurfaceRemoveCollisionModelsFromSystem(ChSystem* msys) {
    assert(msys);
    for (unsigned int j = 0; j < vnodes.size(); j++) {
        msys->GetCollisionSystem()->Remove(this->vnodes[j]->GetCollisionModel());
    }
    for (unsigned int j = 0; j < vnodes_rot.size(); j++) {
        msys->GetCollisionSystem()->Remove(this->vnodes_rot[j]->GetCollisionModel());
    }
}

}  // end namespace fea
}  // end namespace chrono
