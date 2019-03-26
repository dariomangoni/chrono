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
// Authors: Dario Mangoni
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"

#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include <chrono_interiorpoint/ChInteriorPoint.h>
#include <functional>
#include "chrono_mumps/ChSolverMumps.h"
#include "physics/ChSystemSMC.h"

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;


class BuildTetrahedronBeam
{
private:

	// There are 13 different divisions of a cube in tetrahedrons
	// 12 in 6 tetra; 1 in 5 tetra
	// numbers [0-7] in the table refer to node indexes for a generic cube
	// 0  | 0137, 0147, 1237, 1267, 1457, 1567
	// 1  | 0137, 0147, 1235, 1457, 2356, 3567
	// 2  | 0137, 0147, 1237, 1257, 1457, 2567
	// 3  | 0137, 0147, 1236, 1367, 1457, 1567
	// 4  | 0137, 0157, 0457, 1237, 1257, 2567
	// 5  | 0137, 0157, 0457, 1236, 1367, 1567
	// 6  | 0134, 1237, 1257, 1347, 1457, 2567
	// 7  | 0134, 1236, 1347, 1367, 1457, 1567
	// 9  | 0137, 0147, 1235, 1457, 2357, 2567
	// 10 | 0137, 0157, 0457, 1235, 2357, 2567
	// 11 | 0134, 1235, 1347, 1457, 2356, 3567
	// 12 | 0134, 1235, 1347, 1457, 2357, 2567
	// 13 | 0134, 1236, 1346, 1456, 3467

	ChVector<int> num_elements = { 1,1,1 };
	ChVector<double> position_offset = { 0.0, 0.0, 0.0 };
	ChQuaternion<double> rotation_offset = { 1.0, 0.0, 0.0, 0.0 };
	double cube_edge = 1;

	static const std::array<ChVector<int>, 8> cube_nodes_position;
	static const std::array<std::array<int, 4>, 6> tetahedra_nodes_order;

	std::vector<std::shared_ptr<ChNodeFEAxyz>> nodes_list; // not needed if mesh returns a non-const vector of nodes...

	std::shared_ptr<ChNodeFEAxyz>& getNode(ChVector<int> element_position, int node_local_id)
	{
		auto num_nodes = num_elements + ChVector<int>{1, 1, 1};
		element_position += cube_nodes_position[node_local_id];
		return nodes_list[element_position.x() + element_position.y()*num_nodes.x() + element_position.z()*num_nodes.x()*num_nodes.y()];
	}

public:
	BuildTetrahedronBeam() {}

	BuildTetrahedronBeam(double cubes_edge_length, int num_elements_x, int num_elements_y, int num_elements_z)
	{
		SetBeamDimensions(cubes_edge_length, num_elements_x, num_elements_y, num_elements_z);
	}

	BuildTetrahedronBeam(std::shared_ptr<ChMesh> mesh, std::shared_ptr<ChContinuumElastic> tetra_material, double cubes_edge_length, int num_elements_x, int num_elements_y, int num_elements_z, ChVector<double> pos = ChVector<>(), ChQuaternion<double> rot = ChQuaternion<>())
	{
		SetBeamDimensions(cubes_edge_length, num_elements_x, num_elements_y, num_elements_z, pos);
		CreateBeam(mesh, tetra_material);
	}

	void SetBeamDimensions(double cubes_edge_length, int num_elements_x, int num_elements_y, int num_elements_z, ChVector<double> pos = ChVector<double>(), ChQuaternion<double> rot = ChQuaternion<>(1.0, 0.0, 0.0, 0.0))
	{
		assert(num_elements_x > 0 && num_elements_y > 0 && num_elements_z > 0 && "Wrong beam dimensions");
		num_elements = { num_elements_x, num_elements_y, num_elements_z };
		cube_edge = cubes_edge_length;
		position_offset = pos;
		rotation_offset = rot;
	}

	std::list<std::shared_ptr<ChNodeFEAbase>> GetNodes(std::function<bool(int, int, int)> choose_fnc) const
	{
		auto num_nodes = num_elements + ChVector<int>{1, 1, 1};
		std::list<std::shared_ptr<ChNodeFEAbase>> nodes_list_export;
		// create nodes
		for (auto node_sel_z = 0; node_sel_z< num_nodes.z(); ++node_sel_z)
		{
			for (auto node_sel_y = 0; node_sel_y < num_nodes.y(); ++node_sel_y)
			{
				for (auto node_sel_x = 0; node_sel_x < num_nodes.x(); ++node_sel_x)
				{
					if (choose_fnc(node_sel_x, node_sel_y, node_sel_z))
						nodes_list_export.push_back(nodes_list[node_sel_x + node_sel_y*num_nodes.x() + node_sel_z*num_nodes.x()*num_nodes.y()]);
				}
			}
		}

		return nodes_list_export;
	}


	void CreateBeam(std::shared_ptr<ChMesh> mesh, std::shared_ptr<ChContinuumElastic> tetra_material)
	{
		auto num_nodes = num_elements + ChVector<int>{1, 1, 1};

		// create nodes
		nodes_list.resize(num_nodes.x()*num_nodes.y()*num_nodes.z());
		for (auto node_sel_z = 0; node_sel_z< num_nodes.z(); ++node_sel_z)
		{
			for (auto node_sel_y = 0; node_sel_y < num_nodes.y(); ++node_sel_y)
			{
				for (auto node_sel_x = 0; node_sel_x < num_nodes.x(); ++node_sel_x)
				{
					nodes_list[node_sel_x + node_sel_y*num_nodes.x() + node_sel_z*num_nodes.x()*num_nodes.y()] = std::make_shared<ChNodeFEAxyz>(rotation_offset.Rotate(ChVector<double>(node_sel_x*cube_edge, node_sel_y*cube_edge, node_sel_z*cube_edge)) + position_offset);
					mesh->AddNode(nodes_list[node_sel_x + node_sel_y*num_nodes.x() + node_sel_z*num_nodes.x()*num_nodes.y()]);
				}
			}
		}


		// select the element
		ChVector<int> elem_pos;
		for (elem_pos.z() = 0; elem_pos.z() < num_elements.z(); ++elem_pos.z())
		{
			for (elem_pos.y() = 0; elem_pos.y() < num_elements.y(); ++elem_pos.y())
			{
				for (elem_pos.x() = 0; elem_pos.x() < num_elements.x(); ++elem_pos.x())
				{
					// associate nodes to tetrahedra
					for (auto tetra_sel = 0; tetra_sel < 6; ++tetra_sel)
					{
						auto tetra_element = std::make_shared<ChElementTetra_4>();
						tetra_element->SetNodes(getNode(elem_pos, tetahedra_nodes_order[tetra_sel][0]),
							getNode(elem_pos, tetahedra_nodes_order[tetra_sel][1]),
							getNode(elem_pos, tetahedra_nodes_order[tetra_sel][2]),
							getNode(elem_pos, tetahedra_nodes_order[tetra_sel][3]));
						tetra_element->SetMaterial(tetra_material);
						mesh->AddElement(tetra_element);
					}
				}
			}
		}
	}




};

const std::array<ChVector<int>, 8> BuildTetrahedronBeam::cube_nodes_position = { { { 0,0,1 },{ 1,0,1 },{ 1,0,0 },{ 0,0,0 },{ 0,1,1 },{ 1,1,1 },{ 1,1,0 },{ 0,1,0 } } }; // the 1st node is on the origin
const std::array<std::array<int, 4>, 6> BuildTetrahedronBeam::tetahedra_nodes_order = { { { 0,3,1,7 },{ 0,1,4,7 },{ 1,3,2,7 },{ 1,2,6,7 },{ 1,5,4,7 },{ 1,6,5,7 } } };

#define USE_NSC

int main(int argc, char* argv[]) {
#ifdef USE_NSC
	bool contact_model_NSC = true;
#else
	bool contact_model_NSC = false;
#endif

	std::cout << "Interior-Point test bench using " << (contact_model_NSC ? "NSC" : "SMC") << " contact model" << std::endl;

	double cubes_edge = 0.05;
	int cubes_x = 12;
	int cubes_y = 3;
	int cubes_z = 28;
	double support_cylinder_radius = 0.05;
	double clearance = 0.01;

	// Create a Chrono::Engine physical system
#ifdef USE_NSC
	ChSystemNSC my_system;
	typedef ChMaterialSurfaceNSC material_surface_class;
#else
	ChSystemSMC my_system;
	typedef ChMaterialSurfaceSMC material_surface_class;
#endif

	auto timestep = contact_model_NSC ? 0.075 : 0.005;
	auto mat_surface_type = contact_model_NSC ? ChMaterialSurface::NSC : ChMaterialSurface::SMC;


	// Irrlicht setup
	ChIrrApp application(&my_system, L"IP FEA contacts", core::dimension2d<u32>(1440, 1080), false, true, true, irr::video::EDT_OPENGL);
	application.AddTypicalLogo();
	application.AddTypicalSky();
	application.AddTypicalLights();
	application.AddTypicalCamera(core::vector3df(1.2f, 0.5f, 0.8f));
	application.AddLightWithShadow(core::vector3df(2.5, 3, 2.5), core::vector3df(0, 0, 0), 3, 1, 7.2, 40, 512, video::SColorf(1, 1, 1));


	//application.SetContactsDrawMode(ChIrrTools::CONTACT_NORMALS);
	//application.SetPlotLinkFrames(true);
	application.SetSymbolscale(0.25);

	// Physic system
	collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
	collision::ChCollisionModel::SetDefaultSuggestedMargin(0.001);

	auto mysurfmaterial = std::make_shared<material_surface_class>();
	mysurfmaterial->SetFriction(0.4f);
#ifdef USE_NSC
	mysurfmaterial->SetCompliance(0.0000005f);
	mysurfmaterial->SetComplianceT(0.0000005f);
	mysurfmaterial->SetDampingF(0.2f);
#else
	mysurfmaterial->SetYoungModulus(6e4);
	mysurfmaterial->SetFriction(0.3f);
	mysurfmaterial->SetRestitution(0.2f);
	mysurfmaterial->SetAdhesion(0);
#endif

	// FEA tetahedron beam
	auto tetra_mesh = std::make_shared<ChMesh>();

	auto tetra_material = std::make_shared<ChContinuumElastic>();
	tetra_material->Set_E(0.01e9);  // rubber 0.01e9, steel 200e9
	tetra_material->Set_v(0.3);
	tetra_material->Set_RayleighDampingK(0.003);
	tetra_material->Set_density(1000);

	BuildTetrahedronBeam(tetra_mesh, tetra_material, cubes_edge, cubes_x, cubes_y, cubes_z, ChVector<double>(-cubes_edge*static_cast<double>(cubes_x) / 2.0, support_cylinder_radius + clearance, -cubes_edge * static_cast<double>(cubes_z) / 2.0));

	auto mcontactsurf = std::make_shared<ChContactSurfaceMesh>();
	tetra_mesh->AddContactSurface(mcontactsurf);
	mcontactsurf->AddFacesFromBoundary(0.002);
	mcontactsurf->SetMaterialSurface(mysurfmaterial);

	my_system.Add(tetra_mesh);

	// Add cylinders

	auto cylinder_left = std::make_shared<ChBodyEasyCylinder>(support_cylinder_radius, 1, 1000, true, true, mat_surface_type);
	cylinder_left->SetMaterialSurface(mysurfmaterial);
	cylinder_left->SetPos(ChVector<>(0, 0, -0.5));
	cylinder_left->SetRot(Q_from_AngZ(-CH_C_PI_2));
	cylinder_left->SetBodyFixed(true);
	my_system.Add(cylinder_left);

	auto cylinder_right = std::make_shared<ChBodyEasyCylinder>(support_cylinder_radius, 1, 1000, true, true, mat_surface_type);
	cylinder_right->SetMaterialSurface(mysurfmaterial);
	cylinder_right->SetPos(ChVector<>(0, 0, 0.5));
	cylinder_right->SetRot(Q_from_AngZ(-CH_C_PI_2));
	cylinder_right->SetBodyFixed(true);
	my_system.Add(cylinder_right);

	auto cylinder_up = std::make_shared<ChBodyEasyCylinder>(support_cylinder_radius, 1, 1000, true, true, mat_surface_type);
	cylinder_up->SetMaterialSurface(mysurfmaterial);
	cylinder_up->SetPos(ChVector<>(0, 2 * (support_cylinder_radius + clearance) + cubes_y*cubes_edge, 0));
	cylinder_up->SetRot(Q_from_AngZ(-CH_C_PI_2));
	my_system.Add(cylinder_up);

	auto fixed_truss = std::make_shared<ChBodyEasyBox>(0.1, 0.1, 0.1, 1000, false, false, mat_surface_type);
	fixed_truss->SetBodyFixed(true);
	fixed_truss->SetPos(cylinder_up->GetPos() + ChVector<>(0.0, clearance, 0.0));
	my_system.Add(fixed_truss);

	auto cyl_motion = std::make_shared<ChLinkLockLock>();
	cyl_motion->Initialize(cylinder_up, fixed_truss, ChCoordsys<>(VNULL, Q_from_AngX(-CH_C_PI_2)));
	auto cyl_motion_function = std::make_shared<ChFunction_Ramp>();
	cyl_motion_function->Set_ang(-1);
	cyl_motion_function->Set_y0(0.0);
	cyl_motion->SetMotion_Z(cyl_motion_function);
	//ChLinkLimit link_lim;
	//link_lim.SetMax(-2.0);
	//cyl_motion->SetLimit_Z(&link_lim);
	my_system.AddLink(cyl_motion);

	// Visualization settings
	auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*tetra_mesh.get());
	mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_STRESS_VONMISES);
	mvisualizemesh->SetColorscaleMinMax(0.0, 2e6);
	mvisualizemesh->SetSmoothFaces(true);
	tetra_mesh->AddAsset(mvisualizemesh);

	auto mvisualizemeshcoll = std::make_shared<ChVisualizationFEAmesh>(*tetra_mesh.get());
	mvisualizemeshcoll->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_CONTACTSURFACES);
	mvisualizemeshcoll->SetWireframe(true);
	mvisualizemeshcoll->SetDefaultMeshColor(ChColor(1, 0.5, 0));
	tetra_mesh->AddAsset(mvisualizemeshcoll);

	// Application settings
	application.AssetBindAll();
	application.AssetUpdateAll();
	application.AddShadowAll();
	my_system.SetupInitial();

	// Timestepper and solver
	my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

#ifdef USE_NSC
	auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
	auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
	my_system.SetStabSolver(ip_solver_stab);
	my_system.SetSolver(ip_solver_speed);
	ip_solver_speed->RecordHistory(false);
	//ip_solver_speed->SetNullPivotDetection(true, 1e-18);
	ip_solver_speed->SetVerbose(true);
	ip_solver_speed->SetUseSymmetry(false);
	application.GetSystem()->Update();
#else
	my_system.SetSolverType(ChSolver::Type::MINRES);
	my_system.SetSolverWarmStarting(true);
	my_system.SetMaxItersSolverSpeed(200);
	my_system.SetTolForce(1e-10);
	my_system.GetSolver()->SetVerbose(true);
#endif


	// Run simulation
	//application.SetPaused(true);
	application.SetTimestep(timestep);

	ChTimer<> timer;
	timer.start();
	while (application.GetDevice()->run()) {
		application.BeginScene();

		application.DrawAll();

		application.DoStep();

		application.EndScene();

		//if (my_system.GetChTime() > 0.1)
		//	break;

	}
	timer.stop();
#ifdef USE_NSC
	std::cout << "Time spent: " << ip_solver_speed->GetTimeSolve_Assembly() + ip_solver_speed->GetTimeSolve_SolverCall() << "; IP calls: " << ip_solver_speed->GetIPSolverCalls() << "; IP iterations: " << ip_solver_speed->GetIPIterations() << std::endl;
	std::cout << "Mass matrix size: " << ip_solver_speed->GetMassMatrixDimension() << "; Inequalities: " << ip_solver_speed->GetUnilateralConstraintsMatrixRows() << "; Equalities: " << ip_solver_speed->GetBilateralConstraintsMatrixRows() << std::endl;
	std::cout << "Time spent (external timer): " << timer.GetTimeSeconds() << std::endl;
#else
	std::cout << "Time spent (external timer): " << timer.GetTimeSeconds() << std::endl;

#endif

	getchar();

	return 0;
}
