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
// Authors: Alessandro Tasora
// =============================================================================
//
// FEA for 3D beams
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/collision/ChCCollisionSystemBullet.h"

#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono_irrlicht/ChIrrApp.h"

#include "chrono_interiorpoint/ChInteriorPoint.h"
#include <iomanip>
#include "utils/ChUtilsInputOutput.h"
#include <chrono_thirdparty/filesystem/path.h>
#include <functional>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

bool include_plasticity = false;

double scaleLengthUnit = 1.0;

void drawReferenceFrame(std::shared_ptr<ChSystem> my_system)
{
    double refFrameScale = 0.1;
    double axesLength = 1.0*refFrameScale;
    double axesRadius = 0.1*refFrameScale;

    auto axisX = std::make_shared<ChBodyEasyCylinder>(axesRadius, axesLength, false, true);
    auto axisY = std::make_shared<ChBodyEasyCylinder>(*axisX);
    auto axisZ = std::make_shared<ChBodyEasyCylinder>(*axisX);
    axisX->SetRot(Q_from_AngAxis(-CH_C_PI_2, VECT_Z));
    axisX->SetPos(ChVector<>(+0.5*axesLength, 0.0, 0.0));

    axisY->SetRot(Q_from_AngAxis(0.0, VECT_Z));
    axisY->SetPos(ChVector<>(0.0, +0.5*axesLength, 0.0));

    axisZ->SetRot(Q_from_AngAxis(+CH_C_PI_2, VECT_X));
    axisZ->SetPos(ChVector<>(0.0, 0.0, +0.5*axesLength));

    axisX->SetBodyFixed(true);
    axisY->SetBodyFixed(true);
    axisZ->SetBodyFixed(true);
    axisX->AddAsset(std::make_shared<ChColorAsset>(1.0, 0.0, 0.0));
    axisY->AddAsset(std::make_shared<ChColorAsset>(0.0, 1.0, 0.0));
    axisZ->AddAsset(std::make_shared<ChColorAsset>(0.0, 0.0, 1.0));
    my_system->Add(axisX);
    my_system->Add(axisY);
    my_system->Add(axisZ);
}


class _label_reporter_class : public ChContactContainer::ReportContactCallback {
public:
    virtual bool OnReportContact(const ChVector<>& pA,
        const ChVector<>& pB,
        const ChMatrix33<>& plane_coord,
        const double& distance,
        const double& eff_radius,
        const ChVector<>& react_forces,
        const ChVector<>& react_torques,
        ChContactable* modA,
        ChContactable* modB) override {


        contact_force.push_back(react_forces.x());
        distance_vect.push_back(distance);
        //std::cout << "Contact force " << react_forces.x();

        return true;  // to continue scanning contacts
    }

    void Reset()
    {
        distance_vect.clear();
        contact_force.clear();
    }

    std::vector<double> distance_vect;
    std::vector<double> contact_force;

    bool isContact() const
    {
        return distance_vect.size() > 0;
    }
};




int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    bool use_NSC;
    double timestep;


    if (argc == 1)
    {
        use_NSC = true;
        timestep = 0.025;
        GetLog() << "Run with two arguments: \n - [NSC|SMC] choose which contact method; \n - [dt] choose timestep\n";
    }
    else
    {
        if (argc != 3)
            throw std::exception("demo_IP_FEA_netIGA called with wrong number of arguments (either zero or two).");

        std::string contact_method(argv[1]);
        if (!contact_method.compare("NSC"))
            use_NSC = true;
        else if (!contact_method.compare("SMC"))
            use_NSC = false;
        else
            throw std::exception("demo_IP_FEA_netIGA called with unknown contact method.");

        std::string timestep_string(argv[2]);
        timestep = std::stod(timestep_string);
    }


    std::shared_ptr<ChSystem> my_system;
    std::shared_ptr<ChMaterialSurface> mysurfmaterial;
    ChMaterialSurface::ContactMethod mat_surface_type;
    if (use_NSC)
    {
        my_system = std::make_shared<ChSystemNSC>();
        mysurfmaterial = std::make_shared<ChMaterialSurfaceNSC>();

        mat_surface_type = ChMaterialSurface::NSC;
    }
    else
    {
        my_system = std::make_shared<ChSystemSMC>();
        mysurfmaterial = std::make_shared<ChMaterialSurfaceSMC>();
        mat_surface_type = ChMaterialSurface::SMC;

    }



    // Create a section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.

    double wire_diameter = 0.002*scaleLengthUnit;
    double youngModulus = 1e8/scaleLengthUnit;
    auto melasticity = std::make_shared<ChElasticityCosseratSimple>();
    melasticity->SetYoungModulus(youngModulus);
    melasticity->SetGshearModulus(melasticity->GetYoungModulus() * 0.5);

    auto mdamping = std::make_shared<ChDampingCosseratLinear>();
    mdamping->SetDampingCoefficientsRe((1e-3)/scaleLengthUnit*ChVector<>(1, 1, 1));
    mdamping->SetDampingCoefficientsRk((1e-4)/scaleLengthUnit*ChVector<>(1, 1, 1));

    auto mplasticity = std::make_shared<ChPlasticityCosseratLumped>();
    mplasticity->n_yeld_Mx = std::make_shared<ChFunction_Ramp>(1*scaleLengthUnit, 0.01*scaleLengthUnit);
    mplasticity->n_yeld_My = std::make_shared<ChFunction_Ramp>(0.2*scaleLengthUnit, 0.001*scaleLengthUnit);
    mplasticity->n_yeld_Mz = std::make_shared<ChFunction_Ramp>(0.2*scaleLengthUnit, 0.001*scaleLengthUnit);

    auto msection = std::make_shared<ChBeamSectionCosserat>(melasticity, include_plasticity ? mplasticity : nullptr, mdamping);
    msection->SetDensity(1000.0/scaleLengthUnit/scaleLengthUnit/scaleLengthUnit);
    msection->SetAsCircularSection(wire_diameter);

    if (use_NSC)
    {
        auto mysurfmaterial_NSC = std::dynamic_pointer_cast<ChMaterialSurfaceNSC>(mysurfmaterial);
        mysurfmaterial_NSC->SetRestitution(0.1f);
        mysurfmaterial_NSC->SetFriction(0.2f);
        mysurfmaterial_NSC->SetCompliance(0.0000005f*scaleLengthUnit);
        mysurfmaterial_NSC->SetComplianceT(0.0000005f*scaleLengthUnit);
        mysurfmaterial_NSC->SetDampingF(0.2f);
    }
    else
    {
        auto mysurfmaterial_SMC = std::dynamic_pointer_cast<ChMaterialSurfaceSMC>(mysurfmaterial);
        mysurfmaterial_SMC->SetYoungModulus(youngModulus);
        mysurfmaterial_SMC->SetFriction(0.3f);
        mysurfmaterial_SMC->SetRestitution(0.2f);
        mysurfmaterial_SMC->SetAdhesion(0);
    }
    

    
    ///////////////
    // BEAM BUNDLE
    ///////////////
    double wavelength = 0.03*scaleLengthUnit;
    double contact_radius = wire_diameter/2.0;

    int spline_order = 3;
    const bool wave_beams = true;
    const int beam_num = 7;
    const int net_layers = 2;
    const int elements_foreachbeam = 64;
    //const double grid_spacing = wavelength/2.0;
    const double beam_length0 = (beam_num+1)/2.0*wavelength;
    double start_position = -0.5*beam_length0;
    double amplitude = 1.5*0.5*wire_diameter;
    //const double layer_spacing = 0;
    double omega = CH_C_2PI/wavelength;

    auto fixedAsset = std::make_shared<ChColorAsset>(1.0, 0.0, 0.0, 0.5); // color for fixed objects
    auto movingAsset = std::make_shared<ChColorAsset>(1.0, 1.0, 1.0); // color for moving objects

    //drawReferenceFrame(my_system);


    auto fixedBody = std::make_shared<ChBody>(mat_surface_type);
    fixedBody->SetBodyFixed(true);
    my_system->Add(fixedBody);


    //double supportThickness = 0.01;

    //auto supportAlongZ = std::make_shared<ChBodyEasyBox>(supportThickness, supportThickness, beam_length0, false, true, mat_surface_type);
    //auto supportAlongX = std::make_shared<ChBodyEasyBox>(*supportAlongZ);
    //supportAlongZ->SetPos(ChVector<>(-0.5*supportThickness-0.5*beam_length0, -0.5*supportThickness, 0.0));
    //supportAlongZ->SetBodyFixed(true);
    //supportAlongZ->AddAsset(fixedAsset);
    //my_system->Add(supportAlongZ);

    //if (net_layers>1)
    //{
    //    supportAlongX->SetPos(ChVector<>(0.0, -0.5*supportThickness, -0.5*supportThickness-0.5*beam_length0));
    //    supportAlongX->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_Y));
    //    supportAlongX->SetBodyFixed(true);
    //    supportAlongX->AddAsset(fixedAsset);
    //    my_system->Add(supportAlongX);
    //}

    
    ChBuilderBeamIGA beam_build;
    ChBuilderBeamANCF builder;

    //std::array<std::shared_ptr<ChMesh>, beam_num> my_mesh_beams;
    //std::array<std::shared_ptr<ChContactSurfaceNodeCloud>, beam_num> mcontactcloud;

    //std::array<std::shared_ptr<ChLinkMateGeneric>, beam_num_x> test_force_constraint;
    //std::array<std::shared_ptr<ChMesh>, beam_num> my_mesh_beams;


    //my_system->Set_G_acc(-VECT_Y*9.81*scaleLengthUnit);
    my_system->Set_G_acc(0.0);

        // Here set the inward-outward margins for collision shapes: should make sense in the scale of the model
    //collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0005);
    //collision::ChCollisionModel::SetDefaultSuggestedMargin(0.001);
    //collision::ChCollisionSystemBullet::SetContactBreakingThreshold(0.0001);

    collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.2*contact_radius);
    collision::ChCollisionModel::SetDefaultSuggestedMargin(0.2*contact_radius);
    collision::ChCollisionSystemBullet::SetContactBreakingThreshold(1e-4);
    my_system->SetMaxPenetrationRecoverySpeed(0.5*contact_radius);

    double loadFactor = 1.0;
    double forceField = 0.1;

    std::shared_ptr<ChNodeFEAxyzrot> refNode;

    // IGA beams between flanges
    for (auto layer_sel = 0; layer_sel<net_layers; ++layer_sel)
    {    
        //if (layer_sel == 1 )
        //    continue;

        for (int beam_sel = 0; beam_sel<beam_num; ++beam_sel)
        {    

            //if (beam_sel > 4 )
            //    continue;
            
            double offset_phase = (layer_sel+beam_sel)*CH_C_PI;


            ChMatrix33<> rotMat;
            rotMat.Set_A_AngAxis(-layer_sel*CH_C_PI_2, VECT_Y);

            auto my_mesh_beams = std::make_shared<ChMesh>();
            // create beams at given position;
            // position before rotation: along Z, offset along Z
            //double z_position = (beam_sel+1)*grid_spacing-0.5*(beam_num+1)*grid_spacing;
            //double z_position = -0.5*beam_length0 + (beam_sel-beam_num+1)*beam_length0;
            double z_position = (static_cast<double>(beam_sel)-(static_cast<double>(beam_num) - 1.0)*0.5) * beam_length0 / (static_cast<double>(beam_num)+1.0);
            ChVector<> beam_pos_start(-0.5*beam_length0, pow(-1, beam_sel+layer_sel)*amplitude, z_position);
            ChVector<> beam_pos_end(+0.5*beam_length0, pow(-1, beam_sel+layer_sel)*amplitude, z_position);

            //ChVector<> beam_pos_start(0, 0, 0);
            //ChVector<> beam_pos_end(1, 0, 0);

            beam_pos_start = rotMat*beam_pos_start;
            beam_pos_end = rotMat*beam_pos_end;

            beam_build.BuildBeam(my_mesh_beams, msection, elements_foreachbeam, beam_pos_start, beam_pos_end, VECT_Y, spline_order);
            auto firstNode = beam_build.GetLastBeamNodes().front();
            auto endNode = beam_build.GetLastBeamNodes().back();

            double steplength;




            //firstNode->SetFixed(true);
            //endNode->SetFixed(true);



            if (wave_beams)
            {
                // MAKE WAVES
                ChMatrix33<> invRotMat(rotMat);
                invRotMat.MatrTranspose();

                auto& nodes = beam_build.GetLastBeamNodes();
                steplength = beam_length0/(nodes.size()-1);
                steplength *= 1.01;


                //std::cout << "Beam " << beam_sel << "\n";
                for (auto node_sel = 1; node_sel<nodes.size(); ++node_sel)
                {
                    auto& prevFrameGlob = nodes[node_sel-1]->Frame();
                    auto prevPosLoc = invRotMat.Matr_x_Vect(prevFrameGlob.GetPos());

                    // derivative step
                    ChVector<> derUnitaryLoc = ChVector<>(1.0, -amplitude*omega*sin(omega*(prevPosLoc.x()+steplength*0.5-start_position)+offset_phase), 0.0).GetNormalized();
                    double beamAngle = atan2(derUnitaryLoc.y(), derUnitaryLoc.x());
                    ChVector<> derUnitaryGlob = rotMat*derUnitaryLoc;
                    ChVector<> newPos = prevFrameGlob.GetPos() + steplength*derUnitaryGlob;

         
                    //// secant step
                    //double x_old = prevPosLoc.x();
                    //double x_oldold;
                    //if (node_sel>1)
                    //{
                    //    x_oldold = invRotMat.Matr_x_Vect(nodes[node_sel-2]->Frame().GetPos()).x();
                    //} else
                    //{
                    //    x_oldold = x_old - 0.01;
                    //}

                    //auto fun = [&](double pos_x){ return cos(omega*(pos_x-start_position)+offset_phase); };
                    //double x_new = secant_method(fun, x_old, x_oldold);
                    //double beamAngle = atan2(fun(x_new)- fun(x_old), x_new - x_old);
                    //ChVector<> newStep = rotMat*
                    //ChVector<> newPos = prevFrameGlob.GetPos() + 


                    ChVector<> rotAxis = rotMat*VECT_Z;
                    auto& prevQuat = nodes[node_sel]->Frame().GetRot();

                    ChQuaternion<> changeQuat = Q_from_AngAxis(beamAngle, rotAxis);

                    ChQuaternion<> newQuat;
                    newQuat.Cross(changeQuat, prevQuat);

                    ChVector<double> newAxis;
                    double newAngle;
                    Q_to_AngAxis(newQuat, newAngle, newAxis);

                    //ChVector<double> oldAxis;
                    //double oldAngle;
                    //Q_to_AngAxis(prevQuat, oldAngle, oldAxis);

                    //std::cout << "\n OLD Angle: "<< oldAngle << " Axis: ["<< oldAxis.x() << ", "<< oldAxis.y()<< ", "<< oldAxis.z() << "]";
                    //std::cout << "\n CNG Angle: "<< beamAngle << " Axis: ["<< rotAxis.x() << ", "<< rotAxis.y()<< ", "<< rotAxis.z() << "]";
                    //std::cout << "\n NEW Angle: "<< newAngle << " Axis: ["<< newAxis.x() << ", "<< newAxis.y()<< ", "<< newAxis.z() << "]";


                    nodes[node_sel]->Frame() = ChFrame<>(newPos, newQuat);

                    nodes[node_sel]->SetX0(ChFrame<>(newPos, newQuat));

                    nodes[node_sel]->SetForce(-VECT_Y*steplength*forceField*loadFactor);


                    //GetLog() << "\n Angle: "<< nodes[node_sel]->Frame().GetRotAngle() << "\n Axis: "<< nodes[node_sel]->Frame().GetRotAxis();
                    //std::cout << "\n";

                    refNode = nodes[node_sel];
                }
                //std::cout << "\n";
                //std::cout << "START: " << nodes[0]->Frame().GetPos().x() << ", " << nodes[0]->Frame().GetPos().y() << ", " << nodes[0]->Frame().GetPos().z() << "\n";
                //std::cout << "END:   " << nodes[nodes.size()-1]->Frame().GetPos().x() << ", " << nodes[nodes.size()-1]->Frame().GetPos().y() << ", " << nodes[nodes.size()-1]->Frame().GetPos().z() << "\n";

                //std::cout << "Axis(:, :, " << beam_sel+1 << ") = [" ;
                //for (auto node_sel = 0; node_sel<nodes.size(); ++node_sel)
                //{
                //    auto quat = nodes[node_sel]->Frame().GetRot();
                //    ChVector<double> newAxis;
                //    double newAngle;
                //    Q_to_AngAxis(quat, newAngle, newAxis);

                //     std:: cout << newAxis.x() << ", "<< newAxis.y()<< ", "<< newAxis.z() << ";\n";
                //}
                //std::cout << "]';\n" ;
                //std::cout << std::flush;

                //std::cout << "Angle(:, " << beam_sel+1 << ") = [" ;
                //for (auto node_sel = 0; node_sel<nodes.size(); ++node_sel)
                //{
                //    auto quat = nodes[node_sel]->Frame().GetRot();
                //    ChVector<double> newAxis;
                //    double newAngle;
                //    Q_to_AngAxis(quat, newAngle, newAxis);

                //     std:: cout << newAngle << ", ";
                //}
                //std::cout << "]';\n" ;
                //std::cout << std::flush;

                //std::cout << "X0(:,:," << beam_sel+1 << ") = [" ;
                //for (auto node_sel = 0; node_sel<nodes.size(); ++node_sel)
                //{
                //     std:: cout << nodes[node_sel]->GetX0().GetPos().x() << ", "
                //                << nodes[node_sel]->GetX0().GetPos().y() << ", "
                //                << nodes[node_sel]->GetX0().GetPos().z() << ";\n ";
                //}
                //std::cout << "]';\n" ;

                //std::cout << "X(:, :, " << beam_sel+1 <<", "<< layer_sel+1 << ") = [" ;
                //for (auto node_sel = 0; node_sel<nodes.size(); ++node_sel)
                //{
                //     std:: cout << nodes[node_sel]->Frame().GetPos().x() << ", "
                //                << nodes[node_sel]->Frame().GetPos().y() << ", "
                //                << nodes[node_sel]->Frame().GetPos().z() << ";\n ";
                //}
                //std::cout << "]';\n" ;

            }
            
            

            // Create contact cloud
            auto mcontactcloud = std::make_shared<ChContactSurfaceNodeCloud>();
            my_mesh_beams->AddContactSurface(mcontactcloud);
            //mcontactcloud->AddAllNodes_Spheres(contact_radius, layer_sel);  // use larger point size to match beam section radius
            mcontactcloud->AddAllNodes_Spheres(contact_radius, layer_sel+1);  // use larger point size to match beam section radius
            //mcontactcloud->AddAllNodes_Cylinders(cyl_radius, cyl_length, layer_sel+1);  // use larger point size to match beam section radius
            mcontactcloud->SetMaterialSurface(mysurfmaterial);

            my_mesh_beams->SetAutomaticGravity(false);

            // fix beam to start flange
            auto startConstr = std::make_shared<ChLinkMateGeneric>();
            startConstr->Initialize(firstNode, fixedBody, false, firstNode->Frame(), firstNode->Frame());
            startConstr->SetConstrainedCoords(true, true, true, true, true, true);
            my_system->Add(startConstr);

            my_system->Add(my_mesh_beams);


            // VISUALIZATION
            // beams
            auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
            mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
            mvisualizebeamA->SetColorscaleMinMax(-100, +100);
            mvisualizebeamA->SetSmoothFaces(true);
            mvisualizebeamA->SetWireframe(false);
            mvisualizebeamA->SetDefaultMeshColor(ChColor(1.0,1.0,1.0));
            my_mesh_beams->AddAsset(mvisualizebeamA);

            auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
            mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NONE);
            mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
            mvisualizebeamC->SetSymbolsThickness(0.006);
            mvisualizebeamC->SetSymbolsScale(0.01);
            mvisualizebeamC->SetZbufferHide(false);
            my_mesh_beams->AddAsset(mvisualizebeamC);
        }
    }




    // Irrlicht application
    ChIrrApp application(my_system.get(), L"IGA Net", core::dimension2d<u32>(1200, 900), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    //application.AddTypicalLogo();
    application.AddTypicalSky();
    //application.AddTypicalLights();
    application.AddLightWithShadow(100 * irr::core::vector3df(0.50867, 0.306668, 0.264132)*scaleLengthUnit, 100 * irr::core::vector3df(0.178601, 0.081, 0.275477)*scaleLengthUnit, 150, 50, 150, 90);
    //auto targ = core::vector3df(0.172294, 0.180592, 0.038451 | -0.0228889, 4.4704e-10, 0.177341);
    //application.AddTypicalCamera(irr::core::vector3df(0.172294, 0.180592, 0.038451), irr::core::vector3df(-0.0228889, 4.4704e-10, 0.177341));
    application.AddTypicalCamera(irr::core::vector3df(0.141734, 0.0488902, -0.00461925)*scaleLengthUnit, irr::core::vector3df(0.0,0.0,0.0)*scaleLengthUnit);
    //application.AddTypicalLights(targ + core::vector3df(+0.1f, +0.2f, -0.2f), targ);

    application.AddShadowAll();

    // Initial setup
    application.AssetBindAll();
    application.AssetUpdateAll();
    my_system->SetupInitial();



    // Solver and timestepper
    my_system->SetSolverType(ChSolver::Type::MINRES);
    my_system->SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system->SetMaxItersSolverSpeed(1000);
    my_system->SetMaxItersSolverStab(1000);
    my_system->SetTolForce(1e-8);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system->GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);

if (use_NSC)
{
    auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
    auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
    my_system->SetStabSolver(ip_solver_stab);
    my_system->SetSolver(ip_solver_speed);
    // ip_solver_speed->SetVerbose(true);
    application.GetSystem()->Update();
}
else
{
    //my_system->SetTimestepperType(ChTimestepper::Type::HHT);
    //auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(my_system->GetTimestepper());
    //mystepper->SetAlpha(-0.2);
    //mystepper->SetMaxiters(10);
    //mystepper->SetStepControl(false);
    //mystepper->SetAbsTolerances(1e-5);
    //mystepper->SetScaling(true);
    //mystepper->SetVerbose(true);
}

    
    utils::CSV_writer csv_problem;
    //application.SetPlotLinkFrames(true);
    //application.SetContactsDrawMode(ChIrrTools::eCh_ContactsDrawMode::CONTACT_NORMALS);
    application.SetPlotCollisionShapes(false);

    ChTimer<> tim;

    GetLog() << "demo_IP_FEA_netIGA; Contact method: " << (use_NSC ? "NSC" : "SMC") << "; dt=" << timestep << "s \n\n";

    _label_reporter_class reporter;

    application.SetTimestep(timestep);
    application.SetVideoframeSaveInterval(1);
    application.SetVideoframeSave(false);
    application.SetPaused(true);

    int step_count = 0;

    double old_refNodePosY = 1000;

    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        tim.start();
        application.DoStep();
        tim.stop();


        if (!application.GetPaused()) {

            step_count++;

            //double reactForceTot = 0;
            //for (auto link_sel = 0; link_sel < test_force_constraint.size(); ++link_sel)
            //{
            //    reactForceTot += test_force_constraint[link_sel]->Get_react_force().Length();
            //}

            //std::cout << std::setprecision(6)
            //    << "RotAngle: " << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << "; "
            //    << "Link Force: " << reactForceTot << "; ";
            //    //<< std::endl;

            std::cout << "Step: " << step_count << std::endl;
            //my_system->GetContactContainer()->ReportAllContacts(&reporter);
            //if (reporter.isContact())
            //{
            //    auto contact_force_max_ptr = std::max_element(reporter.contact_force.begin(), reporter.contact_force.end());
            //    auto contact_force_max = contact_force_max_ptr != reporter.contact_force.end() ? *contact_force_max_ptr : 0;
            //    std::cout << "Contact | Force Max: " << contact_force_max << " | ";
            //    auto distance_min_ptr = std::min_element(reporter.distance_vect.begin(), reporter.distance_vect.end());
            //    auto distance_min = distance_min_ptr != reporter.distance_vect.end() ? *distance_min_ptr : 1;
            //    std::cout << "Distance Min: " << distance_min << "; ";

            //    std::cout << std::endl;
            //}

       

            //reporter.Reset();



            //csv_problem << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << reactForceTot << contact_force_max << distance_min << tim();
            //csv_problem << std::endl;
        }

        //std::cout << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().X << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().Y << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().Z << " | "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().X << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().Y << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().Z << std::endl;

        std::cout << "Pos refNode y: " << refNode->GetPos().y() << "; delta: " << std::fabs(old_refNodePosY-(refNode->GetPos().y())) << std::endl;

        //if (std::fabs(old_refNodePosY-refNode->GetPos().y())<1e-7)
        //    application.SetPaused(true);

        old_refNodePosY = refNode->GetPos().y();

        application.EndScene();
    }

    std::ostringstream outfile;
    outfile << "iplog_" << (use_NSC ? "NSC" : "SMC") << std::setw(6) << std::setfill('0') << std::floor(timestep*1e6) << ".txt";
    csv_problem.write_to_file(outfile.str(), "rotangle, linkforce, time\n");



    std::cout << "Elapsed time: " << tim() << "s" << std::endl;


    return 0;
}

