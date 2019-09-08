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

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

bool include_plasticity = false;


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

    std::vector<double> distance_vect;
    std::vector<double> contact_force;
};




int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    bool use_NSC;
    double timestep;


    if (argc == 1)
    {
        use_NSC = true;
        timestep = 0.025;
        GetLog() << "Run with two arguments: \n - [NSC|SMC] choose which contact method; \n - [dt] choose timestep";
    }
    else
    {
        if (argc != 3)
            throw std::exception("demo_IP_FEA_torqueIGA called with wrong number of arguments (either zero or two).");

        std::string contact_method(argv[1]);
        if (!contact_method.compare("NSC"))
            use_NSC = true;
        else if (!contact_method.compare("SMC"))
            use_NSC = false;
        else
            throw std::exception("demo_IP_FEA_torqueIGA called with unknown contact method.");

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

    // Here set the inward-outward margins for collision shapes: should make sense in the scale of the model
    collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
    collision::ChCollisionModel::SetDefaultSuggestedMargin(0.002);
    collision::ChCollisionSystemBullet::SetContactBreakingThreshold(0.0001);

    // Create a section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.

    double wire_diameter = 0.010;
    auto melasticity = std::make_shared<ChElasticityCosseratSimple>();
    melasticity->SetYoungModulus(0.5e9);
    melasticity->SetGshearModulus(0.5e9 * 0.7);

    auto mdamping = std::make_shared<ChDampingCosseratLinear>();
    mdamping->SetDampingCoefficientsRe((1e-3)*ChVector<>(1, 1, 1));
    mdamping->SetDampingCoefficientsRk((1e-4)*ChVector<>(1, 1, 1));

    auto mplasticity = std::make_shared<ChPlasticityCosseratLumped>();
    mplasticity->n_yeld_Mx = std::make_shared<ChFunction_Ramp>(1, 0.01);
    mplasticity->n_yeld_My = std::make_shared<ChFunction_Ramp>(0.2, 0.001);
    mplasticity->n_yeld_Mz = std::make_shared<ChFunction_Ramp>(0.2, 0.001);

    auto msection = std::make_shared<ChBeamSectionCosserat>(melasticity, include_plasticity ? mplasticity : nullptr, mdamping);
    msection->SetDensity(1000);
    msection->SetAsCircularSection(wire_diameter);

    if (use_NSC)
    {
        auto mysurfmaterial_NSC = std::dynamic_pointer_cast<ChMaterialSurfaceNSC>(mysurfmaterial);
        mysurfmaterial_NSC->SetRestitution(0.1f);
        mysurfmaterial_NSC->SetFriction(0.2f);
        mysurfmaterial_NSC->SetCompliance(0.0000005f);
        mysurfmaterial_NSC->SetComplianceT(0.0000005f);
        mysurfmaterial_NSC->SetDampingF(0.2f);
    }
    else
    {
        auto mysurfmaterial_SMC = std::dynamic_pointer_cast<ChMaterialSurfaceSMC>(mysurfmaterial);
        mysurfmaterial_SMC->SetYoungModulus(1.2e8);
        mysurfmaterial_SMC->SetFriction(0.3f);
        mysurfmaterial_SMC->SetRestitution(0.2f);
        mysurfmaterial_SMC->SetAdhesion(0);
    }
    

    
    ///////////////
    // BEAM BUNDLE
    ///////////////
    const int beam_num = 8;
    const int elements_foreachbeam = 56;
    const double bundle_radius = 0.025;
    const double bundle_length0 = 0.5;
    const double flange_thickness = 0.05;

    auto fixedAsset = std::make_shared<ChColorAsset>(1.0, 0.0, 0.0); // color for fixed objects
    auto movingAsset = std::make_shared<ChColorAsset>(0.0, 0.0, 1.0); // color for moving objects

    auto mFlangeStart = std::make_shared<ChBodyEasyCylinder>((bundle_radius + wire_diameter)*1.1, flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeStart->SetPos(-0.5*VECT_Z*flange_thickness);
    mFlangeStart->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
    mFlangeStart->SetBodyFixed(true);
    mFlangeStart->AddAsset(fixedAsset);
    my_system->Add(mFlangeStart);

    auto mFlangeEnd = std::make_shared<ChBodyEasyCylinder>((bundle_radius+ wire_diameter)*1.1, flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeEnd->SetPos(VECT_Z*(bundle_length0+0.5*flange_thickness));
    mFlangeEnd->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
    mFlangeEnd->AddAsset(movingAsset);
    my_system->Add(mFlangeEnd);

    auto mFlangeEndFixed = std::make_shared<ChBodyEasyBox>(2.5*(bundle_radius + wire_diameter), 2.5*(bundle_radius + wire_diameter), flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeEndFixed->SetPos(mFlangeEnd->GetPos() + VECT_Z* flange_thickness);
    mFlangeEndFixed->SetBodyFixed(true);
    mFlangeEndFixed->AddAsset(fixedAsset);
    my_system->Add(mFlangeEndFixed);


    auto central_cyl = std::make_shared<ChBodyEasyCylinder>(bundle_radius*0.5, bundle_length0, 1000, true, true, mat_surface_type);
    central_cyl->SetBodyFixed(true);
    central_cyl->SetPos(ChVector<>(VECT_Z)*bundle_length0 / 2.0);
    central_cyl->SetRot(Q_from_AngX(CH_C_PI_2));
    central_cyl->SetMaterialSurface(mysurfmaterial);
    central_cyl->AddAsset(fixedAsset);
    my_system->Add(central_cyl);

    auto rotMot = std::make_shared<ChLinkMotorRotationSpeed>();
    auto speedFun = std::make_shared<ChFunction_Const>(1);
    rotMot->SetSpeedFunction(speedFun);
    //rotMot->SetAngleOffset(30.0*CH_C_DEG_TO_RAD); // not working as expected
    rotMot->Initialize(mFlangeEndFixed, mFlangeEnd, ChFrame<>(mFlangeEnd->GetPos(), QUNIT));
    my_system->Add(rotMot);

    double contact_radius = wire_diameter/2.0;

    
    ChBuilderBeamIGA beam_build;
    ChBuilderBeamANCF builder;

    //std::array<std::shared_ptr<ChMesh>, beam_num> my_mesh_beams;
    //std::array<std::shared_ptr<ChContactSurfaceNodeCloud>, beam_num> mcontactcloud;

    std::array<std::shared_ptr<ChLinkMateGeneric>, beam_num> test_force_constraint;
    //std::array<std::shared_ptr<ChMesh>, beam_num> my_mesh_beams;

    // IGA beams between flanges
    for (auto beam_sel = 0; beam_sel<beam_num; ++beam_sel)
    {

        auto my_mesh_beams = std::make_shared<ChMesh>();
        // create beams at given position
        ChVector<> beam_pos_start(bundle_radius*cos(beam_sel*CH_C_2PI / beam_num), bundle_radius*sin(beam_sel*CH_C_2PI / beam_num), 0);
        ChVector<> beam_pos_end(bundle_radius*cos(beam_sel*CH_C_2PI / beam_num), bundle_radius*sin(beam_sel*CH_C_2PI / beam_num), bundle_length0);
        beam_build.BuildBeam(my_mesh_beams, msection, elements_foreachbeam, beam_pos_start, beam_pos_end, VECT_Y, 3);
        auto firstNode = beam_build.GetLastBeamNodes().front();
        auto endNode = beam_build.GetLastBeamNodes().back();

        // fix beam to start flange
        auto startConstr = std::make_shared<ChLinkMateGeneric>();
        startConstr->Initialize(firstNode, mFlangeStart, false, firstNode->Frame(), firstNode->Frame());
        startConstr->SetConstrainedCoords(true, true, true, true, true, true);
        my_system->Add(startConstr);

        test_force_constraint[beam_sel] = startConstr;

        // fix beam to end flange
        auto endConstr = std::make_shared<ChLinkMateGeneric>();
        endConstr->Initialize(endNode, mFlangeEnd, false, endNode->Frame(), endNode->Frame());
        endConstr->SetConstrainedCoords(true, true, true, true, true, true);
        my_system->Add(endConstr);

        auto mcontactcloud = std::make_shared<ChContactSurfaceNodeCloud>();
        my_mesh_beams->AddContactSurface(mcontactcloud);
        mcontactcloud->AddAllNodes(contact_radius);  // use larger point size to match beam section radius
        mcontactcloud->SetMaterialSurface(mysurfmaterial);
        my_system->Add(my_mesh_beams);

        my_mesh_beams->SetAutomaticGravity(false);


        // VISUALIZATION
        // beams
        auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
        mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
        mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
        mvisualizebeamA->SetSmoothFaces(true);
        mvisualizebeamA->SetWireframe(false);
        my_mesh_beams->AddAsset(mvisualizebeamA);

        auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
        mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NONE);
        mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizebeamC->SetSymbolsThickness(0.006);
        mvisualizebeamC->SetSymbolsScale(0.01);
        mvisualizebeamC->SetZbufferHide(false);
        my_mesh_beams->AddAsset(mvisualizebeamC);
    }



    // Irrlicht application
    ChIrrApp application(my_system.get(), L"Torquing IGA beams", core::dimension2d<u32>(1200, 900), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    //application.AddTypicalLogo();
    application.AddTypicalSky();
    //application.AddTypicalLights();
    application.AddLightWithShadow(100 * irr::core::vector3df(0.50867, 0.306668, 0.264132), 100 * irr::core::vector3df(0.178601, 0.081, 0.275477), 150, 50, 150, 90);
    //auto targ = core::vector3df(0.172294, 0.180592, 0.038451 | -0.0228889, 4.4704e-10, 0.177341);
    application.AddTypicalCamera(irr::core::vector3df(0.172294, 0.180592, 0.038451), irr::core::vector3df(-0.0228889, 4.4704e-10, 0.177341));
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
    application.SetContactsDrawMode(ChIrrTools::eCh_ContactsDrawMode::CONTACT_NORMALS);

    ChTimer<> tim;

    GetLog() << "demo_IP_FEA_torqueIGA; Contact method: " << (use_NSC ? "NSC" : "SMC") << "; dt=" << timestep << "s \n\n";

    _label_reporter_class reporter;

    application.SetTimestep(timestep);
    application.SetVideoframeSaveInterval(1);
    application.SetVideoframeSave(false);
    //application.SetPaused(true);
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        tim.start();
        application.DoStep();
        tim.stop();


        if (!application.GetPaused()) {

            double reactForceTot = 0;
            for (auto link_sel = 0; link_sel < test_force_constraint.size(); ++link_sel)
            {
                reactForceTot += test_force_constraint[link_sel]->Get_react_force().Length();
            }

            std::cout << std::setprecision(6)
                << "RotAngle: " << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << "; "
                << "Link Force: " << reactForceTot << "; ";
                //<< std::endl;

            my_system->GetContactContainer()->ReportAllContacts(&reporter);
            auto contact_force_max_ptr = std::max_element(reporter.contact_force.begin(), reporter.contact_force.end());
            auto contact_force_max = contact_force_max_ptr != reporter.contact_force.end() ? *contact_force_max_ptr : 0;
            std::cout << "Contact Force Max: " << contact_force_max << "; ";
            auto distance_min_ptr = std::min_element(reporter.distance_vect.begin(), reporter.distance_vect.end());
            auto distance_min = distance_min_ptr != reporter.distance_vect.end() ? *distance_min_ptr : 1;
            std::cout << "Contact Distance Min: " << distance_min << "; ";

            std::cout << std::endl;



            csv_problem << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << reactForceTot << contact_force_max << distance_min << tim();
            csv_problem << std::endl;
        }

        //std::cout << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().X << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().Y << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getAbsolutePosition().Z << " | "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().X << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().Y << ", "
        //          << application.GetSceneManager()->getActiveCamera()->getTarget().Z << std::endl;


        if (rotMot->GetMotorRot()*CH_C_RAD_TO_DEG > 360.0 || rotMot->GetMotorRot() != rotMot->GetMotorRot())
        {
            filesystem::create_directory(filesystem::path("video_capture"));
            irr::video::IImage* image = application.GetVideoDriver()->createScreenShot();
            std::ostringstream outframe;
            outframe << "frame_" << (use_NSC ? "NSC" : "SMC") << std::setw(6) << std::setfill('0') << std::floor(timestep*1e6) << ".bmp";
            if (image)
                application.GetDevice()->getVideoDriver()->writeImageToFile(image, outframe.str().c_str());
            image->drop();

            break;
        }


        application.EndScene();
    }

    std::ostringstream outfile;
    outfile << "iplog_" << (use_NSC ? "NSC" : "SMC") << std::setw(6) << std::setfill('0') << std::floor(timestep*1e6) << ".txt";
    csv_problem.write_to_file(outfile.str(), "rotangle, linkforce, time\n");



    std::cout << "Elapsed time: " << tim() << "s" << std::endl;


    return 0;
}

