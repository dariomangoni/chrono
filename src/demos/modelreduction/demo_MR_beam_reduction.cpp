//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be 
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

//
//   Demo code about 
//
//     - FEA for 3D beams



// Include some headers used by this tutorial...

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChLinkDirFrame.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_modelreduction/ChEigenAnalysis.h"
#include <typeinfo>

#include <unsupported/Eigen/SparseExtra>

#include <set>

//#include "chrono_matlab/ChMatlabEngine.h"
//#include "chrono_matlab/ChSolverMatlab.h"

// Remember to use the namespace 'chrono' because all classes 
// of Chrono::Engine belong to this namespace and its children...



using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

bool enable_eigenanalysis = true;
bool use_explicit_constraints = false;
bool use_realtime = false;
bool save_matrix = false;

int main(int argc, char* argv[])
{
    // Create a Chrono::Engine physical system
    ChSystemSMC my_system;


    // Create the Irrlicht visualization (open the Irrlicht device, 
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Beam modes", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(-0.1f, 0.2f, -0.2f));



    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();


    // Create a section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.

    double beam_wy = 0.04;                  //[m]
    double beam_wz = 0.05;                  //[m]
    double beam_L = 4;                      //[m]
    double E = 21e+10;                    //[Pa]
    double G = 6.3e+10;                     //[Pa]
    double density = 7800;                  //[kg/m^3]

    double Area = beam_wy * beam_wz;
    double Izz = (1.0 / 12.0) * beam_wz * pow(beam_wy, 3);
    double Iyy = (1.0 / 12.0) * beam_wy * pow(beam_wz, 3);
    double Ixx = Iyy + Izz;

    // use Roark's formulas for torsion of rectangular sect:
    double t = ChMin(beam_wy, beam_wz);
    double b = ChMax(beam_wy, beam_wz);
    double J = b * pow(t, 3) * ((1.0 / 3.0) - 0.210 * (t / b) * (1.0 - (1.0 / 12.0) * pow((t / b), 4)));


    auto msection = std::make_shared<ChBeamSectionEulerAdvancedGeneric>();

    msection->SetBeamRaleyghDamping(0.000);

    msection->SetAxialRigidity(Area*E);
    msection->SetXtorsionRigidity(J*G);
    msection->SetYbendingRigidity(Iyy*E);
    msection->SetZbendingRigidity(Izz*E);
    msection->SetMassPerUnitLength(Area*density);
    msection->SetInertiaJxxPerUnitLength(Ixx*density);

    //msection->SetArtificialJyyJzzFactor(1.0/500.0);
    ////msection->SetSectionRotation(45*CH_C_RAD_TO_DEG);


    //
    // Add some EULER-BERNOULLI BEAMS (the fast way!)
    //


    int num_elements = 32;

    ChBuilderBeamEuler builder;

    builder.BuildBeam(my_mesh,		// the mesh where to put the created nodes and elements 
                      msection,		// the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      num_elements,				// the number of ChElementBeamEuler to create
                      ChVector<>(0, 0, 0),		// the 'A' point in space (beginning of beam)
                      ChVector<>(beam_L, 0, 0),	// the 'B' point in space (end of beam)
                      ChVector<>(0, 1, 0));			// the 'Y' up direction of the section for the beam

                                                    // After having used BuildBeam(), you can retrieve the nodes used for the beam,
                                                    // For example say you want to fix the A end and apply a force to the B end:

    // DeadBody object
    auto my_deadbody = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_deadbody);
    my_deadbody->SetBodyFixed(true);
    my_deadbody->SetName("DeadBody");

    // Cube (rigid body)
    //auto my_cube = chrono_types::make_shared<ChBodyEasyBox>(0.1, 0.2, 0.3, 2700);
    //ChVector<double> mpos = (1, 0.5, 0.5);
    //my_cube->SetPos(mpos);
    //my_system.AddBody(my_cube);

    if(!use_explicit_constraints)
        builder.GetLastBeamNodes()[0]->SetFixed(true);
    else{
        //Fix link between beam and deadboy
        auto my_FixLink = chrono_types::make_shared<ChLinkMateFix>();
        my_FixLink->SetName("FixLink");
        my_FixLink->Initialize(my_deadbody, builder.GetLastBeamNodes()[0]);
        my_system.AddLink(my_FixLink);

    }

    builder.GetLastBeamNodes().back()->SetForce(ChVector<>(0, -1.0, 0));
    
    // Modes only on the vertical plane
    for (int i = 1; i < builder.GetLastBeamNodes().size()-1; i++) {
        auto abearing = chrono_types::make_shared<ChLinkMateGeneric>(true, false, true, true, true, false);
        abearing->Initialize(builder.GetLastBeamNodes()[i], my_deadbody, ChFrame<>(builder.GetLastBeamNodes()[i]->GetPos()));
        my_system.Add(abearing);
    }
    

    // Final touches..
    // 


    // We do not want gravity effect on FEA elements in this demo
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);




    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChTriangleMeshShape
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    // postprocessor that can handle a coloured ChTriangleMeshShape).
    // Do not forget AddAsset() at the end!


    /*
    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizebeamA->SetSmoothFaces(true);
    my_mesh->AddAsset(mvisualizebeamA);
    */

    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    //mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
    //mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizebeamA->SetDefaultMeshColor(ChColor(0,1,0,0));
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddAsset(mvisualizebeamA);

    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddAsset(mvisualizebeamC);


    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
    // in the system. These ChIrrNodeAsset assets are 'proxies' to the Irrlicht meshes.
    // If you need a finer control on which item really needs a visualization proxy in 
    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.

    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!

    application.AssetUpdateAll();

    //// Mark completion of system construction
    //my_system.SetupInitial();


    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(600);
    solver->SetTolerance(1e-13);
    solver->EnableDiagonalPreconditioner(true);
    //solver->SetVerbose(true);

    application.SetTimestep(0.01);

    application.DoStep();

    // Save matrix
    if (save_matrix) {
    
        if (!use_explicit_constraints) {
            ChSparseMatrix matM;
            application.GetSystem()->GetMassMatrix(&matM);
            matM.makeCompressed();
            //std::cout << matM << std::endl;
            Eigen::saveMarket(matM, "matM.mat");

            ChSparseMatrix matK;
            application.GetSystem()->GetStiffnessMatrix(&matK);
            matK.makeCompressed();
            //std::cout << matK << std::endl;
            Eigen::saveMarket(matK, "matK.mat");

            ChSparseMatrix matCq;
            application.GetSystem()->GetConstraintJacobianMatrix(&matCq);
            matCq.makeCompressed();
            //std::cout << matCq << std::endl;
            Eigen::saveMarket(matCq, "matCq.mat");
        
        }

        else {

            ChSparseMatrix matKaug;
            application.GetSystem()->KRMmatricesLoad(1.0, 0, 0);
            application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
            application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matKaug, nullptr, true);
            //std::cout << matKaug << std::endl;
            Eigen::saveMarket(matKaug, "matKaug.mat");

            ChSparseMatrix matMaug;
            application.GetSystem()->KRMmatricesLoad(0, 0, 1.0);
            application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
            application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matMaug, nullptr, false);
            //std::cout << matMaug << std::endl;
            Eigen::saveMarket(matMaug, "matMaug.mat");

            ChSparseMatrix matM_expl;
            application.GetSystem()->GetMassMatrix(&matM_expl);
            matM_expl.makeCompressed();
            //std::cout << matM_expl << std::endl;
            Eigen::saveMarket(matM_expl, "matM_expl.mat");

            ChSparseMatrix matK_expl;
            application.GetSystem()->GetStiffnessMatrix(&matK_expl);
            matK_expl.makeCompressed();
            //std::cout << matK_expl << std::endl;
            Eigen::saveMarket(matK_expl, "matK_expl.mat");

            ChSparseMatrix matCq_expl;
            application.GetSystem()->GetConstraintJacobianMatrix(&matCq_expl);
            matCq_expl.makeCompressed();
            //std::cout << matCq_expl << std::endl;
            Eigen::saveMarket(matCq_expl, "matCq_expl.mat");

        }
    }

    if (enable_eigenanalysis) {
        GetLog() << "\n\n===========EIGENPROBLEM======== \n";


        ChEigenAnalysis eig_analysis(application);
        eig_analysis.SetVerbose(true);
        eig_analysis.EigenAnalysis(10, 0.1, false);

        std::cout << "Time for assembly:" << eig_analysis.GetTime_Assembly() << std::endl;
        std::cout << "Time for the eigensolver:" << eig_analysis.GetTime_Eigensolve() << std::endl;

        ChVectorDynamic<> residualsNorm;
        eig_analysis.GetResidualsNorm(residualsNorm);
        std::cout << "Residual norm rev: " << residualsNorm.reverse() << std::endl;

        ChVectorDynamic<> frequencies;
        eig_analysis.GetFrequencies(frequencies);
        std::cout << "frequencies rev: " << frequencies.reverse() << std::endl;

        std::cout << "Time for assembly:" << eig_analysis.GetTime_Assembly() << std::endl;
        std::cout << "Time for the eigensolver:" << eig_analysis.GetTime_Eigensolve() << std::endl;

        Eigen::saveMarket(eig_analysis.GetEigenValues(), "eig_val.mat");
        Eigen::saveMarket(eig_analysis.GetEigenVectors(), "eig_vect.mat");

        // Strict realtime timer
        std::function<void()> simadvance_fun = std::bind(&ChEigenAnalysis::UpdateEigenMode, &eig_analysis, 1);
        // std::function<void()> render_fun = std::bind(&ChIrrApp::DrawAll, application);
        std::function<void(double)> render_fun = [&application](double dummy) { application.DrawAll(); };
        std::function<void()> pre_sim_fun = std::bind(&ChIrrApp::BeginScene, &application, true, true, irr::video::SColor(255, 0, 0, 0));

        // std::function<void(double)> render_fun = [&application, &eig_analysis](double ratio) {
        // eig_analysis.UpdateEigenMode(ratio); application.DrawAll(); };
        ChRealtimeDualStepTimer drealtime = {application, simadvance_fun, render_fun};
        drealtime.SetPreAdvanceSimulationFunction(pre_sim_fun);

        if (use_realtime){
             while (application.GetDevice()->run()) {
                //application.BeginScene(); // already embedded in realtime timer
                //application.DrawAll(); // already embedded in realtime timer
                drealtime.DoRealtimeStep();
                application.EndScene();
            }
        } else {
            while (application.GetDevice()->run()) {
                application.BeginScene();
                application.DrawAll();
                eig_analysis.UpdateEigenMode();
                application.EndScene();
            }
        }

    }
    else {
        while (application.GetDevice()->run()) {
            application.BeginScene();
            application.DrawAll();
            application.DoStep();
            application.EndScene();
        }
    }

}


