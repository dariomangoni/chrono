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
#include <ChModalAnalysisArpack.h>

//using ChSparseMatrixRef = Eigen::Ref<Eigen::SparseMatrix<double, Eigen::RowMajor, int>>;
//using ChSparseMatrixConstRef = const Eigen::Ref<const Eigen::SparseMatrix<double, Eigen::RowMajor, int>>;


//#include "chrono_matlab/ChMatlabEngine.h"
//#include "chrono_matlab/ChSolverMatlab.h"

// Remember to use the namespace 'chrono' because all classes 
// of Chrono::Engine belong to this namespace and its children...



using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

std::string outPath = "C:/workspace/chrono_build/bin/";

bool enable_eigenanalysis = true;
bool fix_by_explicit_constraint = true;
bool block_out_of_plane = true;
bool use_realtime = false;
bool save_matrix = true;
bool use_eigen_eigensolver = true;

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

    double beam_wy = 0.06;                  //[m]
    double beam_wz = 0.04;                  //[m]
    double beam_L = 2;                      //[m]
    double E = 210e9;                      //[Pa]
    double nu = 0.3;
    double G = E/2/(1+nu);                     //[Pa]
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

    msection->SetBeamRaleyghDamping(0.013);

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


    int num_elements = 16;

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

    if(fix_by_explicit_constraint){
        //Fix link between beam and deadboy
        auto my_FixLink = chrono_types::make_shared<ChLinkMateFix>();
        my_FixLink->SetName("FixLink");
        my_FixLink->Initialize(my_deadbody, builder.GetLastBeamNodes()[0]);
        my_system.AddLink(my_FixLink);
        }
    else{
        builder.GetLastBeamNodes()[0]->SetFixed(true);
    }

    builder.GetLastBeamNodes().back()->SetForce(ChVector<>(0, -1.0, 0));
    
    if (block_out_of_plane){
        // Modes only on the vertical plane
        for (int i = 1; i < builder.GetLastBeamNodes().size()-1; i++) {
            auto abearing = chrono_types::make_shared<ChLinkMateGeneric>(true, false, true, true, true, false);
            abearing->Initialize(builder.GetLastBeamNodes()[i], my_deadbody, ChFrame<>(builder.GetLastBeamNodes()[i]->GetPos()));
            my_system.Add(abearing);
        }
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

        //ChSparseMatrix matKaug;
        //application.GetSystem()->KRMmatricesLoad(1.0, 0, 0);
        //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
        //application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matKaug, nullptr, true);
        ////std::cout << matKaug << std::endl;
        //Eigen::saveMarket(matKaug, outPath+"matKaug.mat");

        //ChSparseMatrix matMaug;
        //application.GetSystem()->KRMmatricesLoad(0, 0, 1.0);
        //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
        //application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matMaug, nullptr, false);
        ////std::cout << matMaug << std::endl;
        //Eigen::saveMarket(matMaug, outPath+"matMaug.mat");

        ChSparseMatrix matM;
        application.GetSystem()->GetMassMatrix(&matM);
        matM.makeCompressed();
        //std::cout << matM << std::endl;
        Eigen::saveMarket(matM, outPath+"matM.mat");

        ChSparseMatrix matR;
        application.GetSystem()->GetDampingMatrix(&matR);
        matR.makeCompressed();
        //std::cout << matR << std::endl;
        Eigen::saveMarket(matR, outPath+"matR.mat");

        ChSparseMatrix matK;
        application.GetSystem()->GetStiffnessMatrix(&matK);
        matK.makeCompressed();
        //std::cout << matK << std::endl;
        Eigen::saveMarket(matK, outPath+"matK.mat");

        ChSparseMatrix matCq;
        application.GetSystem()->GetConstraintJacobianMatrix(&matCq);
        matCq.makeCompressed();
        //std::cout << matCq << std::endl;
        Eigen::saveMarket(matCq, outPath+"matCq.mat");

    }

    if (enable_eigenanalysis) {
        GetLog() << "\n\n===========EIGENPROBLEM======== \n";
        //{
        //    ChSparseMatrix matA;
        //    ChSparseMatrix matB;

        //    double preshift = 0.1;

        //    auto m_system = application.GetSystem();

        //    int n = m_system->GetNcoords_w();
        //    int m = m_system->GetNconstr();

        //    matA.resize(2*(n+m), 2*(n+m));
        //    matB.resize(2*(n+m), 2*(n+m));

        //    double preconditioner = 1;

        //    // A = [K_,    0; 0, I]
        //    // B = [-R_, -M_; I, 0]
        //    // where:
        //    //       K_ = [K, Cq'; Cq, 0]
        //    //       R_ = [R, 0;    0, 0]
        //    //       M_ = [M, 0;    0, 0]

        //    // loading Kaug matrix
        //    m_system->KRMmatricesLoad(preconditioner, preshift, preshift*preshift);
        //    m_system->GetSystemDescriptor()->SetMassFactor(0.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matA, 0, 0, false);

        //    // loading constraint Jacobian
        //    m_system->ConstraintsLoadJacobians();
        //    m_system->GetSystemDescriptor()->BuildCq(&matA, n, 0, false, false);
        //    m_system->GetSystemDescriptor()->BuildCq(&matA, 0, n, false, true);

        //    // loading Raug matrix
        //    m_system->KRMmatricesLoad(0, 1.0, 2*preshift);
        //    m_system->GetSystemDescriptor()->SetMassFactor(0.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matB, 0, 0, false);

        //    // loading Maug matrix
        //    m_system->KRMmatricesLoad(0, 0, 1.0);
        //    m_system->GetSystemDescriptor()->SetMassFactor(1.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matB, 0, m+n, false);

        //    matB *= -1;
        //
        //    for (auto i = 0; i<n+m; ++i){
        //        matB.coeffRef(i+n+m, i) = 1.0;
        //        matA.coeffRef(i+n+m, i+n+m) = 1.0;
        //    }


        //    matA.makeCompressed();
        //    matB.makeCompressed();

        //    Eigen::saveMarket(matA, "matA.mat");
        //    Eigen::saveMarket(matB, "matB.mat");

        //}
        ////////////////////////////
        //{
        //    ChSparseMatrix matA;
        //    ChSparseMatrix matB;

        //    double preshift = 1;

        //    auto m_system = application.GetSystem();

        //    int n = m_system->GetNcoords_w();
        //    int m = m_system->GetNconstr();

        //    matA.resize(2*(n+m), 2*(n+m));
        //    matB.resize(2*(n+m), 2*(n+m));

        //    double preconditioner = 1;


        //    //m_system->KRMmatricesLoad(1, 0, 0);
        //    //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
        //    //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
        //    //std::cout << "K: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

        //    //m_system->KRMmatricesLoad(0, 1, 0);
        //    //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
        //    //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
        //    //std::cout << "R: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

        //    //m_system->KRMmatricesLoad(0, 0, 1);
        //    //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
        //    //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
        //    //std::cout << "M: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

        //    // loading Kaug matrix
        //    m_system->KRMmatricesLoad(preconditioner, preshift, preshift*preshift);
        //    m_system->GetSystemDescriptor()->SetMassFactor(0.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, 0, false);
        //    

        //    // loading constraint Jacobian
        //    m_system->ConstraintsLoadJacobians();
        //    m_system->GetSystemDescriptor()->BuildCq(&matA, m+n+n, 0, false, false);
        //    m_system->GetSystemDescriptor()->BuildCq(&matA, m+n  , n, false, true);

        //    // loading Raug matrix
        //    m_system->KRMmatricesLoad(0, 1.0, 2*preshift);
        //    m_system->GetSystemDescriptor()->SetMassFactor(0.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);

        //    matA*=-1;

        //    // loading Maug matrix
        //    m_system->KRMmatricesLoad(0, 0, 1.0);
        //    m_system->GetSystemDescriptor()->SetMassFactor(1.0);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matA,   0, m+n, false);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matB,   0,   0, false);
        //    m_system->GetSystemDescriptor()->BuildKRM(&matB, m+n, m+n, false);

        //    matA.makeCompressed();
        //    matB.makeCompressed();

        //    std::cout << "Writing matrices: " << std::endl;
        //    //Eigen::saveMarket(matA, outPath+"matA_dndrv2.mat");
        //    Eigen::saveMarket(matB, outPath+"matB_dndrv2.mat");

        //}

        {
            ChSparseMatrix matA;
            ChSparseMatrix matB;

            double preshift = 1;

            auto m_system = application.GetSystem();

            int n = m_system->GetNcoords_w();
            int m = m_system->GetNconstr();

            matA.resize(2*(n+m), 2*(n+m));
            matB.resize(2*(n+m), 2*(n+m));

            double preconditioner = 1;


            //m_system->KRMmatricesLoad(1, 0, 0);
            //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
            //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
            //std::cout << "K: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

            //m_system->KRMmatricesLoad(0, 1, 0);
            //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
            //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
            //std::cout << "R: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

            //m_system->KRMmatricesLoad(0, 0, 1);
            //application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
            //m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
            //std::cout << "M: matA.coeff(m+n, m+n): " << matA.coeff(m+n, m+n) << std::endl;

            // loading Kaug matrix
            m_system->KRMmatricesLoad(preconditioner, preshift, preshift*preshift);
            
            m_system->GetSystemDescriptor()->BuildKRM(&matA, 0, 0, false);
            

            // loading constraint Jacobian
            m_system->ConstraintsLoadJacobians();
            m_system->GetSystemDescriptor()->BuildCq(&matA, n, 0, false, false);
            m_system->GetSystemDescriptor()->BuildCq(&matA, 0, n, false, true);

            // loading Raug matrix
            m_system->KRMmatricesLoad(0, 1.0, 2*preshift);
            m_system->GetSystemDescriptor()->SetMassFactor(0.0);
            m_system->GetSystemDescriptor()->BuildKRM(&matA, 0, m+n, false);

            matA*=-1;

            // loading Maug matrix
            m_system->KRMmatricesLoad(0, 0, 1.0);
            m_system->GetSystemDescriptor()->SetMassFactor(1.0);
            m_system->GetSystemDescriptor()->BuildKRM(&matA, m+n, m+n, false);
            m_system->GetSystemDescriptor()->BuildKRM(&matB, 0, m+n, false);
            m_system->GetSystemDescriptor()->BuildKRM(&matB, m+n, 0, false);

            matA.makeCompressed();
            matB.makeCompressed();

            std::cout << "Writing matrices: " << std::endl;
            Eigen::saveMarket(matA, outPath+"matA_dndrv2.mat");
            Eigen::saveMarket(matB, outPath+"matB_dndrv2.mat");

        }


        ////////////////////////////


        //void ModeSolve::DoEigenSolve(	std::vector<Eigen::MatrixXd>& system_matrix, int mode_count, std::vector<OutputEigenResults>& res_eigen_values_and_vectors)
    
    double freq_shift = 100;

    ChSparseMatrix matM;
    application.GetSystem()->GetMassMatrix(&matM);

    ChSparseMatrix matR;
    application.GetSystem()->GetDampingMatrix(&matR);

    ChSparseMatrix matK;
    application.GetSystem()->GetStiffnessMatrix(&matK);

    ChSparseMatrix matCq;
    application.GetSystem()->GetConstraintJacobianMatrix(&matCq);
	
	Eigen::MatrixXd Mall = matM;
	Eigen::MatrixXd Kall = matK;
	Eigen::MatrixXd Rall = matR;
	Eigen::MatrixXd Cqall = matCq;

	// Cq matrix null space（kernel）
	Eigen::MatrixXd Cq_null_space = Cqall.fullPivLu().kernel();
	Eigen::MatrixXd M_hat = Cq_null_space.transpose() * Mall * Cq_null_space;
	Eigen::MatrixXd K_hat = Cq_null_space.transpose() * Kall * Cq_null_space;
	Eigen::MatrixXd R_hat = Cq_null_space.transpose() * Rall * Cq_null_space;
	
	// frequency-shift，matrix singularity issue does not need to be faced here, can be set to zero
	//REAL freq_shift = 0.0;  // any value from 1 to 10 (?)
	Eigen::MatrixXd M_bar = M_hat;
	Eigen::MatrixXd K_bar = pow(freq_shift, 2) * M_hat + freq_shift * R_hat + K_hat;
	Eigen::MatrixXd R_bar = 2*freq_shift*M_hat + R_hat;
	Eigen::MatrixXd M_bar_inv = M_bar.inverse();  //it is not too slow to directly use the dense matrix to find the inverse
	
	// Generate the A matrix of the state equation, whose eigenvalues ​​are the modal frequencies
	int dim = M_bar.rows();
	Eigen::MatrixXd A_tilde(2 * dim, 2 * dim);  // Combine the system matrix A matrix, dense matrix.
	A_tilde << Eigen::MatrixXd::Zero(dim, dim), Eigen::MatrixXd::Identity(dim, dim), -M_bar_inv * K_bar, -M_bar_inv * R_bar;
	
	// Call EIGEN3, dense matrix to directly solve the eigenvalues ​​and eigenvectors
	Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver(A_tilde);
	//Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues() + freq_shift;
    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues(); // TODO: add freq_shift
	Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();

    //std::vector<double> undamped_pulsation;
    //std::vector<double> damped_pulsation;
    //std::vector<double> damping_ratio;

    //undamped_pulsation.resize(eigen_values.size());
    //damped_pulsation.resize(eigen_values.size());
    //damping_ratio.resize(eigen_values.size());

    //for (int i = 0; i < eigen_values.size(); i++) {
    //    undamped_pulsation[i] = std::abs(eigen_values(i));
    //    damped_pulsation[i] = eigen_values(i).imag();
    //    damping_ratio[i] = -eigen_values(i).real() / std::abs(eigen_values(i));
    //}


	//std::vector<EigenResults> all_eigen_values_and_vectors;
	//for (int i = 0; i < eigen_values.size(); i++) {
	//	EigenResults vec_and_val(
	//		std::abs(eigen_values(i)) / (2 * PI),                   // index：0     Absolute value of eigenvalue, undamped modal frequency
	//		eigen_values(i).imag() / (2 * PI),                      // index：1     The imaginary part of the eigenvalue, with damped modal frequencies
	//		-eigen_values(i).real() / std::abs(eigen_values(i)),    // index：2     The real part of the eigenvalue represents the damping ratio
	//		eigen_values(i),                                        // index：3     Plural eigenvalues
	//		eigen_vectors.col(i)                                    // index：4     Plural mode shapes
	//	);
	//	if (std::abs(std::get<1>(vec_and_val)) > 0.001 &&  // Filter out the eigenvalues whose imaginary part is greater than zero
	//		std::abs(std::get<2>(vec_and_val)) < 0.95 ) {  // Filter out eigenvalues with damping ratio less than 0.95
	//		all_eigen_values_and_vectors.push_back(vec_and_val);
	//	}
	//}
    
	//// Sorting
	// TODO

	// Organize the results and return the calculation results
	int middle_number = static_cast<int>(eigen_values.size()/2);  // The eigenvalues filtered out are conjugate complex roots, just take half
	int DOF_counts = static_cast<int>(eigen_vectors.rows()/2);  // The number of degrees of freedom in the model, used when extracting the mode shape. NOTE: The number of degrees of freedom is not equal to that of the blade beam model

	//int jj = 0; int i = 0; int index_mode = 0;
 //   Eigen::MatrixXd vel_shapes;
 //   Eigen::MatrixXd pos_shapes;
 //   int mode_count = 5;
	//while (i < mode_count) {  
	//	index_mode = jj + middle_number;
	//	if (damping_ratio[i] < 0.95) { // Only extract the modes with a modal damping ratio less than 0.95
	//		//std::get<0>(output_temp) = std::get<0>(all_eigen_values_and_vectors.at(index_mode)); // Undamped modal frequency
	//		//std::get<1>(output_temp) = std::get<1>(all_eigen_values_and_vectors.at(index_mode)); // Damped modal frequency
	//		//std::get<2>(output_temp) = std::get<2>(all_eigen_values_and_vectors.at(index_mode)); // Damping ratio
	//		//std::get<3>(output_temp) = std::get<3>(all_eigen_values_and_vectors.at(index_mode)); // Plural eigenvalues
	//		//vector_temp = std::get<4>(all_eigen_values_and_vectors.at(index_mode));
	//		pos_shapes = Cq_null_space * eigen_vectors.col(i).head(DOF_counts); // Mode shape, displacement term
	//		vel_shapes = Cq_null_space * eigen_vectors.col(i).tail(DOF_counts); // Mode shape, velocity term
	//		i++;
	//	} 
	//	jj++;
	//}

    std::cout << eigen_values << std::endl;
    Eigen::saveMarket(eigen_values, outPath+"Eigen_eigval.mat");
    Eigen::saveMarket(eigen_vectors, outPath+"Eigen_eigvect.mat");


        //////////


        ChEigenAnalysis eig_analysis(application);
        eig_analysis.SetVerbose(true);

        //if (use_eigen_eigensolver)
        //    {eig_analysis.InjectEigenData(ChMatrixDynamic<>(undamped_pulsation), ChMatrixDynamic<>(eigen_vectors));}
        //else {
                    eig_analysis.EigenAnalysis(16, 100, false);

        //}

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

        Eigen::saveMarket(eig_analysis.GetEigenValues(), outPath+"eig_val.mat");
        Eigen::saveMarket(eig_analysis.GetEigenVectors(), outPath+"eig_vect.mat");

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


