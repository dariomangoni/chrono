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

#ifndef CHEIGENANALYSIS_H
#define CHEIGENANALYSIS_H

#include "core/ChVectorDynamic.h"
#include "timestepper/ChState.h"
#include "timestepper/ChIntegrable.h"
#include "physics/ChSystem.h"
#include "ChSolverMR.h"

namespace chrono {

/// Base class for static analysis

class ChEigenAnalysis {
  protected:
    ChIntegrableIIorder* integrable;
    ChSystem* m_system;

    ChState X0;
    ChState X;
    ChStateDelta V;
    ChStateDelta A;
    ChVectorDynamic<> L;

    //double magnitude_amplification = 1;
    //double time_magnification = 1;
    double current_time = 0;

    ChMatrixDynamic<double> eig_val;
    ChMatrixDynamic<double> eig_vect;
    ChMatrixDynamic<double> eig_vect_col;

  public:
    /// Constructor
    ChEigenAnalysis(ChSystem& msystem) : m_system(&msystem){
        integrable = static_cast<ChIntegrableIIorder*>(m_system);
        L.Reset(0);
        X0.Reset(1, integrable);
        X.Reset(1, integrable);
        V.Reset(1, integrable);
        A.Reset(1, integrable);
    };

    /// Destructor
    virtual ~ChEigenAnalysis() {}

    /// Performs the static analysis
    virtual void EigenAnalysis(int total_modes)
    {
        m_system->Setup();
        m_system->Update();
        ChCSR3Matrix matK, matM;
        m_system->GetStiffnessMatrix(&matK);
        m_system->GetMassMatrix(&matM);

        assert(total_modes < matK.GetNumRows() && "ChEigenAnalyis: the requested modes exceed matrix dimension");

        matK.Compress();
        matM.Compress();

        //matK.VerifyMatrix();
        //matM.VerifyMatrix();

        //matK.ExportToDatFile("C:/K", 6);
        //matM.ExportToDatFile("C:/M", 6);

        ChSymGEigsSolver eig_solver(matK, matM, eig_val, eig_vect);
        eig_solver.compute(total_modes);
        eig_vect_col.Reset(matK.GetNumRows(), 1);

        // setup main vectors
        integrable->StateSetup(X0, V, A);
        integrable->StateSetup(X, V, A);
        L.Reset(integrable->GetNconstr());

        L.FillElem(0);// set Lagrangians to zero
        integrable->StateScatterReactions(L);  // -> system auxiliary data 

        // store the initial configuration of the system
        double T;
        integrable->StateGather(X0, V, T);

    }

    void UpdateEigenMode(int selected_mode, double magnitude_amplification = 1, double time_amplification = 1)
    {
        // get the current time
        double T;
        integrable->StateGather(X, V, T);
        current_time += T;

        // evaluate the state for the selected_mode a thte current timestep
        eig_vect_col.PasteClippedMatrix(&eig_vect, 0, selected_mode-1, eig_vect.GetRows(), 1, 0, 0);
        double pulse = sqrt(eig_val(selected_mode-1));
        double amplitude_corrector = magnitude_amplification*sin(time_amplification*pulse*T);
        eig_vect_col.MatrScale(amplitude_corrector);

        // update the state
        ChStateDelta Dx;
        Dx.Reset(integrable->GetNcoords_v(), GetIntegrable());
        Dx.CopyFromMatrix(eig_vect_col);
        auto prova = X0.GetIntegrable();
        X = X0 + Dx;

        V.FillElem(0); // set V speed to zero
        integrable->StateScatter(X, V, T);  // state -> system

    }

    /// Access the lagrangian multipliers, if any
    virtual ChVectorDynamic<>& get_L() { return L; }

    /// Get the integrable object
    ChIntegrable* GetIntegrable() const { return integrable; }

    /// Access the state, position part, at current analysis
    virtual ChState& get_X() { return X; }

    /// Access the state, speed part, at current analysis
    virtual ChStateDelta& get_V() { return V; }

    /// Access the acceleration, at current analysis
    virtual ChStateDelta& get_A() { return A; }

};

}  // END_OF_NAMESPACE____
#endif
