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
//
// Irrlicht-based visualization wrapper for modal. This class is a derived
// from ChVisualSystemIrrlicht and provides functionality to draw mode shapes
// =============================================================================

#ifndef CH_MODAL_VISUAL_SYSTEM_IRRLICHT_H
#define CH_MODAL_VISUAL_SYSTEM_IRRLICHT_H

#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtils.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono_modal/ChModalSolver.h"

namespace chrono {
namespace modal {

/// @addtogroup modal_vis
/// @{
///

template <typename ScalarType>
class ChModalVisualSystemIrrlicht : public irrlicht::ChVisualSystemIrrlicht {
  public:
    ChModalVisualSystemIrrlicht() {
        SetWindowSize(1024, 768);
        SetWindowTitle("Chrono::Modal");
    }

    virtual ~ChModalVisualSystemIrrlicht();

    /// Initialize the visualization system.
    virtual void Initialize() override;

    virtual void BeginScene(bool backBuffer, bool zBuffer, ChColor color) override;

    using irrlicht::ChVisualSystemIrrlicht::BeginScene;

    /// Render the Irrlicht scene and additional visual elements.
    virtual void Render() override;

    void SetMode(int mode) { m_selected_mode = ChClamp(mode, 0, m_eigvects->cols()); }

    void ResetTimer() { m_timer.reset(); }

    void SetAmplitude(double amplitude) { m_amplitude = amplitude; }

    void SetTimeFactor(double time_factor) { m_time_factor = time_factor; }

    void ResetInitialState();

    void AttachAssembly(const ChAssembly& assembly,
                        const ChMatrixDynamic<ScalarType>& eigvects,
                        const ChVectorDynamic<double>& freq) {
        m_assembly = const_cast<ChAssembly*>(&assembly);

        UpdateModes(eigvects, freq);
    }

    void UpdateModes(const ChMatrixDynamic<ScalarType>& eigvects, const ChVectorDynamic<double>& freq) {
        m_eigvects = &eigvects;
        m_freq = &freq;
        m_selected_mode = 0;
    }

  protected:
    ChAssembly* m_assembly = nullptr;
    const ChMatrixDynamic<ScalarType>* m_eigvects = nullptr;
    const ChVectorDynamic<double>* m_freq = nullptr;

    mutable int m_selected_mode = 0;
    ChTimer m_timer;

    ChState m_assembly_initial_state;
    double m_amplitude = 1.0;
    double m_time_factor = 1.0;

    template <typename U = ScalarType>
    typename std::enable_if<std::is_same<U, double>::value, ChVectorDynamic<double>>::type GetModeShape(
        const ChVectorDynamic<double>& eigv,
        double angle) {
        return eigv.stableNormalized() * sin(angle);
    }

    template <typename U = ScalarType>
    typename std::enable_if<std::is_same<U, std::complex<double>>::value, ChVectorDynamic<double>>::type GetModeShape(
        const ChVectorDynamic<std::complex<double>>& eigv,
        double angle) {
        return eigv.stableNormalized().real() * sin(angle) + eigv.stableNormalized().imag() * cos(angle);
    }
};

template <typename ScalarType>
ChModalVisualSystemIrrlicht<ScalarType>::~ChModalVisualSystemIrrlicht() {}

template <typename ScalarType>
void ChModalVisualSystemIrrlicht<ScalarType>::Initialize() {
    ChVisualSystemIrrlicht::Initialize();

    m_timer.reset();
    ResetInitialState();
}

template <typename ScalarType>
inline void ChModalVisualSystemIrrlicht<ScalarType>::BeginScene(bool backBuffer, bool zBuffer, ChColor color) {
    ChVisualSystemIrrlicht::BeginScene(backBuffer, zBuffer, color);

    // check if the timer is running, if not, start it

    if (m_timer.GetTimeMilliseconds() == 0.0)
        m_timer.start();

    double angle = m_time_factor * CH_2PI * m_freq->coeff(m_selected_mode) * m_timer.GetTimeSeconds();

    ChState assembly_state_new;
    ChStateDelta assembly_state_delta;
    ChStateDelta assembly_v_dummy;

    double time_dummy = 0;

    assembly_state_new.setZero(m_assembly->GetNumCoordsPosLevel(), nullptr);
    assembly_state_delta.setZero(m_assembly->GetNumCoordsVelLevel(), nullptr);
    assembly_v_dummy.setZero(m_assembly->GetNumCoordsVelLevel(), nullptr);

    assembly_state_delta = m_amplitude * GetModeShape<>(m_eigvects->col(m_selected_mode), angle);

    m_assembly->IntStateIncrement(0, assembly_state_new, m_assembly_initial_state, 0, assembly_state_delta);
    m_assembly->IntStateScatter(0, assembly_state_new, 0, assembly_v_dummy, time_dummy, true);

    m_assembly->Update();

    OnUpdate(m_systems[0]);
}

template <typename ScalarType>
void ChModalVisualSystemIrrlicht<ScalarType>::Render() {
    ChVisualSystemIrrlicht::Render();
}

template <typename ScalarType>
void ChModalVisualSystemIrrlicht<ScalarType>::ResetInitialState() {
    if (!m_assembly)
        return;

    m_assembly->Setup();
    m_assembly_initial_state.setZero(m_assembly->GetNumCoordsPosLevel(), nullptr);

    ChStateDelta dummy_state_delta;
    dummy_state_delta.setZero(m_assembly->GetNumCoordsVelLevel(), nullptr);

    double time;
    m_assembly->IntStateGather(0, m_assembly_initial_state, 0, dummy_state_delta, time);
}

// @} modal_vis

}  // end namespace modal
}  // end namespace chrono

#endif
