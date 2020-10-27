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

#include <iostream>
#include <algorithm>

#include "timestepper/ChState.h"
#include "timestepper/ChIntegrable.h"
#include "physics/ChSystem.h"
#include "ChModelReduction.h"

#include <unsupported/Eigen/SparseExtra>

#ifdef CHRONO_IRRLICHT
#include "chrono_irrlicht/ChIrrApp.h"
#include <irrlicht.h>
#endif

namespace chrono {

#ifdef CHRONO_IRRLICHT

class ChEigenAnalysis; // forward decl

    /// GUI class for the modal analysis
class ModalAnalysisGUIReceiver : public irr::IEventReceiver {
public:
    ModalAnalysisGUIReceiver(){}
    ~ModalAnalysisGUIReceiver(){}

    bool OnEvent(const irr::SEvent& event) override;

    void Initialize(irrlicht::ChIrrApp& application, ChEigenAnalysis& eigen_analysis);
    static const int buffer_size = 100;
    wchar_t buffer[buffer_size];

protected:
    template<typename... T>
    wchar_t* fwd_swprintf(wchar_t* ws, size_t len, const wchar_t* format, T... vararg)
    {
        swprintf(ws, len, format, vararg...);
        return ws;
    }


private:
    irrlicht::ChIrrApp* app = nullptr;
    irr::IrrlichtDevice* mdevice = nullptr;
    ChEigenAnalysis* eig_analysis = nullptr;

    irr::gui::IGUIStaticText* mode_freq_text = nullptr;
    irr::gui::IGUIStaticText* mode_selmax_text = nullptr;

    static const int modes_timestep_ID = 102;
    static const int modes_selection_ID = 103;
    static const int modes_magn_ampl_ID = 104;
    static const int modes_time_ampl_ID = 105;

};

#endif


/// Base class for static analysis

class ChEigenAnalysis {
  protected:
    ChIntegrableIIorder* integrable;
    ChSystem* m_system;

    ChTimer<> m_timer_assembly;  ///< timer for matrix assembly
    ChTimer<> m_timer_eigensolve;  ///< timer for eigen-solve

    bool verbose = false;


#ifdef CHRONO_IRRLICHT
    ModalAnalysisGUIReceiver gui_receiver;
    irrlicht::ChIrrApp* irrlicht_app = nullptr;
#endif


    ChState X0;
    ChState X;
    ChStateDelta V;
    ChStateDelta A;
    ChVectorDynamic<> L;

    double current_time = 0;
    int current_selected_mode = 1;
    double magnitude_amplification = 1;
    double time_amplification = 1;

    // Eigen analysis matrices
    ChSparseMatrix matKaug;
    ChSparseMatrix matMaug;
    ChMatrixDynamic<double> eig_val;
    ChMatrixDynamic<double> eig_vect;
    ChMatrixDynamic<double> eig_vect_col; //TODO: check if ChVectorDynamic can be used here

  public:
    /// Constructor
    ChEigenAnalysis(ChSystem& msystem) : m_system(&msystem){
        integrable = static_cast<ChIntegrableIIorder*>(m_system);
        L.resize(0);
        X0 = ChState(1, integrable);
        X = ChState(1, integrable);
        V = ChStateDelta(1, integrable);
        A = ChStateDelta(1, integrable);
    }
#ifdef CHRONO_IRRLICHT
    ChEigenAnalysis(irrlicht::ChIrrApp& application) : m_system(application.GetSystem()), irrlicht_app(&application) {
        integrable = static_cast<ChIntegrableIIorder*>(m_system);
        L.resize(0);
        X0 = ChState(1, integrable);
        X = ChState(1, integrable);
        V = ChStateDelta(1, integrable);
        A = ChStateDelta(1, integrable);
        AttachToChIrrApp(application);
    }
#endif

    /// Destructor
    virtual ~ChEigenAnalysis() {}

    /// Performs the static analysis
    virtual void EigenAnalysis(int requested_eigval = -1, double sigma = 10e-2, bool restart_sigma = false, bool Cq_preconditioning = false) {
        m_system->Setup();
        m_system->Update();

        m_timer_assembly.start();

        // loading Kaug matrix
        m_system->KRMmatricesLoad(1.0, 0, 0);
        m_system->GetSystemDescriptor()->SetMassFactor(0.0);
        m_system->GetSystemDescriptor()->ConvertToMatrixForm(&matKaug, nullptr, true);

        // scaling Cq
        if (Cq_preconditioning) {
             int constr = m_system->GetSystemDescriptor()->CountActiveConstraints();
             double diagKmean = matKaug.diagonal().mean();  // TODO: pick only the top left corner
             std::cout << "Cq scaling: " << diagKmean << std::endl;
             matKaug.bottomRows(constr) *= diagKmean;
             matKaug.rightCols(constr) *= diagKmean;
        }


        // loading Maug matrix
        m_system->KRMmatricesLoad(0, 0, 1.0);
        m_system->GetSystemDescriptor()->SetMassFactor(1.0);
        m_system->GetSystemDescriptor()->ConvertToMatrixForm(&matMaug, nullptr, false);

        assert(requested_eigval < matKaug.rows() && "ChEigenAnalyis: the requested modes exceed matrix dimension");

        // requested_eigval = min(m-1,n); m = augmented matrix dimension, n = number of active variables -> if n=m it computes n-1 eigenvalues, if n<m it computes all n eigenvalues
        if (requested_eigval == -1) 
            requested_eigval = std::min(static_cast<int>(matKaug.rows()) - 1, m_system->GetSystemDescriptor()->CountActiveVariables());

        matKaug.makeCompressed();
        matMaug.makeCompressed();

        m_timer_assembly.stop();

        m_timer_eigensolve.start();

        ChSymGEigsSolver eig_solver(matKaug, matMaug, eig_val, eig_vect);

        if (restart_sigma) {
            eig_solver.computeShift(1, sigma);
            sigma = 0.9*eig_val(0);
            std::cout << "new sigma: " << sigma << std::endl;
        }

        eig_solver.SetVerbose(verbose);
        
        eig_solver.computeShift(requested_eigval, sigma);
        //eig_solver.computeRegularInverse(requested_eigval);
        //eig_solver.computeCholesky(requested_eigval);

        m_timer_eigensolve.stop();

        eig_vect_col.resize(matKaug.rows(), 1);

        // setup main vectors
        integrable->StateSetup(X0, V, A);
        integrable->StateSetup(X, V, A);
        L.resize(integrable->GetNconstr());

        L.fill(0);// set Lagrangians to zero
        integrable->StateScatterReactions(L);  // -> system auxiliary data 

        // store the initial configuration of the system
        double T;
        integrable->StateGather(X0, V, T);

#ifdef CHRONO_IRRLICHT
        if (irrlicht_app) gui_receiver.Initialize(*irrlicht_app, *this);
#endif

    }

    void SetVerbose(bool val) { verbose = val; }

    void InjectEigenData(ChMatrixDynamic<>& eig_val_in, ChMatrixDynamic<>& eig_vect_in) {  

        eig_val = eig_val_in;
        eig_vect = eig_vect_in;

         // setup main vectors
        integrable->StateSetup(X0, V, A);
        integrable->StateSetup(X, V, A);
        L.resize(integrable->GetNconstr());

        L.fill(0);                             // set Lagrangians to zero
        integrable->StateScatterReactions(L);  // -> system auxiliary data

        // store the initial configuration of the system
        double T;
        integrable->StateGather(X0, V, T);

    #ifdef CHRONO_IRRLICHT
        if (irrlicht_app)
            gui_receiver.Initialize(*irrlicht_app, *this);
    #endif

    }

    /// Get cumulative time for assembly operations in Setup phase.
    double GetTime_Assembly() const { return m_timer_assembly(); }
    /// Get cumulative time for assembly operations in Setup phase.
    double GetTime_Eigensolve() const { return m_timer_eigensolve(); }

    /// Get eigenvalues.
    const ChMatrixDynamic<double>& GetEigenValues() const { return eig_val; }
    /// Get eigenvectors.
    const ChMatrixDynamic<double>& GetEigenVectors() const { return eig_vect; }



    void UpdateEigenMode(double step_percentage = 1)
    {
        // get the current time
        double systemdT;
        integrable->StateGather(X, V, systemdT);
        double dt = step_percentage*systemdT;
        current_time += dt;

        assert(current_selected_mode <= eig_val.rows());
            
        // evaluate the state for the selected_mode at the current timestep
        eig_vect_col = eig_vect.col(current_selected_mode-1); // TODO: chech if Eigen::Map can be used
        double pulse = sqrt(eig_val(current_selected_mode -1));
        double amplitude_corrector = magnitude_amplification*sin(time_amplification*pulse*current_time);

        eig_vect_col *= amplitude_corrector;

        // update the state
        ChStateDelta Dx;
        Dx = ChStateDelta(integrable->GetNcoords_v(), integrable);
        Dx = eig_vect_col;
        auto prova = X0.GetIntegrable();
        X = X0 + Dx;

        V.fill(0); // set V speed to zero
        integrable->StateScatter(X, V, dt, true);  // state -> system

    }

    void SetVisualizedMode(int mode_selection) { current_selected_mode = mode_selection; }
    int GetVisualizedMode() const { return current_selected_mode; }

    double GetVisualizedFrequency() const { return sqrt(eig_val(current_selected_mode - 1)) / 2 / CH_C_PI; }
    double GetFrequency(int mode) const { return sqrt(eig_val(mode - 1)) / 2 / CH_C_PI; }

    void SetTimeAmplification(double time_ampl) { time_amplification = time_ampl; }
    double GetTimeAmplification() const { return time_amplification; }
    void SetMagnitudeAmplification(double magnitude_ampl) { magnitude_amplification = magnitude_ampl; }
    double GetMagnitudeAmplification() const { return magnitude_amplification; }

    int GetComputedModes() const { return eig_val.rows(); }

    void GetResidualsNorm(ChVectorDynamic<double>& residualNorms) const {
        residualNorms.resize(eig_val.rows(), 1);
        for (int i = 0; i < eig_val.rows(); i++) {
            residualNorms(i) = ((matKaug - eig_val(i) * matMaug) * eig_vect.col(i)).norm();
        }
    }

    void GetFrequencies(ChVectorDynamic<double>& frequencies) const {
    frequencies.resize(eig_val.rows(), 1);
    for (int i = 0; i < eig_val.rows(); i++) {
        frequencies(i) = sqrt(eig_val(i)) / 2 / CH_C_PI;
        }
    }


#ifdef CHRONO_IRRLICHT
    void AttachToChIrrApp(irrlicht::ChIrrApp& application)
    {
        irrlicht_app = &application;
        irrlicht_app->SetUserEventReceiver(&gui_receiver);
    }
#endif

};

#ifdef CHRONO_IRRLICHT

    inline void ModalAnalysisGUIReceiver::Initialize(irrlicht::ChIrrApp& application, ChEigenAnalysis& eigen_analysis) {
        // store pointer to physical system & other stuff so we can tweak them by user keyboard
        app = &application;
        eig_analysis = &eigen_analysis;
        mdevice = application.GetDevice();

        auto win_size = mdevice->getVideoDriver()->getScreenSize();

        auto mode_tab = mdevice->getGUIEnvironment()->addTabControl(irr::core::rect<irr::s32>(win_size.Width - 250, 10, win_size.Width - 20, win_size.Height / 2), 0, true, true);
        mode_tab->addTab(L"Modes");

        auto left_margin = 10;
        auto text_width = 120;
        auto box_width = 80;
        auto item_height = 15;
        auto cum_y_pos = 35;

        mdevice->getGUIEnvironment()->addStaticText(L"Modes timestep", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_timestep_box = mdevice->getGUIEnvironment()->addEditBox(L"",
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width, cum_y_pos + item_height), true, mode_tab, modes_timestep_ID); cum_y_pos += item_height;

        mdevice->getGUIEnvironment()->addStaticText(L"Mode selected", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_sel_box = mdevice->getGUIEnvironment()->addEditBox(L"", 
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width*0.6, cum_y_pos + item_height), true, mode_tab, modes_selection_ID);
        mode_selmax_text = mdevice->getGUIEnvironment()->addStaticText(L"/",
                                                    irr::core::rect<irr::s32>(left_margin + text_width + box_width*0.6, cum_y_pos, left_margin + text_width + box_width, cum_y_pos + item_height), false, true, mode_tab);
        cum_y_pos += item_height;

        mdevice->getGUIEnvironment()->addStaticText(L"Magnitude amplification", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_magn_ampl_box = mdevice->getGUIEnvironment()->addEditBox(L"", 
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width, cum_y_pos + item_height), true, mode_tab, modes_magn_ampl_ID); cum_y_pos += item_height;

        mdevice->getGUIEnvironment()->addStaticText(L"Time amplification", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_time_ampl_box = mdevice->getGUIEnvironment()->addEditBox(L"",
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width, cum_y_pos + item_height), true, mode_tab, modes_time_ampl_ID); cum_y_pos += item_height;

        mode_freq_text = mdevice->getGUIEnvironment()->addStaticText(L"Frequency",
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);

        mode_timestep_box->setText(fwd_swprintf(buffer, buffer_size, L"%.2e", app->GetSystem()->GetStep()));
        mode_sel_box->setText(fwd_swprintf(buffer, buffer_size, L"%d", eig_analysis->GetVisualizedMode()));
        mode_magn_ampl_box->setText(fwd_swprintf(buffer, buffer_size, L"%.1g", eig_analysis->GetMagnitudeAmplification()));
        mode_time_ampl_box->setText(fwd_swprintf(buffer, buffer_size, L"%.1g", eig_analysis->GetTimeAmplification()));
        mode_freq_text->setText(fwd_swprintf(buffer, buffer_size, L"Frequency: %.1fHz", eig_analysis->GetVisualizedFrequency()));
        mode_selmax_text->setText(fwd_swprintf(buffer, buffer_size, L"/%d", eig_analysis->GetComputedModes()));

    }

    inline bool ModalAnalysisGUIReceiver::OnEvent(const irr::SEvent& event) {

        if (event.EventType == irr::EET_GUI_EVENT) {
            irr::s32 id = event.GUIEvent.Caller->getID();
            irr::gui::IGUIEnvironment* env = mdevice->getGUIEnvironment();

            switch (event.GUIEvent.EventType) {
                case irr::gui::EGET_EDITBOX_ENTER:
                    switch (id)
                    {
                        case modes_timestep_ID:
                        {
                            double dt = atof(irr::core::stringc(static_cast<irr::gui::IGUIEditBox*>(event.GUIEvent.Caller)->getText()).c_str());
                            if (dt > 0)
                                app->GetSystem()->SetStep(dt);
                            break;
                        }

                        case modes_selection_ID:
                        {
                            int mode = atoi(irr::core::stringc(static_cast<irr::gui::IGUIEditBox*>(event.GUIEvent.Caller)->getText()).c_str());
                            eig_analysis->SetVisualizedMode(std::min(eig_analysis->GetComputedModes(), mode));
                            mode_freq_text->setText(fwd_swprintf(buffer, buffer_size, L"Frequency: %.1fHz", eig_analysis->GetVisualizedFrequency()));
                            mode_selmax_text->setText(fwd_swprintf(buffer, buffer_size, L"/%d", eig_analysis->GetComputedModes()));
                            break;
                        }

                        case modes_magn_ampl_ID:
                        {
                            double magn_ampl = atof(irr::core::stringc(static_cast<irr::gui::IGUIEditBox*>(event.GUIEvent.Caller)->getText()).c_str());
                            eig_analysis->SetMagnitudeAmplification(magn_ampl);
                            std::cout << "new amplitude is " << eig_analysis->GetMagnitudeAmplification() << std::endl;
                            break;
                        }
                            

                        case modes_time_ampl_ID:
                        {
                            double time_ampl = atof(irr::core::stringc(static_cast<irr::gui::IGUIEditBox*>(event.GUIEvent.Caller)->getText()).c_str());
                            eig_analysis->SetTimeAmplification(time_ampl);
                            break;
                        }
                    }
                    break;


            }
        }

        return false;
    }

#endif

class ChRealtimeDualStepTimer
{
public:
    
    ChRealtimeDualStepTimer(ChSystem& chrono_system, std::function<void()>& simulation_advance_fun, std::function<void(double)>& rendering_fun) : 
        sim_adv_fun(simulation_advance_fun),
        render_fun(rendering_fun),
        m_system(chrono_system) {}
    
    ChRealtimeDualStepTimer(ChSystem& chrono_system, std::function<void()>& simulation_advance_fun, std::function<void()>& rendering_fun) : 
        sim_adv_fun(simulation_advance_fun),
        m_system(chrono_system)
    {
        SetRenderFunction(rendering_fun);
    }

#ifdef CHRONO_IRRLICHT
    ChRealtimeDualStepTimer(irrlicht::ChIrrApp& irrlicht_app, std::function<void()>& simulation_advance_fun, std::function<void(double)>& rendering_fun) :
        sim_adv_fun(simulation_advance_fun),
        render_fun(rendering_fun),
        m_system(*irrlicht_app.GetSystem()),
        irr_app(&irrlicht_app) {}
    
    ChRealtimeDualStepTimer(irrlicht::ChIrrApp& irrlicht_app, std::function<void()>& simulation_advance_fun, std::function<void()>& rendering_fun) :
        sim_adv_fun(simulation_advance_fun),
        m_system(*irrlicht_app.GetSystem()),
        irr_app(&irrlicht_app)
    {
        SetRenderFunction(rendering_fun);
    }
#endif
    
    void DoRealtimeStep()
    {
#ifdef CHRONO_IRRLICHT
        if (irr_app && irr_app->GetPaused())
            return;
#endif

        if (on_realtime && try_realtime) delay_with_wallclock += wallclock_timer.GetTimeSecondsIntermediate();
        wallclock_timer.start();

        // eventually catch user input here
        if (pre_adv_fun && state_updates)
            pre_adv_fun();

        simulation_timestep = m_system.GetStep();
        state_updates = 0;
        //std::cout << delay_with_wallclock << " initial delay" << std::endl;
        while (delay_with_wallclock >= simulation_timestep && try_realtime)
        {
            double start = wallclock_timer.GetTimeSecondsIntermediate();
            sim_adv_fun();
            ++state_updates;
            delay_with_wallclock -= simulation_timestep;
            if (wallclock_timer.GetTimeSecondsIntermediate() - start > simulation_timestep) break;
        }
        on_realtime = delay_with_wallclock < simulation_timestep;
        std::cout << delay_with_wallclock << " < " << simulation_timestep << ": " << on_realtime << std::endl;
        std::cout << "State updates: " << state_updates << std::endl;

        //std::cout << "It is " << (on_realtime ? "" : "not") << " realtime." << std::endl;

        if (state_updates) render_fun(delay_with_wallclock / simulation_timestep);

    }

    void SetAdvanceSimulationFunction(std::function<void()>& simulation_advance_fun) { sim_adv_fun = simulation_advance_fun; }
    void SetPreAdvanceSimulationFunction(std::function<void()>& pre_simulation_advance_fun) { pre_adv_fun = pre_simulation_advance_fun; }
    void SetRenderFunction(std::function<void(double)>& rendering_fun) { render_fun = rendering_fun; }
    void SetRenderFunction(std::function<void()>& rendering_fun) { render_fun = [&rendering_fun](double dummy) {return rendering_fun(); }; }

    void SetRealtime(bool on_off) { try_realtime = on_off; }

    //void SetMaxFrameRate(double framerate) { maxframerate = framerate; }

    bool IsRealtime() const { return on_realtime; }

#ifdef CHRONO_IRRLICHT
    void AttachToGUI(irrlicht::ChIrrApp& irrlicht_app) { irr_app = &irrlicht_app; }
#endif

private:
    ChTimer<double> wallclock_timer;
    std::function<void()> sim_adv_fun;
    std::function<void(double)> render_fun;
    std::function<void()> pre_adv_fun;
    ChSystem& m_system;
#ifdef CHRONO_IRRLICHT
    irrlicht::ChIrrApp* irr_app = nullptr;
#endif
    double previous_wallclock_time = -1;
    double simulation_timestep = 1e-2;
    bool loop_started = false;
    bool on_realtime = false;
    bool try_realtime = true;
    int state_updates = 0;
    //double maxframerate = 60;

    double delay_with_wallclock = 0;

};

}  // END_OF_NAMESPACE____
#endif
