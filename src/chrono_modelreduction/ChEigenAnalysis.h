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
#include "ChModelReduction.h"

#ifdef CHRONO_IRRLICHT
#include <irrlicht.h>
#endif


namespace chrono {

#ifdef CHRONO_IRRLICHT
// The ModalAnalysisGUIReceiver class will be used to manage input
// from the GUI graphical user interface (the interface will
// be created with the basic -yet flexible- platform
// independent toolset of Irrlicht).

class ChEigenAnalysis;

class ModalAnalysisGUIReceiver : public irr::IEventReceiver {
public:
    ModalAnalysisGUIReceiver(){}
    ~ModalAnalysisGUIReceiver(){}

    bool OnEvent(const irr::SEvent& event) override;

    void Initialize(irrlicht::ChIrrApp& application, ChEigenAnalysis& eigen_analysis);

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

#ifdef CHRONO_IRRLICHT
    ModalAnalysisGUIReceiver gui_receiver;
#endif


    ChState X0;
    ChState X;
    ChStateDelta V;
    ChStateDelta A;
    ChVectorDynamic<> L;

    double current_time = 0;
    int current_selected_mode = 1;
    ChRealtimeStepTimer real_timer;
    bool try_realtime = true;
    double magnitude_amplification = 1;
    double time_amplification = 1;

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
    virtual void EigenAnalysis(int total_modes = -1)
    {
        m_system->Setup();
        m_system->Update();
        ChCSR3Matrix matK, matM;
        m_system->GetStiffnessMatrix(&matK);
        m_system->GetMassMatrix(&matM);

        assert(total_modes < matK.GetNumRows() && "ChEigenAnalyis: the requested modes exceed matrix dimension");
        if (total_modes == -1)
            total_modes = matK.GetNumRows() - 1;

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

    void UpdateEigenMode()
    {
        // get the current time
        //double systemdT;
        //integrable->StateGather(X, V, systemdT);
        //double dt = try_realtime ? real_timer.SuggestSimulationStep(systemdT) : systemdT;
        auto dt = try_realtime ? real_timer.SuggestSimulationStep(m_system->GetStep()) : m_system->GetStep();
        current_time += dt;
        std::cout << dt << " s" << std::endl;

        assert(current_selected_mode <= eig_val.GetRows());
            
        // evaluate the state for the selected_mode a thte current timestep
        eig_vect_col.PasteClippedMatrix(&eig_vect, 0, current_selected_mode -1, eig_vect.GetRows(), 1, 0, 0);
        double pulse = sqrt(eig_val(current_selected_mode -1));
        double amplitude_corrector = magnitude_amplification*sin(time_amplification*pulse*current_time);
        //std::cout << amplitude_corrector << " ampl" << std::endl;
        //GetLog() << eig_vect_col << "\n";

        eig_vect_col.MatrScale(amplitude_corrector);

        // update the state
        ChStateDelta Dx;
        Dx.Reset(integrable->GetNcoords_v(), integrable);
        Dx.CopyFromMatrix(eig_vect_col);
        auto prova = X0.GetIntegrable();
        X = X0 + Dx;

        V.FillElem(0); // set V speed to zero
        integrable->StateScatter(X, V, dt);  // state -> system

    }

    void SetVisualizedMode(int mode_selection) { current_selected_mode = mode_selection; }
    int GetVisualizedMode() const { return current_selected_mode; }

    double GetVisualizedFrequency() const { return sqrt(eig_val(current_selected_mode - 1)) / 2 / CH_C_PI; }
    double GetFrequency(int mode) const { return sqrt(eig_val(mode - 1)) / 2 / CH_C_PI; }

    void SetRealtime(bool val) { try_realtime = val; }

    void SetTimeAmplification(double time_ampl) { time_amplification = time_ampl; }
    double GetTimeAmplification() const { return time_amplification; }
    void SetMagnitudeAmplification(double magnitude_ampl) { magnitude_amplification = magnitude_ampl; }
    double GetMagnitudeAmplification() const { return magnitude_amplification; }

    int GetComputedModes() const { return eig_val.GetRows(); }



#ifdef CHRONO_IRRLICHT
    void ActivateModalAnalysisGUI(irrlicht::ChIrrApp& application)
    {
        // ..Finally create the event receiver, for handling all the GUI (user will use
        //   buttons/sliders to modify parameters)
        gui_receiver.Initialize(application, *this);
        // note how to add the custom event receiver to the default interface:
        application.SetUserEventReceiver(&gui_receiver);
    }
#endif

};


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
        auto box_width = 50;
        auto item_height = 15;
        auto cum_y_pos = 35;

        mdevice->getGUIEnvironment()->addStaticText(L"Modes timestep", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_timestep_box = mdevice->getGUIEnvironment()->addEditBox(L"",
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width, cum_y_pos + item_height), true, mode_tab, modes_timestep_ID); cum_y_pos += item_height;

        mdevice->getGUIEnvironment()->addStaticText(L"Mode selected", 
            irr::core::rect<irr::s32>(left_margin, cum_y_pos, left_margin + text_width, cum_y_pos + item_height), false, true, mode_tab);
        auto mode_sel_box = mdevice->getGUIEnvironment()->addEditBox(L"", 
            irr::core::rect<irr::s32>(left_margin + text_width, cum_y_pos, left_margin + text_width + box_width*0.8, cum_y_pos + item_height), true, mode_tab, modes_selection_ID); cum_y_pos += item_height;
        mode_selmax_text = mdevice->getGUIEnvironment()->addStaticText(L"/",
                                                    irr::core::rect<irr::s32>(left_margin + text_width + box_width*0.8, cum_y_pos, left_margin + text_width + box_width*1.0, cum_y_pos + item_height), false, true, mode_tab);


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

        mode_timestep_box->setText(std::to_wstring(app->GetSystem()->GetStep()).c_str());
        mode_sel_box->setText(std::to_wstring(eig_analysis->GetVisualizedMode()).c_str());
        mode_magn_ampl_box->setText(std::to_wstring(eig_analysis->GetMagnitudeAmplification()).c_str());
        mode_time_ampl_box->setText(std::to_wstring(eig_analysis->GetTimeAmplification()).c_str());
        mode_freq_text->setText(std::wstring(L"Frequency: ").append(std::to_wstring(eig_analysis->GetVisualizedFrequency())).append(L"Hz").c_str());
        mode_selmax_text->setText(std::wstring(L"/").append(std::to_wstring(eig_analysis->GetComputedModes())).c_str());

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
                            eig_analysis->SetVisualizedMode(mode);
                            mode_selmax_text->setText(std::wstring(L"/").append(std::to_wstring(eig_analysis->GetComputedModes())).c_str());
                            mode_freq_text->setText(std::wstring(L"Frequency: ").append(std::to_wstring(eig_analysis->GetVisualizedFrequency())).append(L"Hz").c_str());
                            break;
                        }

                        case modes_magn_ampl_ID:
                        {
                            double magn_ampl = atof(irr::core::stringc(static_cast<irr::gui::IGUIEditBox*>(event.GUIEvent.Caller)->getText()).c_str());
                            eig_analysis->SetMagnitudeAmplification(magn_ampl);
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


}  // END_OF_NAMESPACE____
#endif
