// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Irrlicht-based visualization wrapper for track test rig.
// This class extends ChVehicleIrrApp.
//
// =============================================================================

#include "chrono/core/ChMathematics.h"

#include "chrono_vehicle/tracked_vehicle/utils/ChTrackTestRigIrrApp.h"

using namespace irr;

namespace chrono {
namespace vehicle {

ChTrackTestRigIrrApp::ChTrackTestRigIrrApp(ChTrackTestRig* rig,
                                           const std::wstring& title,
                                           const irr::core::dimension2d<irr::u32>& dims,
                                           irr::ELOG_LEVEL log_level)
    : ChVehicleIrrApp(rig, title, dims, log_level), m_rig(rig) {}

// Render contact normals for monitored subsystems
void ChTrackTestRigIrrApp::renderOtherGraphics() {
    bool normals = m_rig->m_contact_manager->m_render_normals;
    bool forces = m_rig->m_contact_manager->m_render_forces;
    double scale_normals = 0.4;
    double scale_forces = m_rig->m_contact_manager->m_scale_forces;

    // Contact normals on left sprocket.
    // Note that we only render information for contacts on the outside gear profile
    for (const auto& c : m_rig->m_contact_manager->m_sprocket_L_contacts) {
        ChVector<> v1 = c.m_point;
        if (normals) {
        ChVector<> v2 = v1 + c.m_csys.Get_A_Xaxis() * scale_normals;
        if (v1.y() > m_rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos().y())
            irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, video::SColor(255, 180, 0, 0), false);
        }
        if (forces) {
            ChVector<> v2 = v1 + c.m_force * scale_forces;
            if (v1.y() > m_rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos().y())
                irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, video::SColor(255, 180, 0, 0), false);
        }
    }

    // Contact normals on rear sprocket.
    // Note that we only render information for contacts on the outside gear profile
    for (const auto& c : m_rig->m_contact_manager->m_sprocket_R_contacts) {
        ChVector<> v1 = c.m_point;
        if (normals) {
        ChVector<> v2 = v1 + c.m_csys.Get_A_Xaxis() * scale_normals;
        if (v1.y() < m_rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos().y())
            irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, video::SColor(255, 180, 0, 0), false);
        }
        if (forces) {
            ChVector<> v2 = v1 + c.m_force * scale_forces;
            if (v1.y() > m_rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos().y())
                irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, video::SColor(255, 180, 0, 0), false);
        }
    }

    // Contact normals on monitored track shoes.
    renderContacts(m_rig->m_contact_manager->m_shoe_L_contacts, video::SColor(255, 180, 180, 0), normals, forces,
                   scale_normals, scale_forces);
    renderContacts(m_rig->m_contact_manager->m_shoe_R_contacts, video::SColor(255, 180, 180, 0), normals, forces,
                   scale_normals, scale_forces);

    // Contact normals on idler wheels.
    renderContacts(m_rig->m_contact_manager->m_idler_L_contacts, video::SColor(255, 0, 0, 180), normals, forces,
                   scale_normals, scale_forces);
    renderContacts(m_rig->m_contact_manager->m_idler_R_contacts, video::SColor(255, 0, 0, 180), normals, forces,
                   scale_normals, scale_forces);
}

// Render normal for all contacts in the specified list, using the given color.
void ChTrackTestRigIrrApp::renderContacts(const std::list<ChTrackContactManager::ContactInfo>& lst,
                                          const video::SColor& col,
                                          bool normals,
                                          bool forces,
                                          double scale_normals,
                                          double scale_forces) {
    for (const auto& c : lst) {
        ChVector<> v1 = c.m_point;
        if (normals) {
            ChVector<> v2 = v1 + c.m_csys.Get_A_Xaxis() * scale_normals;
            irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, col, false);
        }
        if (forces) {
            ChVector<> v2 = v1 + c.m_force * scale_forces;
            irrlicht::tools::drawSegment(GetVideoDriver(), v1, v2, video::SColor(255, 180, 0, 0), false);
        }
    }
}

}  // end namespace vehicle
}  // end namespace chrono