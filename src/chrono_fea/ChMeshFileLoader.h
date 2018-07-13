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
// Authors: Andrea Favali, Alessandro Tasora
// =============================================================================
// Utilities for loading meshes from file
// =============================================================================

#ifndef CHMESH_FILE_LOADER_H
#define CHMESH_FILE_LOADER_H

#include <map>

#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChMesh.h"

#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChElementTetra_4.h"
#include "chrono_fea/ChMeshFileLoader.h"
#include "chrono_fea/ChNodeFEAxyz.h"
#include "ChMaterialShellReissner.h"
#include "ChNodeFEAxyzrot.h"
#include "ChElementShellReissner4.h"

namespace chrono {
namespace fea {

/// @addtogroup fea_module
/// @{

/// Collection of mesh file loader utilities.
class ChApiFea ChMeshFileLoader {
  public:
    /// Load tetrahedrons from .node and .ele files as saved by TetGen.
    /// The file format for .node (with point# starting from 1) is:
    ///   [# of points] [dimension (only 3)] [# of attributes (only 0)] [markers (only 0)]
    ///   [node #] [x] [y] [z]
    ///   [node #] [x] [y] [z]   etc.
    /// The file format for .ele (with tet# starting from 1) is:
    ///   [# of tetrahedrons] [dimension (only 4 supported)] [# of attributes (only 0)]
    ///   [tet #] [node #] [node #] [node #] [node #]
    ///   [tet #] [node #] [node #] [node #] [node #]   etc.
    /// If you pass a material inherited by ChContinuumElastic, nodes with 3D motion are used, and corotational
    /// elements.
    /// If you pass a material inherited by ChContinuumPoisson3D, nodes with scalar field are used (ex. thermal,
    /// electrostatics, etc)
    static void FromTetGenFile(
        std::shared_ptr<ChMesh> mesh,                      ///< destination mesh
        const char* filename_node,                         ///< name of the .node file
        const char* filename_ele,                          ///< name of the .ele  file
        std::shared_ptr<ChContinuumMaterial> my_material,  ///< material for the created tetahedrons
        ChVector<> pos_transform = VNULL,                  ///< optional displacement of imported mesh
        ChMatrix33<> rot_transform = ChMatrix33<>(1)       ///< optional rotation/scaling of imported mesh
    );

    static void FromAbaqusFileMOD(std::string filename,
        std::map<unsigned, std::tuple<std::string, std::vector<unsigned>>>& elements_map,
        std::map<unsigned, std::vector<double>>& nodes_map,
        std::map<std::string, std::vector<unsigned int>>& nset_map,
        std::map<std::string, std::vector<unsigned int>>& elset_map);

    /// Load tetrahedrons, if any, saved in a .inp file for Abaqus.
    template <class material_type>
    static void FromAbaqusFile(
        std::shared_ptr<ChMesh> mesh,                      ///< destination mesh
        const char* filename,                              ///< input file name
        std::shared_ptr<material_type> my_material,  ///< material for the created tetahedrons
        std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase> > >&
            node_sets,                                 ///< vect of vectors of 'marked'nodes
        ChVector<> pos_transform = VNULL,              ///< optional displacement of imported mesh
        ChMatrix33<> rot_transform = ChMatrix33<>(1),  ///< optional rotation/scaling of imported mesh
        bool discard_unused_nodes =
            true  ///< if true, Abaqus nodes that are not used in elements or sets are not imported in C::E
    ) {

        std::map<unsigned int, std::pair < std::shared_ptr<ChNodeFEAbase>, bool >> parsed_nodes;
        std::vector<std::shared_ptr<ChNodeFEAbase>>* current_nodeset_vector = nullptr;

        enum eChAbaqusParserSection {
            E_PARSE_UNKNOWN = 0,
            E_PARSE_NODES_XYZ,
            E_PARSE_TETS_4,
            E_PARSE_TETS_10,
            E_PARSE_CPS4,
            E_PARSE_NODESET
        } e_parse_section = E_PARSE_UNKNOWN;

        std::ifstream fin(filename);
        if (fin.good())
            GetLog() << "Parsing Abaqus INP file: " << filename << "\n";
        else
            throw ChException("ERROR opening Abaqus .inp file: " + std::string(filename) + "\n");

        std::string line;
        while (getline(fin, line)) {
            // trims white space from the beginning of the std::string
            line.erase(line.begin(), find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(isspace))));
            // convert parsed line to uppercase (since std::string::find is case sensitive and Abaqus INP is not)
            std::for_each(line.begin(), line.end(), [](char& c) { c = toupper(static_cast<unsigned char>(c)); });

            // skip empty lines
            if (line[0] == 0)
                continue;

            // check if the current line opens a new section
            if (line[0] == '*') {
                e_parse_section = E_PARSE_UNKNOWN;

                if (line.find("*NODE") == 0) {
                    std::string::size_type nse = line.find("NSET=");
                    if (nse > 0) {
                        std::string::size_type ncom = line.find(",", nse);
                        std::string s_node_set = line.substr(nse + 5, ncom - (nse + 5));
                        GetLog() << "| parsing nodes " << s_node_set << "\n";
                    }
                    e_parse_section = E_PARSE_NODES_XYZ;
                }

                if (line.find("*ELEMENT") == 0) {
                    std::string::size_type nty = line.find("TYPE=");
                    if (nty > 0) {
                        std::string::size_type ncom = line.find(",", nty);
                        std::string s_ele_type = line.substr(nty + 5, ncom - (nty + 5));
                        e_parse_section = E_PARSE_UNKNOWN;
                        if (s_ele_type == "C3D10") {
                            e_parse_section = E_PARSE_TETS_10;
                        }
                        else if (s_ele_type == "DC3D10") {
                            e_parse_section = E_PARSE_TETS_10;
                        }
                        else if (s_ele_type == "C3D4") {
                            e_parse_section = E_PARSE_TETS_4;
                        }
                        else if (s_ele_type == "CPS4") {
                            e_parse_section = E_PARSE_CPS4;
                        }
                        if (e_parse_section == E_PARSE_UNKNOWN) {
                            GetLog() << "| WARNING: cannot import element TYPE=" << s_ele_type << " at line " << line << "\n";
                        }
                        else
                        {
                            std::string::size_type nse = line.find("ELSET=");
                            if (nse > 0) {
                                std::string::size_type ncom = line.find(",", nse);
                                std::string s_ele_set = line.substr(nse + 6, ncom - (nse + 6));
                                GetLog() << "| parsing element set: " << s_ele_set << "\n";
                            }
                        }
                    }

                }

                if (line.find("*NSET") == 0) {
                    std::string::size_type nse = line.find("NSET=", 5);
                    if (nse > 0) {
                        std::string::size_type ncom = line.find(",", nse);
                        std::string s_node_set = line.substr(nse + 5, ncom - (nse + 5));
                        GetLog() << "| parsing nodeset: " << s_node_set << "\n";
                        auto new_node = node_sets.insert(std::pair < std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>>(
                            s_node_set, std::vector<std::shared_ptr<ChNodeFEAbase>>()));
                        if (new_node.second)
                        {
                            current_nodeset_vector = &new_node.first->second;
                        }
                        else
                            throw ChException("ERROR in .inp file, multiple NSET with same name has been specified\n");

                    }
                    e_parse_section = E_PARSE_NODESET;
                }

                continue;
            }

            // node parsing
            if (e_parse_section == E_PARSE_NODES_XYZ) {
                int idnode = 0;
                double x = -10e30;
                double y = -10e30;
                double z = -10e30;
                double tokenvals[20];
                int ntoken = 0;

                std::string token;
                std::istringstream ss(line);
                while (getline(ss, token, ',') && ntoken < 20) {
                    std::istringstream stoken(token);
                    stoken >> tokenvals[ntoken];
                    ++ntoken;
                }

                if (ntoken != 4)
                    throw ChException("ERROR in .inp file, nodes require ID and three x y z coords, see line:\n" + line +
                        "\n");

                x = tokenvals[1];
                y = tokenvals[2];
                z = tokenvals[3];
                if (x == -10e30 || y == -10e30 || z == -10e30)
                    throw ChException("ERROR in .inp file, in parsing x,y,z coordinates of node: \n" + line + "\n");

                // TODO: is it worth to keep a so specific routine inside this function?
                // especially considering that is affecting only some types of elements...
                ChVector<> node_position(x, y, z);
                node_position = rot_transform * node_position;  // rotate/scale, if needed
                node_position = pos_transform + node_position;  // move, if needed

                idnode = static_cast<unsigned int>(tokenvals[0]);
                if (std::dynamic_pointer_cast<ChContinuumElastic>(my_material)) {
                    auto mnode = std::make_shared<ChNodeFEAxyz>(node_position);
                    mnode->SetIndex(idnode);
                    parsed_nodes[idnode] = std::make_pair(mnode, false);
                    if (!discard_unused_nodes)
                        mesh->AddNode(mnode);
                }
                else if (std::dynamic_pointer_cast<ChContinuumPoisson3D>(my_material)) {
                    auto mnode = std::make_shared<ChNodeFEAxyzP>(ChVector<>(x, y, z));
                    mnode->SetIndex(idnode);
                    parsed_nodes[idnode] = std::make_pair(mnode, false);
                    if (!discard_unused_nodes)
                        mesh->AddNode(mnode);
                }
                else if (std::dynamic_pointer_cast<ChMaterialShellReissner>(my_material)) {
                    auto mnode = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector<>(x, y, z), QUNIT));
                    mnode->SetIndex(idnode);
                    parsed_nodes[idnode] = std::make_pair(mnode, false);
                    if (!discard_unused_nodes)
                        mesh->AddNode(mnode);
                }
                else
                    throw ChException("ERROR in .inp generation. Material type not supported. \n");
            }

            // element parsing
            if (e_parse_section == E_PARSE_TETS_10 || e_parse_section == E_PARSE_TETS_4) {
                int idelem = 0;
                unsigned int tokenvals[20];
                int ntoken = 0;

                std::string token;
                std::istringstream ss(line);
                while (std::getline(ss, token, ',') && ntoken < 20) {
                    std::istringstream stoken(token);
                    stoken >> tokenvals[ntoken];
                    ++ntoken;
                }

                if (e_parse_section == E_PARSE_TETS_10) {
                    if (ntoken != 11)
                        throw ChException("ERROR in .inp file, tetrahedrons require ID and 10 node IDs, see line:\n" +
                            line + "\n");
                    // TODO: 'idelem' might be stored in an index in ChElementBase in order to provide consistency with the
                    // INP file
                    idelem = (int)tokenvals[0];

                    for (int in = 0; in < 10; ++in)
                        if (tokenvals[in + 1] == -10e30)
                            throw ChException("ERROR in in .inp file, in parsing IDs of tetrahedron: \n" + line + "\n");
                }
                else if (e_parse_section == E_PARSE_TETS_4) {
                    if (ntoken != 5)
                        throw ChException("ERROR in .inp file, tetrahedrons require ID and 4 node IDs, see line:\n" + line +
                            "\n");

                    // TODO: 'idelem' might be stored in an index in ChElementBase in order to provide consistency with the
                    // INP file
                    idelem = (int)tokenvals[0];

                    for (int in = 0; in < 4; ++in)
                        if (tokenvals[in + 1] == -10e30)
                            throw ChException("ERROR in in .inp file, in parsing IDs of tetrahedron: \n" + line + "\n");
                }

                if (std::dynamic_pointer_cast<ChContinuumElastic>(my_material) ||
                    std::dynamic_pointer_cast<ChContinuumPoisson3D>(my_material)) {
                    std::array<std::shared_ptr<ChNodeFEAbase>, 4> element_nodes;
                    for (auto node_sel = 0; node_sel < 4; ++node_sel) {
                        // check if the nodes required by the current element exist
                        std::pair<std::shared_ptr<ChNodeFEAbase>, bool>& node_found =
                            parsed_nodes.at(static_cast<unsigned int>(tokenvals[node_sel + 1]));

                        element_nodes[node_sel] = node_found.first;
                        node_found.second = true;
                    }

                    if (std::dynamic_pointer_cast<ChContinuumElastic>(my_material)) {
                        auto mel = std::make_shared<ChElementTetra_4>();
                        mel->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(element_nodes[3]),
                            std::dynamic_pointer_cast<ChNodeFEAxyz>(element_nodes[1]),
                            std::dynamic_pointer_cast<ChNodeFEAxyz>(element_nodes[2]),
                            std::dynamic_pointer_cast<ChNodeFEAxyz>(element_nodes[0]));
                        mel->SetMaterial(std::dynamic_pointer_cast<ChContinuumElastic>(my_material));
                        mesh->AddElement(mel);

                    }
                    else if (std::dynamic_pointer_cast<ChContinuumPoisson3D>(my_material)) {
                        auto mel = std::make_shared<ChElementTetra_4_P>();
                        mel->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzP>(element_nodes[0]),
                            std::dynamic_pointer_cast<ChNodeFEAxyzP>(element_nodes[1]),
                            std::dynamic_pointer_cast<ChNodeFEAxyzP>(element_nodes[2]),
                            std::dynamic_pointer_cast<ChNodeFEAxyzP>(element_nodes[3]));
                        mel->SetMaterial(std::dynamic_pointer_cast<ChContinuumPoisson3D>(my_material));
                        mesh->AddElement(mel);
                    }

                }
                else
                    throw ChException("ERROR in .inp generation. Material type not supported.\n");
            }

            if (e_parse_section == E_PARSE_CPS4) {
                int idelem = 0;
                unsigned int tokenvals[20];
                int ntoken = 0;

                std::string token;
                std::istringstream ss(line);
                while (std::getline(ss, token, ',') && ntoken < 20) {
                    std::istringstream stoken(token);
                    stoken >> tokenvals[ntoken];
                    ++ntoken;
                }


                if (ntoken != 5)
                    throw ChException("ERROR in .inp file, tetrahedrons require ID and 4 node IDs, see line:\n" + line +
                        "\n");

                // TODO: 'idelem' might be stored in an index in ChElementBase in order to provide consistency with the
                // INP file
                idelem = static_cast<int>(tokenvals[0]);

                for (int in = 0; in < 4; ++in)
                    if (tokenvals[in + 1] == -10e30)
                        throw ChException("ERROR in in .inp file, in parsing IDs of shells: \n" + line + "\n");


                if (std::dynamic_pointer_cast<ChMaterialShellReissner>(my_material)) {
                    std::array<std::shared_ptr<ChNodeFEAbase>, 4> element_nodes;
                    for (auto node_sel = 0; node_sel < 4; ++node_sel) {
                        // check if the nodes required by the current element exist
                        std::pair<std::shared_ptr<ChNodeFEAbase>, bool>& node_found =
                            parsed_nodes.at(static_cast<unsigned int>(tokenvals[node_sel + 1]));

                        element_nodes[node_sel] = node_found.first;
                        node_found.second = true;
                    }

                    auto mel = std::make_shared<ChElementShellReissner4>();
                    mel->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzrot>(element_nodes[0]),
                        std::dynamic_pointer_cast<ChNodeFEAxyzrot>(element_nodes[1]),
                        std::dynamic_pointer_cast<ChNodeFEAxyzrot>(element_nodes[2]),
                        std::dynamic_pointer_cast<ChNodeFEAxyzrot>(element_nodes[3]));
                    auto thickness = 0.1;
                    mel->AddLayer(thickness, 0 * CH_C_DEG_TO_RAD, std::dynamic_pointer_cast<ChMaterialShellReissner>(my_material));
                    mel->SetAlphaDamp(0.0);
                    mesh->AddElement(mel);

                }
                else
                    throw ChException("ERROR in .inp generation. Material type not supported.\n");
            }

            // parsing nodesets
            if (e_parse_section == E_PARSE_NODESET) {
                unsigned int tokenvals[20]; // strictly speaking, the maximum is 16 nodes for each line
                int ntoken = 0;

                std::string token;
                std::istringstream ss(line);
                while (std::getline(ss, token, ',') && ntoken < 20) {
                    std::istringstream stoken(token);
                    stoken >> tokenvals[ntoken];
                    ++ntoken;
                }

                for (auto node_sel = 0; node_sel < ntoken; ++node_sel) {
                    auto idnode = static_cast<int>(tokenvals[node_sel]);
                    if (idnode > 0) {
                        // check if the nodeset is asking for an existing node
                        std::pair<std::shared_ptr<ChNodeFEAbase>, bool>& node_found =
                            parsed_nodes.at(static_cast<unsigned int>(tokenvals[node_sel]));

                        current_nodeset_vector->push_back(node_found.first);
                        // flag the node to be saved later into the mesh
                        node_found.second = true;

                    }
                    else
                        throw ChException("ERROR in .inp file, negative node ID: " + tokenvals[node_sel]);
                }
            }

        }  // end while

           // only used nodes have been saved in 'parsed_nodes_used' and are now inserted in the mesh node list
        if (discard_unused_nodes) {
            for (auto node_it = parsed_nodes.begin(); node_it != parsed_nodes.end(); ++node_it) {
                if (node_it->second.second)
                    mesh->AddNode(node_it->second.first);
            }
        }
    }

    static void ANCFShellFromGMFFile(
        std::shared_ptr<ChMesh> mesh,                      ///< destination mesh
        const char* filename,                              ///< complete filename
        std::shared_ptr<ChMaterialShellANCF> my_material,  ///< material to be given to the shell
        std::vector<double>& node_ave_area,                ///< output the average area of the nodes
        std::vector<int>& BC_nodes,                        ///< material to be given to the shell
        ChVector<> pos_transform = VNULL,                  ///< optional displacement of imported mesh
        ChMatrix33<> rot_transform = ChMatrix33<>(1),      ///< optional rotation/scaling of imported mesh
        double scaleFactor = 1,                            ///< import scale factor
        bool printNodes = false,                           ///< display the imported nodes
        bool printElements = false                         ///< display the imported elements
    );
};

/// @} fea_module

}  // end namespace fea
}  // end namespace chrono

#endif
