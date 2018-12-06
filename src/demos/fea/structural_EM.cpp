// =============================================================================
//
// FEA analysis for electric motors
//
// =============================================================================

#include <vector>
#include "chrono/physics/ChSystemNSC.h"


#include <set>
#include "chrono_fea/ChElementHexa_8.h"
#include "chrono_fea/ChElementShellReissner4.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChMeshFileLoader.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_mkl/ChSolverMKL.h"
#include "utils/ChUtilsInputOutput.h"


using namespace chrono;
using namespace chrono::fea;

#define USE_MKL
//#define USE_IRRLICHT
#define FULL_STRESS_OUTPUT
//#define EQUAL_ELEMENT_SPACING
#define LOG_OUTPUT true


#ifdef USE_IRRLICHT 
#include "chrono_irrlicht/ChIrrApp.h"
using namespace chrono::irrlicht;
using namespace irr;
#endif

class CSVwriter
{
private:
    std::ofstream myfile;
    std::ostringstream outbuffer;
    bool enabled = true;
    char delim = ',';
public:

    explicit CSVwriter(const std::string& filename, bool enabled = true, char delimiter = ',') : enabled(enabled), delim(delimiter)
    {
        if (myfile.is_open())
            myfile.close();

        if (enabled)
        {
            myfile.open(filename, std::ios_base::out);
            if (!myfile.good())
                throw ChException("File with name: " + filename + "cannot be found.");
        }

    }

    template<typename T, typename... args_t>
    void AppendRow(const T& objects, const args_t&... objects_other)
    {
        if (!enabled)
            return;

        outbuffer << objects << delim;
        this->AppendRow(objects_other...);
    }

    template<typename T>
    void AppendRow(const T& objects)
    {
        if (!enabled)
            return;

        outbuffer << objects;
        myfile << outbuffer.str() << std::endl;
        outbuffer.str("");
        outbuffer.clear();
    }

    void EnableWriting(bool on_off) { enabled = on_off;}


    ~CSVwriter()
    {
        myfile.close();
    }
};

class CSVreader
{
private:
    std::ifstream myfile;
    char delim = ',';
public:

    CSVreader() {}

    explicit CSVreader(const std::string& filename, char delim = ',')
    {
        SetFile(filename, delim);
    }

    void SetFile(const std::string& filename, char delimiter = ',')
    {
        delim = delimiter;

        if (myfile.is_open())
            myfile.close();

        myfile.open(filename, std::ios_base::in);
        if (!myfile.good())
            throw ChException("File with name: " + filename + "cannot be found.");
    }

    template <typename type_t>
    void ParseRow(std::vector<type_t>& vector_out, unsigned int rownum, unsigned int skip_columns = 0)
    {
        GoToLine(rownum);
        this->ParseCurrentRow(vector_out, skip_columns);
    }


    template <typename type_t>
    void ParseCurrentRow(std::vector<type_t>& vector_out, unsigned int skip_columns = 0)
    {
        std::string linepiece;
        std::string line;
        std::getline(myfile, line); // gets line until \n
        std::istringstream line_ss(line);
        type_t val;

        unsigned int current_col = 0;
        while (std::getline(line_ss, linepiece, delim)) {

            if (current_col < skip_columns)
            {
                ++current_col;
                continue;
            }

            std::istringstream stoken(linepiece);
            stoken >> val;
            vector_out.push_back(val);
            ++current_col;

        }
    }

    void GoToLine(unsigned int rownum)
    {
        myfile.seekg(std::ios::beg);
        for (unsigned int line_sel = 0; line_sel < rownum; ++line_sel) {
            myfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }

    template <typename type_t>
    bool ParseElement(type_t& val_out, unsigned int rownum, unsigned int colnum)
    {
        GoToLine(rownum);

        std::string linepiece;
        std::string line;
        std::getline(myfile, line, delim);
        std::istringstream line_ss(line);

        unsigned int current_col = 0;
        while (std::getline(line_ss, linepiece, delim)) {
            if (current_col == colnum)
            {
                std::istringstream stoken(linepiece);
                stoken >> val_out;
                return true;
            }
            ++current_col;

        }

        return false;
    }

    template <typename type_t>
    void ParseFile(std::list<std::vector<type_t>>& vector_out, unsigned int skip_rows = 0, unsigned int skip_columns = 0)
    {

        GoToLine(skip_rows);

        // create vector to store values for the current line
        while(!myfile.eof())
        {
            auto current_vector = vector_out.insert(vector_out.end(), std::vector<type_t>());
            this->ParseCurrentRow(vector_out, skip_columns);
        }

    }

    ~CSVreader()
    {
        myfile.close();
    }
};


void getSigmaGlob(const std::vector<ChVector<>>& sigma_glob_set, ChVector<>& sigma_glob, double angle)
{
    angle = fmod(angle, CH_C_2PI);
    if (angle < 0)
        angle += CH_C_2PI;

    assert(angle <= CH_C_2PI);
    assert(angle >= 0);

    double index = angle / CH_C_2PI * (sigma_glob_set.size()-1);
    auto index_int = static_cast<size_t>(floor(index));
    //int next_index = index_int > sigma_glob_set.size() ? index_int - sigma_glob_set.size() : index_int;
    auto next_index = std::min(index_int + 1, sigma_glob_set.size()-1);

    sigma_glob = sigma_glob_set[index_int] + (index - index_int) * (sigma_glob_set[next_index] - sigma_glob_set[index_int]);
    
}

void getArcLength(const std::set<double>& angles, double& length_previous, double& length_next, double angle, double radius)
{
    angle = fmod(angle, CH_C_2PI);
    if (angle < 0)
        angle += CH_C_2PI;

    assert(angle <= CH_C_2PI);
    assert(angle >= 0);

    auto iter_prev = angles.lower_bound(angle);
    auto iter_next = iter_prev;
    double angle_previous = iter_prev != angles.begin() ? *--iter_prev : *(--angles.end()) - CH_C_2PI;
    double angle_next = (iter_next != angles.end()) && ++iter_next != angles.end() ? *iter_next : *(angles.begin())+CH_C_2PI;

    length_previous = 0.5*radius*(angle - angle_previous);
    length_next = 0.5*radius*(angle_next - angle);

    // OVERRIDE
    //length_next = 0.5*CH_C_2PI * 80.22e-3 / (560.0);
    //length_previous = length_next;

    assert(length_previous >= 0.0);
    assert(length_next >= 0.0);

}


int main(int argc, char* argv[]) {

    ////////////// Parse input arguments //////////////
    std::string meshfilepath;
    std::string sigmafilepath;
    std::string resultspath;
    std::string motor_prefix;
    std::string drivingcyle_prefix;
    unsigned int testindex_first;
    unsigned int testindex_last;

    if (argc>1)
    {
        meshfilepath = argv[1];
        sigmafilepath = argv[2];
        resultspath = argv[3];
        motor_prefix = argv[4];
        drivingcyle_prefix = argv[5];
        std::stringstream testindex_first_ss(argv[6]);
        testindex_first_ss >> testindex_first;
        std::stringstream testindex_last_ss(argv[7]);
        testindex_last_ss >> testindex_last;

        GetLog() << "Running tests from: " << testindex_first << " to " << testindex_last << "\n";
    }
    else
    {
        GetLog() << "structural_EM has to be called with the following syntax:\n";
        GetLog() << argv[0] << " <meshfolder> <stressEMfolder> <resultsfolder> <motor_name> <driving_cycle_name> <testindex_start> <testindex_end>\n";
        GetLog() << "where:\n";
        GetLog() << "<meshfolder> must point to the folder that contains <motor_name>_mesh.INP\n";
        GetLog() << "<stressEMfolder> must point to the folder that contains 'sigma_n' and 'sigma_t' files with the following naming convention:\n";
        GetLog() << "   <motor_name>_sigma_n.csv and <motor_name>_sigma_t.csv\n";
        GetLog() << "<resultsfolder> is the output folder in which results are stored;\n";
        GetLog() << "<motor_name> is the name of the simulated motor;\n";
        GetLog() << "<driving_cycle_name> is the name of the simulated driving cycle;\n";
        GetLog() << "<testindex_start> and <testindex_end> specify the range of test that has to be run\n";
        GetLog() << "    the testindex value refers to the row (of 'sigma_n' and 'sigma_t' files) that holds information about the current test (zero-based);\n";
        GetLog() << "Please mind that folders must end with a / sign. Just put double double-quotes to specify current folder.\n";
        GetLog() << "\n";
        GetLog() << "<motor_name>_sigma_n.csv (and similarly for sigma_t) has multiple lines, each of which\n";
        GetLog() << "    holds information about a specific working point and must have the following structure:\n";
        GetLog() << "    | angularspeed [rpm] | torque [Nm] | rotor position | sigmas...[Pa] |\n";
        GetLog() << "\n";
        GetLog() << "The output will be in:\n";
        GetLog() << "<resultsfolder>/<motor_name>_<driving_cycle_name>_<testindex>_stress_allelements.csv\n";
        GetLog() << "with the following convention:\n";
        GetLog() << "| ElementID | Stress VonMises |\n";
        GetLog() << "<resultsfolder>/<motor_name>_<driving_cycle_name>_<testindex>_stress_bridges.csv\n";
        GetLog() << "<resultsfolder>/<motor_name>_<driving_cycle_name>_<testindex>_stress_forfatigue.csv\n";
        GetLog() << "with the following convention:\n";
        GetLog() << "| ElementID | Stress VonMises | StressXX | StressYY | StressZZ | StressXY | StressYZ | StressXZ |\n";
        GetLog() << "\n";

        return -1;
    }

    auto filename_mesh = meshfilepath + motor_prefix + "_mesh.INP";
    auto filename_meshinfo = meshfilepath + motor_prefix + "_meshinfo.csv";
    auto filename_sigma_n = sigmafilepath + motor_prefix + "_" + drivingcyle_prefix + "_sigma_n.csv";
    auto filename_sigma_t = sigmafilepath + motor_prefix + "_" + drivingcyle_prefix + "_sigma_t.csv";

    GetLog() << "The motor is " << motor_prefix << "\n";
    GetLog() << "The driving cycle is " << drivingcyle_prefix << "\n";
    GetLog() << "Required files are:\n" <<
        " - mesh file: " << filename_mesh << "\n"
        " - sigma_n file: " << filename_sigma_n << "\n"
        " - sigma_t file: " << filename_sigma_t << "\n";
    GetLog() << "Results will be stored in: "<< resultspath <<"\n\n";

    ////////////// Chrono setup //////////////
    ChTimer<> tim;
    // Create a Chrono::Engine physical system
    ChSystemNSC my_system;

#ifdef USE_MKL
    auto mkl_solver = std::make_shared<ChSolverMKL<>>();
    mkl_solver->SetSparsityPatternLock(true);
    mkl_solver->ForceSparsityPatternUpdate(true);
    my_system.SetSolver(mkl_solver);
#else
    my_system.SetSolverType(ChSolver::Type::MINRES);  // <- NEEDED THIS or Matlab or MKL solver
    my_system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(200);
    my_system.SetMaxItersSolverStab(200);
    my_system.SetTolForce(1e-13);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);
#endif

#ifdef USE_IRRLICHT
    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Shells FEA", core::dimension2d<u32>(1024, 768), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.2)),
                                 core::vector3dfCH(ChVector<>(0.0, 0.0, 0.0)));
    // application.SetContactsDrawMode(irr::ChIrrTools::CONTACT_DISTANCES);

    application.AddLightWithShadow(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.5)), core::vector3df(0.0, 0.0, 0.0), 0.1,
                                   0.0, 0.5, 0, 512, video::SColorf((f32)0.8, (f32)0.8, (f32)1.0));

#endif

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    my_mesh->SetAutomaticGravity(false);

    CSVreader csv_utility;

    ////////////// Create Material //////////////
    double rho = 7850;
    double E = 200e9;
    double nu = 0.3;
    auto element_thickness = 5e-3;
    auto element_material = std::make_shared<ChContinuumElastic>(E, nu, rho);


    ////////////// Load mesh from Abaqus file and identify nodes //////////////
    tim.start();
    std::map<unsigned, std::tuple<std::string, std::vector<unsigned>>> elements_map;
    std::map<unsigned, std::vector<double>> nodes_map;
    std::map<std::string, std::vector<unsigned int>> nset_map;
    std::map<std::string, std::vector<unsigned int>> elset_map;
    std::map<unsigned int, std::shared_ptr<ChNodeFEAxyz>> inserted_nodes;
    std::map<std::shared_ptr<ChElementHexa_8>, unsigned int> inserted_elements_ptr_to_ID;
    std::map<unsigned int, std::shared_ptr<ChElementHexa_8>> inserted_elements_ID_to_ptr;
    //std::map<std::shared_ptr<ChNodeFEAxyz>, unsigned int> inserted_nodes_ptr_to_ID;

    
    // parse Abaqus INP file
    try {
        ChMeshFileLoader::FromAbaqusFileMOD(filename_mesh, elements_map, nodes_map, nset_map, elset_map);
    }
    catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    // add node and elements
    std::string element_tag = "C3D8";
    for (auto el_it = elements_map.begin(); el_it != elements_map.end(); ++el_it) {
        if (std::get<0>(el_it->second) == element_tag) {
            auto new_elem = std::make_shared<ChElementHexa_8>();
            auto nodeid_vect = std::get<1>(el_it->second);
            std::array<std::shared_ptr<ChNodeFEAxyz>, 8> nodes;
            for (auto node_sel = 0; node_sel < 8; ++node_sel) {
                // check if the node specified by the current element exists
                auto node = nodes_map.find(nodeid_vect[node_sel]);
                if (node != nodes_map.end()) {
                    // check if the node specified by the current has not been inserted yet
                    auto node_found = inserted_nodes.find(nodeid_vect[node_sel]);
                    if (node_found == inserted_nodes.end()) {
                        nodes[node_sel] = std::make_shared<ChNodeFEAxyz>(ChVector<>(node->second[0], node->second[1], node->second[2]));
                        my_mesh->AddNode(nodes[node_sel]);
                        inserted_nodes.emplace_hint(inserted_nodes.end(), nodeid_vect[node_sel], nodes[node_sel]);
                        //inserted_nodes_ptr_to_ID.emplace_hint(inserted_nodes_ptr_to_ID.end(), nodes[node_sel], nodeid_vect[node_sel]);
                    }
                    else {
                        nodes[node_sel] = node_found->second;
                    }
                }
                else
                    throw ChException("Node not found\n");
            }
            new_elem->SetNodes(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
            new_elem->SetMaterial(element_material);
            
            my_mesh->AddElement(new_elem);
            inserted_elements_ptr_to_ID.emplace_hint(inserted_elements_ptr_to_ID.end(), new_elem, el_it->first);
            inserted_elements_ID_to_ptr.emplace_hint(inserted_elements_ID_to_ptr.end(), el_it->first, new_elem);

        }
    }

    GetLog() << "Added " << inserted_nodes.size() << " nodes over " << nodes_map.size() << ".\n";
    GetLog() << "Added " << my_mesh->GetElements().size() << " elements over " << elements_map.size() << ".\n";



    // Get magnets mass*radius info from Abaqus header section
    std::ifstream fin(filename_mesh);
    if (!fin.good())
        throw ChException("ERROR opening Abaqus INP file: " + std::string(filename_mesh) + "\n");

    // look for mass*radius value
    std::string line;
    double rotor_external_radius = -1e30;
    double rotor_internal_radius = -1e30;
    double magnets_mass_radius_outerext = -1e30;
    double magnets_mass_radius_outerint = -1e30;
    double magnets_mass_radius_median = -1e30;
    while (getline(fin, line)) {
        // trims white space from the beginning of the std::string
        line.erase(line.begin(), find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(isspace))));
        // convert parsed line to uppercase (since std::string::find is case sensitive and Abaqus INP is not)
        std::for_each(line.begin(), line.end(), [](char& c) { c = toupper(static_cast<unsigned char>(c)); });

        std::string string_to_find = "MAGNETS_MASS_RADIUS_OUTEREXT=";
        auto nse = line.find(string_to_find);
        if (nse != std::string::npos) {
            std::string::size_type ncom = line.find(",", nse + string_to_find.size());
            std::string magnets_mass_radius_s = line.substr(nse + string_to_find.size(), ncom - (nse + 5));
            std::stringstream magnets_mass_radius_ss(magnets_mass_radius_s);
            magnets_mass_radius_ss >> magnets_mass_radius_outerext;
            GetLog() << "OuterExt magnets mass*radius: " << magnets_mass_radius_outerext << "\n";
        }

        string_to_find = "MAGNETS_MASS_RADIUS_OUTERINT=";
        nse = line.find(string_to_find);
        if (nse != std::string::npos) {
            std::string::size_type ncom = line.find(",", nse + string_to_find.size());
            std::string magnets_mass_radius_s = line.substr(nse + string_to_find.size(), ncom - (nse + 5));
            std::stringstream magnets_mass_radius_ss(magnets_mass_radius_s);
            magnets_mass_radius_ss >> magnets_mass_radius_outerint;
            GetLog() << "OuterInt magnets mass*radius: " << magnets_mass_radius_outerint << "\n";
        }

        string_to_find = "MAGNETS_MASS_RADIUS_MEDIAN=";
        nse = line.find(string_to_find);
        if (nse != std::string::npos) {
            std::string::size_type ncom = line.find(",", nse + string_to_find.size());
            std::string magnets_mass_radius_s = line.substr(nse + string_to_find.size(), ncom - (nse + 5));
            std::stringstream magnets_mass_radius_ss(magnets_mass_radius_s);
            magnets_mass_radius_ss >> magnets_mass_radius_median;
            GetLog() << "Median magnets mass*radius: " << magnets_mass_radius_median << "\n";
        }

        string_to_find = "INTERNAL_RADIUS=";
        nse = line.find(string_to_find);
        if (nse != std::string::npos) {
            std::string::size_type ncom = line.find(",", nse + string_to_find.size());
            std::string rotor_internal_radius_s = line.substr(nse + string_to_find.size(), ncom - (nse + 5));
            std::stringstream rotor_internal_radius_ss(rotor_internal_radius_s);
            rotor_internal_radius_ss >> rotor_internal_radius;
            GetLog() << "Internal radius:  " << rotor_internal_radius << "\n";
        }

        string_to_find = "EXTERNAL_RADIUS=";
        nse = line.find(string_to_find);
        if (nse != std::string::npos) {
            std::string::size_type ncom = line.find(",", nse + string_to_find.size());
            std::string rotor_external_radius_s = line.substr(nse + string_to_find.size(), ncom - (nse + 5));
            std::stringstream rotor_external_radius_ss(rotor_external_radius_s);
            rotor_external_radius_ss >> rotor_external_radius;
            GetLog() << "External radius: " << rotor_external_radius << "\n";
        }

    }

    if (rotor_external_radius < 0 || rotor_internal_radius < 0 || magnets_mass_radius_outerext < 0 || magnets_mass_radius_median < 0)
    {
        GetLog() << "Radius info contained in mesh file cannot be loaded.\n";
        throw ChException("Radius info contained in mesh file cannot be loaded.");
    }


    // identify internal and external nodes
    //double rotor_external_radius = 80.22e-3;
    //double rotor_internal_radius = 25.5e-3;
    double external_threshold = rotor_external_radius - 1e-4;
    double internal_threshold = rotor_internal_radius + 1e-4;
    auto nodesmesh = my_mesh->GetNodes();
    std::vector<std::shared_ptr<ChNodeFEAxyz>> external_nodes;
    std::vector<std::shared_ptr<ChNodeFEAxyz>> internal_nodes;
    std::set<double> angles;
    for (auto node_sel = 0; node_sel < nodesmesh.size(); ++node_sel) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(nodesmesh[node_sel]);
        auto dist_from_center = sqrt(node->GetPos().x()*node->GetPos().x() + node->GetPos().y()*node->GetPos().y());
        if (dist_from_center > external_threshold) {
            external_nodes.push_back(node);
            double angle = atan2(node->GetPos().y(), node->GetPos().x());
            if (angle < 0)
                angle += CH_C_2PI;
            angles.insert(angle);
        }
        else if (dist_from_center < internal_threshold) {
            internal_nodes.push_back(node);
        }
    }

    // fix internal nodes
    for (auto it = internal_nodes.begin(); it != internal_nodes.end(); ++it) {
        (*it)->SetFixed(true);
    }

    GetLog() << "External nodes: " << external_nodes.size() << "\n";
    GetLog() << "Internal nodes: " << internal_nodes.size() << "\n";

    tim.stop();
    GetLog() << "Load and identify elements: " << tim() << "\n";

#ifdef USE_IRRLICHT

     //auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
     //mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
     //mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
     //mvisualizemesh->SetShrinkElements(true, 0.85);
     //mvisualizemesh->SetSmoothFaces(true);
     //my_mesh->AddAsset(mvisualizemesh);

    auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(false);
    mvisualizemeshref->SetDrawInUndeformedReference(false);
    my_mesh->AddAsset(mvisualizemeshref);

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRESS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_STRESS_VONMISES);
    mvisualizemeshC->SetSymbolsThickness(element_thickness);
    mvisualizemeshC->SetColorscaleMinMax(0,250e6);
    my_mesh->AddAsset(mvisualizemeshC);

    application.AssetBindAll();
    application.AssetUpdateAll();
#endif

    tim.start();
    my_system.SetupInitial();
    tim.stop();
    GetLog() << "SetupInitial time: " << tim() << "\n";


    ///////////////////////////////////////////////////////////////
    ///// Different working points affect only the code below /////
    ///////////////////////////////////////////////////////////////
    GetLog() << "\nTest Loop:\n";
    for (auto test_sel = testindex_first; test_sel <= testindex_last; ++test_sel) {
        tim.start();
        GetLog() << "Running test: " << test_sel << "\n";
        ////////////// Clear state of the nodes //////////////
        for (auto node_it = my_mesh->GetNodes().begin(); node_it != my_mesh->GetNodes().end(); ++node_it)
        {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(*node_it);
            node->SetPos(node->GetX0());
            node->SetForce(VNULL);
            node->SetNoSpeedNoAcceleration();
        }

        ////////////// Acquire EM pressure //////////////
        std::vector<double> sigma_t;
        std::vector<double> sigma_n;
        double omega;
        {
            const unsigned int OMEGA_POSITION_IN_CSV = 0;
            const unsigned int TORQUE_POSITION_IN_CSV = 1;
            const unsigned int POSROT_POSITION_IN_CSV = 2;
            double omega_n, omega_t, torque_n, torque_t, posrot_n, posrot_t;
            const unsigned int skip_columns = 3;
            csv_utility.SetFile(filename_sigma_n);
            csv_utility.ParseRow(sigma_n, test_sel, skip_columns);
            csv_utility.ParseElement(omega_n, test_sel, OMEGA_POSITION_IN_CSV);
            csv_utility.ParseElement(torque_n, test_sel, TORQUE_POSITION_IN_CSV);
            csv_utility.ParseElement(posrot_n, test_sel, POSROT_POSITION_IN_CSV);

            csv_utility.SetFile(filename_sigma_t);
            csv_utility.ParseRow(sigma_t, test_sel, skip_columns);
            csv_utility.ParseElement(omega_t, test_sel, OMEGA_POSITION_IN_CSV);
            csv_utility.ParseElement(torque_t, test_sel, TORQUE_POSITION_IN_CSV);
            csv_utility.ParseElement(posrot_t, test_sel, POSROT_POSITION_IN_CSV);

            if (abs(omega_n - omega_t) > 1e-12 || abs(torque_n - torque_t) > 1e-12 || abs(posrot_n - posrot_t) > 1e-12)
            {
                GetLog() << "Test with row index " << test_sel << " has been skipped because different values of omega, torque or posrot have been found between sigma files\n";
                continue;
            }

            omega = omega_n * CH_C_2PI / 60.0;
            //omega = omega_n;

            //GetLog() << "WARNING: omega 0\n";
            //omega = 0;
        }

        // store EM pressure (in global coordinates)
        CSVwriter sigma_glob_writer(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_sigma_glob.csv", LOG_OUTPUT);
        std::vector<ChVector<>> sigma_glob_set;
        sigma_glob_set.resize(sigma_n.size());
        auto delta_angle = CH_C_2PI / sigma_glob_set.size();
        for (auto sigma_sel = 0; sigma_sel < sigma_n.size(); ++sigma_sel) {

            sigma_glob_set[sigma_sel].x() = sigma_n[sigma_sel] * cos(delta_angle * sigma_sel) - sigma_t[sigma_sel] * sin(delta_angle * sigma_sel);
            sigma_glob_set[sigma_sel].y() = sigma_n[sigma_sel] * sin(delta_angle * sigma_sel) + sigma_t[sigma_sel] * cos(delta_angle * sigma_sel);
            sigma_glob_set[sigma_sel].z() = 0.0;

            sigma_glob_writer.AppendRow(delta_angle * sigma_sel,
                sigma_glob_set[sigma_sel].x(),
                sigma_glob_set[sigma_sel].y(),
                sigma_glob_set[sigma_sel].z());
        }

        ////////////// Apply additional centrifugal forces to emulate magnets //////////////
        {
            auto magnet_nodes = nset_map.find("MAGNETS_OUTEREXT");
            if (magnet_nodes == nset_map.end())
            {
                std::cout << "WARNING: MAGNETS_OUTEREXT nodeset not found." << std::endl;
            }
            double magnets_mass_radius_scattered = magnets_mass_radius_outerext / magnet_nodes->second.size();
            if (magnet_nodes != nset_map.end())
            {
                for (auto node_id_it = magnet_nodes->second.begin(); node_id_it != magnet_nodes->second.end(); ++node_id_it)
                {
                    auto node = inserted_nodes.at(*node_id_it);
                    auto centrifugal_force_vector = ChVector<>(node->GetPos().x(), node->GetPos().y(), 0.0);
                    centrifugal_force_vector.Normalize();
                    auto old_force = node->GetForce();
                    node->SetForce(old_force + magnets_mass_radius_scattered * omega*omega*centrifugal_force_vector);
                }
            }
        }

        {
            auto magnet_nodes = nset_map.find("MAGNETS_OUTERINT");
            if (magnet_nodes == nset_map.end())
            {
                std::cout << "WARNING: MAGNETS_OUTERINT nodeset not found." << std::endl;
            }
            double magnets_mass_radius_scattered = magnets_mass_radius_outerint / magnet_nodes->second.size();
            if (magnet_nodes != nset_map.end())
            {
                for (auto node_id_it = magnet_nodes->second.begin(); node_id_it != magnet_nodes->second.end(); ++node_id_it)
                {
                    auto node = inserted_nodes.at(*node_id_it);
                    auto centrifugal_force_vector = ChVector<>(node->GetPos().x(), node->GetPos().y(), 0.0);
                    centrifugal_force_vector.Normalize();
                    auto old_force = node->GetForce();
                    node->SetForce(old_force + magnets_mass_radius_scattered * omega*omega*centrifugal_force_vector);
                }
            }
        }

        {
            auto magnet_nodes = nset_map.find("MAGNETS_MEDIAN");
            if (magnet_nodes == nset_map.end())
            {
                std::cout << "WARNING: MAGNETS_MEDIAN nodeset not found." << std::endl;
            }
            double magnets_mass_radius_scattered = magnets_mass_radius_median / magnet_nodes->second.size();
            if (magnet_nodes != nset_map.end())
            {
                for (auto node_id_it = magnet_nodes->second.begin(); node_id_it != magnet_nodes->second.end(); ++node_id_it)
                {
                    auto node = inserted_nodes.at(*node_id_it);
                    auto centrifugal_force_vector = ChVector<>(node->GetPos().x(), node->GetPos().y(), 0.0);
                    centrifugal_force_vector.Normalize();
                    auto old_force = node->GetForce();
                    node->SetForce(old_force + magnets_mass_radius_scattered * omega*omega*centrifugal_force_vector);
                }
            }
        }

        ////////////// Apply forces given by EM pressure //////////////
        CSVwriter forces(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_EMforces.csv", LOG_OUTPUT);
        for (auto it = external_nodes.begin(); it != external_nodes.end(); ++it) {
            double angle = atan2((*it)->GetPos().y(), (*it)->GetPos().x());

            // assure angle between [0, 2*pi)
            angle = fmod(angle, CH_C_2PI);
            if (angle < 0)
                angle += CH_C_2PI;

            assert(angle <= CH_C_2PI);
            assert(angle >= 0);

#ifdef EQUAL_ELEMENT_SPACING
            double index = angle / CH_C_2PI * (sigma_n.size() - 1);
            auto index_int = static_cast<size_t>(floor(index));
            // int next_index = index_int > sigma_glob_set.size() ? index_int - sigma_glob_set.size() : index_int;
            int next_index = std::min(index_int + 1, sigma_n.size() - 1);

            ChVector<> sigma_loc;
            sigma_loc[0] = sigma_n[index_int];
            sigma_loc[1] = sigma_t[index_int];
            sigma_loc[0] += (index - index_int) * (sigma_n[next_index] - sigma_n[index_int]);
            sigma_loc[1] += (index - index_int) * (sigma_t[next_index] - sigma_t[index_int]);
            sigma_loc[2] = 0.0;

            ChVector<> sigma_glob;
            sigma_glob[0] = sigma_loc[0] * cos(angle) - sigma_loc[1] * sin(angle);
            sigma_glob[1] = sigma_loc[0] * sin(angle) + sigma_loc[1] * cos(angle);
            sigma_glob[2] = 0;

            ChVector<> forces_glob = sigma_glob;
            forces_glob.Scale(CH_C_2PI * rotor_external_radius * element_thickness / (external_nodes.size()));
            forces.AppendRow(angle, sigma_loc[0], sigma_loc[1], sigma_loc[2], forces_glob[0], forces_glob[1],
                forces_glob[2]);

#else
            auto iter_prev = angles.lower_bound(angle);
            auto iter_next = iter_prev;
            double halfangle_previous = 0.5 * ((iter_prev != angles.begin() ? *--iter_prev : *(--angles.end())) + angle);
            double halfangle_next = 0.5 * (((iter_next != angles.end()) && ++iter_next != angles.end() ? *iter_next : *(angles.begin())) + angle);
            ChVector<> sigma_glob_previous, sigma_glob_next, sigma_glob_center;

            getSigmaGlob(sigma_glob_set, sigma_glob_previous, halfangle_previous);
            getSigmaGlob(sigma_glob_set, sigma_glob_next, halfangle_next);
            getSigmaGlob(sigma_glob_set, sigma_glob_center, angle);
            double archLength_previous, archLength_next;
            getArcLength(angles, archLength_previous, archLength_next, angle, rotor_external_radius);

            auto EMforces_glob = 0.5 * element_thickness *(archLength_previous * sigma_glob_previous + archLength_next * sigma_glob_next);
            
            forces.AppendRow(angle, EMforces_glob[0], EMforces_glob[1], EMforces_glob[2]);

#endif
            auto old_force = (*it)->GetForce();
            (*it)->SetForce(old_force + EMforces_glob);
        }

        ////////////// Apply centrifugal forces //////////////
        // must be done after SetupInitial: volume is evaluated only at that time;
        // we could check if the omega is changed from the previous run and avoid re-evaluation
        CSVwriter b(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_centrifugal_forces.csv", LOG_OUTPUT);
        // std::map<std::shared_ptr<ChNodeFEAxyz>, ChVector<>> centrifugal_forces_map;
        for (auto el_it = my_mesh->GetElements().begin(); el_it != my_mesh->GetElements().end(); ++el_it) {
            auto el = std::dynamic_pointer_cast<ChElementHexa_8>(*el_it);
            ChVector<> mean_point;
            for (auto node_sel = 0; node_sel < 8; ++node_sel) {
                mean_point += std::dynamic_pointer_cast<ChNodeFEAxyz>(el->GetNodeN(node_sel))->GetX0();
            }
            mean_point *= 1.0 / 8.0;
            auto dist_from_center = sqrt(mean_point.x() * mean_point.x() + mean_point.y() * mean_point.y());
            auto centrifugal_force_modulus = el->GetVolume() * element_material->Get_density() * omega * omega * dist_from_center;
            auto centrifugal_force_vector = ChVector<>(mean_point.x(), mean_point.y(), 0.0);
            centrifugal_force_vector.Normalize();
            centrifugal_force_vector *= centrifugal_force_modulus / 8.0;
            for (auto node_sel = 0; node_sel < 8; ++node_sel) {
                auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(el->GetNodeN(node_sel));
                // store the centrifugal force and then apply it
                // centrifugal_forces_map[node] = centrifugal_force_vector;
                auto old_force = node->GetForce();
                node->SetForce(old_force + centrifugal_force_vector);
            }
            b.AppendRow(mean_point.x(), mean_point.y(), centrifugal_force_vector.x(), centrifugal_force_vector.y());
        }

        tim.stop();
        GetLog() << "Forces application time: " << tim() << "\n";

        ////////////// Run simulation //////////////
        tim.start();
        my_system.Setup();
        my_system.Update();
        my_system.DoStaticLinear();
        tim.stop();
        GetLog() << "Simulation time: " << tim() << "\n";


        ////////////// Export stress to file //////////////
        tim.reset();
        tim.start();

        // element stress export: all elements
        CSVwriter stress_file(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_stress_allelements.csv");
        for (auto el_it = my_mesh->GetElements().begin(); el_it != my_mesh->GetElements().end(); ++el_it)
        {
            auto el = std::dynamic_pointer_cast<ChElementHexa_8>(*el_it);
            auto stress = el->GetStress(0, 0, 0);
            //stress_file.AppendRow(inserted_elements_ptr_to_ID.at(el), stress.GetEquivalentVonMises(), stress.XX(), stress.YY(), stress.ZZ(), stress.XY(), stress.YZ(), stress.XZ());
            stress_file.AppendRow(inserted_elements_ptr_to_ID.at(el), stress.GetEquivalentVonMises());
        }
 

        // element stress export: bridge elements (subset of all elements)
        CSVwriter stressbridge_file(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_stress_bridges.csv");
        for (auto el_it = elset_map["BRIDGES"].begin(); el_it != elset_map["BRIDGES"].end(); ++el_it)
        {
            auto el_found = inserted_elements_ID_to_ptr.find(*el_it);
            if (el_found != inserted_elements_ID_to_ptr.end())
            {
                auto el = std::dynamic_pointer_cast<ChElementHexa_8>(el_found->second);
                auto stress = el->GetStress(0, 0, 0);
                stressbridge_file.AppendRow(*el_it, stress.GetEquivalentVonMises(), stress.XX(), stress.YY(), stress.ZZ(), stress.XY(), stress.YZ(), stress.XZ());
            }
        }


        // element stress export: bridge elements (subset of all elements)
        CSVwriter stressfatigue_file(resultspath + motor_prefix + "_" + drivingcyle_prefix + "_" + std::to_string(test_sel) + "_stress_forfatigue.csv");
        for (auto el_it = elset_map["FORFATIGUE"].begin(); el_it != elset_map["FORFATIGUE"].end(); ++el_it)
        {
            auto el_found = inserted_elements_ID_to_ptr.find(*el_it);
            if (el_found != inserted_elements_ID_to_ptr.end())
            {
                auto el = std::dynamic_pointer_cast<ChElementHexa_8>(el_found->second);
                auto stress = el->GetStress(0, 0, 0);
                stressfatigue_file.AppendRow(*el_it, stress.GetEquivalentVonMises(), stress.XX(), stress.YY(), stress.ZZ(), stress.XY(), stress.YZ(), stress.XZ());
            }
        }

        tim.stop();
        GetLog() << "Export time: " << tim() << "\n";
        GetLog() << "Done\n\n";



#ifdef USE_IRRLICHT
        ////////////// Rendering //////////////
        while (application.GetDevice()->run()) {
            application.BeginScene();

            application.DrawAll();

            if (!application.GetPaused()) {

            }

            application.EndScene();
        }
#endif

    }

    GetLog() << "All requested tests have been consumed.\n";

    return 0;
}