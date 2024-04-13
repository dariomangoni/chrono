// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2023 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// External project template for loading an FMU (v2.0) as Cosimulation
// Targeting the "FmuComponentChrono" FMU generated by the same project.
// =============================================================================

#include <iostream>
#include <cstddef>

#include "chrono_fmi/fmi2/ChFmuToolsImport.h"

std::string unzipped_fmu_folder = FMU_UNPACK_DIRECTORY;
// std::string unzipped_fmu_folder = FMU_MAIN_DIRECTORY; // for debug
int main(int argc, char* argv[]) {
    FmuUnit my_fmu;

    try {
        // my_fmu.LoadUnzipped(unzipped_fmu_folder);
        my_fmu.Load(FMU_FILENAME, FMU_UNPACK_DIRECTORY);  // make sure the user has appropriate privileges to
                                                          // remove/create FMU_UNPACK_DIRECTORY
        // my_fmu.Load(FMU_FILENAME); // will go in TEMP/_fmu_temp

    } catch (std::exception& my_exception) {
        std::cout << "ERROR loading FMU: " << my_exception.what() << "\n";
    }

    std::cout << "FMU version:  " << my_fmu._fmi2GetVersion() << "\n";
    std::cout << "FMU platform: " << my_fmu._fmi2GetTypesPlatform() << "\n";

    my_fmu.Instantiate("FmuComponent", false, true);

    std::vector<std::string> categoriesVector = {"logAll"};
    my_fmu.SetDebugLogging(fmi2True, categoriesVector);

    // alternatively, with native interface:
    // std::vector<const char*> categoriesArray;
    // for (const auto& category : categoriesVector) {
    //    categoriesArray.push_back(category.c_str());
    //}

    // my_fmu._fmi2SetDebugLogging(my_fmu.component, fmi2True, categoriesVector.size(), categoriesArray.data());

    double start_time = 0;
    double stop_time = 2;
    my_fmu._fmi2SetupExperiment(my_fmu.component,
                                fmi2False,  // tolerance defined
                                0.0,        // tolerance
                                start_time,
                                fmi2False,  // use stop time
                                stop_time);

    my_fmu.EnterInitializationMode();
    // alternatively, with native interface: my_fmu._fmi2EnterInitializationMode(my_fmu.component);

    my_fmu.ExitInitializationMode();
    // alternatively, with native interface: my_fmu._fmi2ExitInitializationMode(my_fmu.component);

    // test a simulation loop:
    double time = 0;
    double dt = 0.001;

    for (unsigned int i = 0; i < 10000; ++i) {
        fmi2Status readStatus;
        my_fmu._fmi2DoStep(my_fmu.component, time, dt, fmi2True);

        double x = 0;
        readStatus = my_fmu.GetVariable("x", x, FmuVariable::Type::Real);
        std::cout << "x: " << x << std::endl;

        double theta = 0;
        readStatus = my_fmu.GetVariable("theta", theta, FmuVariable::Type::Real);
        std::cout << "theta: " << theta << std::endl;

        time += dt;
    }

    return 0;
}