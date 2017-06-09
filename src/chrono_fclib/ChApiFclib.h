// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni
// =============================================================================

#ifndef CHAPIFCLIB_H
#define CHAPIFCLIB_H

#include "chrono/core/ChPlatform.h"

// Chrono::Engine version
//
// This is an integer, as 0xaabbccdd where
// for example version 1.2.0 is 0x00010200

#define CH_VERSION_FCLIB_MODULE 0x00000100

// When compiling this library, remember to define CH_API_COMPILE_FCLIB
// (so that the symbols with 'ChApiFclib' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE_FCLIB)
#define ChApiFclib ChApiEXPORT
#else
#define ChApiFclib ChApiIMPORT
#endif

/**
    @defgroup fclib_module FCLIB module
    @brief Module for FCLIB library

    For additional information, see:
    - the [installation guide](@ref module_fclib_installation)
*/

#endif
