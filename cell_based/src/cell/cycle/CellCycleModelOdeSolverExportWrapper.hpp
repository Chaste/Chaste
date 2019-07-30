/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CELLCYCLEMODELODESOLVEREXPORTWRAPPER_HPP_
#define CELLCYCLEMODELODESOLVEREXPORTWRAPPER_HPP_

#include "CellCycleModelOdeSolver.hpp"

// Possible ODE solvers.
// We might be able to just do "class CvodeAdaptor;" etc rather than #include here,
// but it probably wouldn't make much difference to build speed.
#include "CvodeAdaptor.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

#endif /*CELLCYCLEMODELODESOLVEREXPORTWRAPPERHPP_*/

// The following comes after the include guard, because it will expand to something
// different in .hpp and .cpp contexts

#ifdef EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER
// Avoid re-definition error in .cpp
#undef EXPORT_CCM_INTERNAL
#undef EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER
#endif // EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER

#define EXPORT_CCM_INTERNAL(CCM_CLASS, ODE_SOLVER) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CCM_CLASS, ODE_SOLVER)

#ifdef CHASTE_CVODE
#ifdef EXPORT_CCM_CVODE
// Avoid re-definition error in .cpp
#undef EXPORT_CCM_CVODE
#endif // EXPORT_CCM_CVODE
#define EXPORT_CCM_CVODE(CCM_CLASS) EXPORT_CCM_INTERNAL(CCM_CLASS, CvodeAdaptor)
#else
#define EXPORT_CCM_CVODE(CCM_CLASS)
#endif // CHASTE_CVODE

#define EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CCM_CLASS)         \
    EXPORT_CCM_CVODE(CCM_CLASS)                               \
    EXPORT_CCM_INTERNAL(CCM_CLASS, BackwardEulerIvpOdeSolver) \
    EXPORT_CCM_INTERNAL(CCM_CLASS, EulerIvpOdeSolver)         \
    EXPORT_CCM_INTERNAL(CCM_CLASS, HeunIvpOdeSolver)          \
    EXPORT_CCM_INTERNAL(CCM_CLASS, RungeKutta2IvpOdeSolver)   \
    EXPORT_CCM_INTERNAL(CCM_CLASS, RungeKutta4IvpOdeSolver)
