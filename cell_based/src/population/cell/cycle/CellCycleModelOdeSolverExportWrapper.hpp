/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
