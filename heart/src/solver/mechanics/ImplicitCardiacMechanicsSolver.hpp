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


#ifndef IMPLICITCARDIACMECHANICSSOLVER_HPP_
#define IMPLICITCARDIACMECHANICSSOLVER_HPP_

#include "AbstractCardiacMechanicsSolver.hpp"
#include "AbstractCardiacMechanicsSolverInterface.hpp"
#include "LogFile.hpp"
#include <cfloat>


/**
 *  Implicit Cardiac Mechanics Solver
 *
 *  The first template parameter should be either IncompressibleNonlinearElasticitySolver
 *  or CompressibleNonlinearElasticityAssembler; this will be the class that this class
 *  ultimately inherits from.
 *
 *  Solves cardiac mechanics implicitly (together with the contraction
 *  models for determining the active tension), taking in the intracellular
 *  calcium concentration. See CardiacElectroMechanicsProblem documentation
 *  for more detail.
 */
template<class ELASTICITY_SOLVER, unsigned DIM>
class ImplicitCardiacMechanicsSolver : public AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>
{
friend class TestImplicitCardiacMechanicsSolver;
friend class TestExplicitCardiacMechanicsSolver;//for comparison test
private:
    /** @return true if this solver is an implicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return true;
    }

    /**
     *  A method called by AbstractCardiacMechanicsSolver::AssembleOnElement(), providing
     *  the active tension (and other info) at a particular quadrature point. This version uses C to
     *  determine the current stretch and stretch rate, and integrates the contraction model ODEs to determine
     *  the active tension, and derivatives of the active tension with respect to stretch and
     *  stretch rate.
     *
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex Quadrature point integrand currently being evaluated at in AssembleOnElement.
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate
     */
    void GetActiveTensionAndTensionDerivs(double currentFibreStretch,
                                          unsigned currentQuadPointGlobalIndex,
                                          bool assembleJacobian,
                                          double& rActiveTension,
                                          double& rDerivActiveTensionWrtLambda,
                                          double& rDerivActiveTensionWrtDLambdaDt);

public:
    /**
     * Constructor
     *
     * @param rQuadMesh A reference to the mesh.
     * @param rProblemDefinition Object defining body force and boundary conditions
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     */
    ImplicitCardiacMechanicsSolver(QuadraticMesh<DIM>& rQuadMesh,
                                   ElectroMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                   std::string outputDirectory);

    /**
     *  Destructor
     */
    virtual ~ImplicitCardiacMechanicsSolver();

    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated implicitly
     *  over the time range using the ODE timestep provided, as part of the solve,
     *  and updated at the end once the solution has been found, as is lambda.
     *
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    void Solve(double time, double nextTime, double odeTimestep);
};

#endif /*IMPLICITCARDIACMECHANICSSOLVER_HPP_*/
