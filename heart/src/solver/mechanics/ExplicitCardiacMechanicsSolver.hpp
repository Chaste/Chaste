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

#ifndef EXPLICITCARDIACMECHANICSSOLVER_HPP_
#define EXPLICITCARDIACMECHANICSSOLVER_HPP_

#include "AbstractCardiacMechanicsSolver.hpp"
#include "AbstractCardiacMechanicsSolverInterface.hpp"
#include "AbstractContractionModel.hpp"


/**
 *  Explicit cardiac mechanics solver for solving electromechanic problems where the
 *  contraction model is not stretch-rate-dependent (for those the implicit solver is
 *  needed).
 *
 *  The first template parameter should be either IncompressibleNonlinearElasticitySolver
 *  or CompressibleNonlinearElasticityAssembler; this will be the class that this class
 *  ultimately inherits from.
 *
 *  The general explicit solution procedure is to do, each timestep:
 *  (0) [solve the electrics and interpolate Ca and voltage onto quad points
 *  (i) pass Ca and voltage to the contraction models
 *  (ii) pass the fibre stretch to the contraction models in case this is needed.
 *  (iii) integrate the contraction models in order to get the active tension
 *  (iv) solve for the deformation using this active tension.
 */
template<class ELASTICITY_SOLVER,unsigned DIM>
class ExplicitCardiacMechanicsSolver : public AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>
{
friend class TestExplicitCardiacMechanicsSolver;

private:
    /** @return true if this solver is an explicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return false;
    }

    /**
     *  Get the active tension and other info at the given quadrature point. This is an explicit
     *  solver so just sets the active tension, it doesn't set the derivatives. It stores the
     *  stretch for the next timestep.
     *
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex Quadrature point integrand currently being evaluated at in AssembleOnElement.
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension.
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch. Unset in this explicit solver.
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate. Unset in this explicit solver.
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
    ExplicitCardiacMechanicsSolver(QuadraticMesh<DIM>& rQuadMesh,
                                   ElectroMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                   std::string outputDirectory);


    /**
     *  Destructor
     */
    virtual ~ExplicitCardiacMechanicsSolver();


    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated explicitly
     *  over the time range using the ODE timestep provided then the active tension
     *  used to solve for the deformation
     *
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    void Solve(double time, double nextTime, double odeTimestep);
};

#endif /*EXPLICITCARDIACMECHANICSSOLVER_HPP_*/
