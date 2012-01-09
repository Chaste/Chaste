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
    /** This solver is an explicit solver (overloaded pure method) */
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

    /**
     * Initialise contraction models for each quadrature point
     * @param contractionModel The name of the contraction model (from the enumeration ContractionModel
     * defined in AbstractContactionModel)
     */
    void InitialiseContractionModels(ContractionModelName contractionModel);

public:
    /**
     * Constructor
     *
     * @param contractionModel The contraction model.
     * @param rQuadMesh A reference to the mesh.
     * @param rProblemDefinition Object defining body force and boundary conditions
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     */
    ExplicitCardiacMechanicsSolver(ContractionModelName contractionModel,
                                   QuadraticMesh<DIM>& rQuadMesh,
                                   SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
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
