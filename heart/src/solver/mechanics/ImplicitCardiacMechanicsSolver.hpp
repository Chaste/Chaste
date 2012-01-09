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

private:
    /** The stretch in the fibre direction at the last timestep, in order
     *  to compute the stretch rate.
     *  Note the indexing: the i-th entry corresponds to the i-th global
     *  quad point, when looping over elements and then quad points
     */
    std::vector<double> mStretchesLastTimeStep;

    /** The current stretch ratio (in the fibre direction). Note the indexing:
     *  the i-th entry corresponds to the i-th global quad point, when looping
     *  over elements and then quad points
     */

    /** This solver is an implicit solver (overloaded pure method) */
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
    ImplicitCardiacMechanicsSolver(ContractionModelName contractionModel,
                                   QuadraticMesh<DIM>& rQuadMesh,
                                   SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                   std::string outputDirectory);

    /**
     *  Destructor
     */
    virtual ~ImplicitCardiacMechanicsSolver();


    /**
     *  Get lambda (the stretch ratio).
     *  NOTE: the i-th entry of this vector is assumed to be the i-th quad point
     *  obtained by looping over cells in the obvious way and then looping over
     *  quad points. These quad points, in the same order, can be obtained by
     *  using the QuadraturePointsGroup class.
     */
    std::vector<double>& rGetFibreStretches();

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
