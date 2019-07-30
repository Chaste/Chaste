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

#ifndef ABSTRACTNONLINEARELASTICITYSOLVER_HPP_
#define ABSTRACTNONLINEARELASTICITYSOLVER_HPP_

#include <vector>
#include <cmath>
#include "AbstractContinuumMechanicsSolver.hpp"
#include "LinearSystem.hpp"
#include "LogFile.hpp"
#include "MechanicsEventHandler.hpp"
#include "ReplicatableVector.hpp"
#include "FourthOrderTensor.hpp"
#include "CmguiDeformedSolutionsWriter.hpp"
#include "AbstractMaterialLaw.hpp"
#include "QuadraticBasisFunction.hpp"
#include "SolidMechanicsProblemDefinition.hpp"
#include "Timer.hpp"
#include "AbstractPerElementWriter.hpp"
#include "petscsnes.h"

//#define MECH_USE_HYPRE    // uses HYPRE to solve linear systems, requires PETSc to be installed with HYPRE

/**
 *  Three options for which type of strain to write
 *  F = dx/dX, C = F^T F, E = 1/2 (C-I)
 */
typedef enum StrainType_
{
    DEFORMATION_GRADIENT_F = 0,
    DEFORMATION_TENSOR_C,
    LAGRANGE_STRAIN_E
} StrainType;

// Bizarrely PETSc 2.2 has this, but doesn't put it in the petscksp.h header...
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
extern PetscErrorCode KSPInitialResidual(KSP,Vec,Vec,Vec,Vec,Vec);
#endif

//////////////////////////////////////////////////////////////
//  Globals functions used by the SNES solver
//////////////////////////////////////////////////////////////

/**
 *  Global function that will be called by the SNES solver
 *
 *  @param snes snes solver
 *  @param currentGuess current guess
 *  @param residualVector residual vector to be computed
 *  @param pContext this pointer can be converted to a ptr to the original caller, ie the
 *  AbstractNonlinearElasticitySolver class
 */
template<unsigned DIM>
PetscErrorCode AbstractNonlinearElasticitySolver_ComputeResidual(SNES snes,
                                                                 Vec currentGuess,
                                                                 Vec residualVector,
                                                                 void* pContext);

#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5))
    /**
     *  Global function that will be called by the SNES solver
     *
     *  Note - this changed in PETSC 3.5.
     *
     *  @param snes snes solver
     *  @param currentGuess current guess
     *  @param globalJacobian jacobian matrix to be computed
     *  @param preconditioner preconditioner matrix to be computed
     *  @param pContext this pointer can be converted to a ptr to the original caller, ie the
     *  AbstractNonlinearElasticitySolver class
     */
    template<unsigned DIM>
    PetscErrorCode AbstractNonlinearElasticitySolver_ComputeJacobian(SNES snes,
                                                                     Vec currentGuess,
                                                                     Mat globalJacobian,
                                                                     Mat preconditioner,
                                                                     void* pContext);
#else
/**
 *  Global function that will be called by the SNES solver
 *
 *  @param snes snes solver
 *  @param currentGuess current guess
 *  @param pGlobalJacobian jacobian matrix to be computed
 *  @param pPreconditioner preconditioner matrix to be computed
 *  @param pMatStructure  The PETSc matrix structure description.
 *  @param pContext this pointer can be converted to a ptr to the original caller, ie the
 *  AbstractNonlinearElasticitySolver class
 */
template<unsigned DIM>
PetscErrorCode AbstractNonlinearElasticitySolver_ComputeJacobian(SNES snes,
                                                                 Vec currentGuess,
                                                                 Mat* pGlobalJacobian,
                                                                 Mat* pPreconditioner,
                                                                 MatStructure* pMatStructure,
                                                                 void* pContext);
#endif

template <unsigned DIM>
class  AbstractNonlinearElasticitySolver; //Forward declaration

/**
 * Helper class for dumping the stresses to file.
 *
 * Currently located here so that it's easy to feed a pointer to the main
 * class AbstractNonlinearElasticitySolver
 */
template<unsigned DIM>
class StressPerElementWriter : public AbstractPerElementWriter<DIM, DIM, DIM*DIM>
{
private:
    AbstractNonlinearElasticitySolver<DIM>* mpSolver; /**< Pointer to the parent class, set in constructor*/
public:

    /** Constructor
     * @param pMesh  A pointer to the mesh whose elements we want to calculate data on.
     * @param pSolver  A pointer to the parent class, used to access data
     */
    StressPerElementWriter(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                           AbstractNonlinearElasticitySolver<DIM>* pSolver)
     : AbstractPerElementWriter<DIM, DIM, DIM*DIM>(pMesh),
       mpSolver(pSolver)
    {
    }

    /**
     * How to associate an element with the stress data
     *
     * @param pElement  a locally-owned element for which to calculate or lookup some data
     * @param localIndex the index (in mElements) of pElement
     * @param rData  the double-precision data to write to file (output from the method)
     */
    void Visit(Element<DIM, DIM>* pElement, unsigned localIndex, c_vector<double, DIM*DIM>& rData)
    {
        ///\todo #2223 Note that the localIndex is intended for DistributedQuadraticMesh, but we might need a map instead.
        /// For now, make sure that it's not a distributed mesh i.e. GetNumLocalElements == GetNumElements
        /// and the indexing is always global
        assert(localIndex == pElement->GetIndex()); //Will fail when we move to DistributedQuadraticMesh
        //Flatten the matrix
        c_matrix<double, DIM, DIM> data = mpSolver->GetAverageStressPerElement(localIndex);
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                rData[i*DIM+j] = data(i,j);
            }
        }
    }
};

/**
 * Abstract nonlinear elasticity solver. IncompressibleNonlinearElasticityAssembler and
 * CompressibleNonlinearElasticityAssembler inherit from this class.
 *
 * The class is both a solver AND a assembler: the AssembleOnElement(), AssembleSystem() methods
 * are hardcoded into this class. In principle something like AbstractContinuumMechanicsAssembler
 * could have been used as a member variable [that class is used for assembling fluids matrices etc]
 * but nonlinear elasticity is too complex for the that class to be used, as things like stress and
 * stress-derivative need to be computed at the AssembleOnElement level, things like
 * pressure-on-deformed-surface are deformation dependent boundary conditions, etc. *
 */
template<unsigned DIM>
class AbstractNonlinearElasticitySolver : public AbstractContinuumMechanicsSolver<DIM>
{
    friend class StressRecoveror<DIM>;
    friend class VtkNonlinearElasticitySolutionWriter<DIM>;


protected:

    /** Number of vertices per element. */
    static const size_t NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const size_t NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2;      // assuming quadratic

    /** Number of nodes per boundary element. */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2; // assuming quadratic

    /** Boundary stencil size. Note this is just the number of spatial unknowns on the boundary
     *  element, because the boundary integral terms (in either compressible or incompressible
     *  formulations) (i) do not involve pressure and (ii) do not appear in the pressure
     *  equations (the constraint equations). */
    static const size_t BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT;

    /**
     * Maximum absolute tolerance for Newton solve. The Newton solver uses the absolute tolerance
     * corresponding to the specified relative tolerance, but has a max and min allowable absolute
     * tolerance. (ie: if max_abs = 1e-7, min_abs = 1e-10, rel=1e-4: then if the norm of the
     * initial_residual (=a) is 1e-2, it will solve with tolerance 1e-7; if a=1e-5, it will solve
     * with tolerance 1e-9; a=1e-9, it will solve with tolerance 1e-10.
     */
    static double MAX_NEWTON_ABS_TOL;

    /** Minimum absolute tolerance for Newton solve. See documentation for MAX_NEWTON_ABS_TOL. */
    static double MIN_NEWTON_ABS_TOL;

    /** Relative tolerance for Newton solve. See documentation for MAX_NEWTON_ABS_TOL. */
    static double NEWTON_REL_TOL;

    /**
     *  This class contains all the information about the problem (except the material law):
     *  body force, surface tractions, fixed nodes, density
     */
    SolidMechanicsProblemDefinition<DIM>& mrProblemDefinition;

    /**
     * Reference to the matrix 'mLhsSystemMatrix' in the parent class, named as the jacobian,
     * just to make code clearer.
     */
    Mat& mrJacobianMatrix;

    /**
     *  Matrix to be passed to material law, in case the material law is anisotropic and depends on
     *  a local coordinate system (eg cardiac problems). Defaults to the identity matrix.
     *  The cardiac mechanics solvers override SetupChangeOfBasisMatrix() which re-sets this
     *  for every element and quad point.
     */
    c_matrix<double,DIM,DIM> mChangeOfBasisMatrix;

    /**
     * Absolute tolerance for linear systems. Can be set by calling
     * SetKspAbsoluteTolerances(), but default to -1, in which case
     * a relative tolerance is used.
     */
    double mKspAbsoluteTol;

    /**
     * By default only the initial and final solutions are written. However, we may
     * want to write the solutions after every Newton iteration, in which case the
     * following should be set to true.
     */
    bool mWriteOutputEachNewtonIteration;

    /**
     * Storage space for a 4th order tensor used in assembling the Jacobian (to avoid repeated memory allocation).
     */
    FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

    /** Number of Newton iterations taken in last solve. */
    unsigned mNumNewtonIterations;

    /**
     * This solver is for static problems, however the body force or surface tractions
     * could be a function of time. The user should call SetCurrentTime() if this is
     * the case.
     */
    double mCurrentTime;


    /** Whether the boundary elements of the mesh have been checked for whether
     *  the ordering if such that the normals are outward-facing
     */
    bool mCheckedOutwardNormals;

    /**
     * If SetUsingSnesSolver() is called on the problem definition class, or the command line argument
     *  "-mech_use_snes" is given the Petsc SNES solver (nonlinear solver) will be used.
     */
    bool mUseSnesSolver;

    /**
     *  The last damping value used in the current nonlinear (non-snes) solve. If near 1.0, this indicates
     *  that the current guess is near the solution. Initialised to 0.0 at the beginning of each nonlinear
     *  solve.
     */
    double mLastDampingValue;

    /**
     *  Whether this is the first Newton iteration or not
     */
    bool mFirstStep;

    /**
     *  Whether to take a full first Newton step or not - see documentation for
     *  SetTakeFullFirstNewtonStep()
     */
    bool mTakeFullFirstNewtonStep;

    /**
     *  Get Petsc to do a (almost) direct solve
     */
    bool mPetscDirectSolve;

    /**
     * Whether to call AddActiveStressAndStressDerivative() when computing stresses or not.
     *
     * Subclasses, such as the cardiac mechanics solvers, may implement the above method to add
     * active contributions to the stress. However, sometimes we might want to switch this off.
     * Defaults to true.
     */
    bool mIncludeActiveTension;

    /**
     * The user may request that the stress for each element (averaged over quadrature point stresses)
     * are saved during the Solve; this bool states this (defaults to false).
     */
    bool mSetComputeAverageStressPerElement;

    /**
     * If the stresses for each element (averaged over quadrature point stresses)
     * are to be stored, they are stored in this variable.
     * Note to save memory we just don't store the lower half of the stress,
     * as the stress is symmetric, hence this is a vector of 6 (in 3d) variables
     * rather than a 3d matrix.
     */
    std::vector<c_vector<double,DIM*(DIM+1)/2> > mAverageStressesPerElement;

    /**
     *  Add the given stress tensor to the store of average stresses.
     *  mSetComputeAverageStressPerElement must be true
     *
     *  @param rT 2nd PK stress (matrix is assumed symmetric)
     *  @param elementIndex element index
     */
    void AddStressToAverageStressPerElement(c_matrix<double,DIM,DIM>& rT, unsigned elementIndex);

    /**
     * Set the KSP type (CG, GMRES, etc) and the preconditioner type (ILU, ICC etc). Depends on
     * incompressible or not, and other factors.
     *
     * ///\todo #2057 Make better choices here...
     *
     * @param solver KSP solver object (Petsc object)
     */
    virtual void SetKspSolverAndPcType(KSP solver);


    /**
     * Assemble the residual vector and/or Jacobian matrix (using the current solution stored
     * in mCurrentSolution, output going to mResidualVector and/or mrJacobianMatrix).
     *
     * Must be overridden in concrete derived classes.
     *
     * @param assembleResidual     A bool stating whether to assemble the residual vector.
     * @param assembleLinearSystem A bool stating whether to assemble the Jacobian matrix and the RHS
     *  vector of the linear system (which is based on the residual but could be slightly different
     *  due to the way dirichlet boundary conditions are applied to the linear system - see comments in
     *  ApplyDirichletBoundaryConditions).
     */
    virtual void AssembleSystem(bool assembleResidual, bool assembleLinearSystem)=0;

    /**
     * To be called at the end of AssembleSystem. Calls (Petsc) assemble methods on the
     * Vecs and Mat, and calls ApplyDirichletBoundaryConditions.
     *
     * @param assembleResidual see documentation for AssembleSystem
     * @param assembleLinearSystem see documentation for AssembleSystem
     */
    virtual void FinishAssembleSystem(bool assembleResidual, bool assembleLinearSystem);

    ////////////////////////////////////////////
    //
    // Element level methods
    //
    ////////////////////////////////////////////

    /**
     * Compute the strain at the centroid at an element
     * @param strainType Which strain to compute, should be one of: DEFORMATION_GRADIENT_F, DEFORMATION_TENSOR_C, or LAGRANGE_STRAIN_E
     * @param rElement The element
     * @param rDeformationGradient Reference to a matrix, which will be filled in
     * by this method.
     */
    void GetElementCentroidStrain(StrainType strainType,
                                  Element<DIM,DIM>& rElement,
                                  c_matrix<double,DIM,DIM>& rDeformationGradient);

    /**
     * Add on the active component to the stress (and maybe also to the stress-derivative).
     * This is called after getting the passive stress and stress-derivative from
     * a material law.
     *
     * This method does nothing but can be overloaded by the other solvers, see eg cardiac
     * mechanics solvers.
     *
     * @param rC The Lagrangian deformation tensor (F^T F)
     * @param elementIndex Index of the current element
     * @param currentQuadPointGlobalIndex The index (assuming an outer loop over elements and an inner
     *   loop over quadrature points), of the current quadrature point.
     * @param rT The stress to be added to
     * @param rDTdE the stress derivative to be added to, assuming
     *   the final parameter is true
     * @param addToDTdE A boolean flag saying whether the stress derivative is
     *   required or not.
     */
    virtual void AddActiveStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                                    unsigned elementIndex,
                                                    unsigned currentQuadPointGlobalIndex,
                                                    c_matrix<double,DIM,DIM>& rT,
                                                    FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                    bool addToDTdE)
    {
        // does nothing - subclass needs to overload this if there are active stresses
    }

    /**
     * Should re-set the variable mChangeOfBasisMatrix if this is not to be the identity.
     * Here, does nothing, but over-ridden in cardiac mechanics solvers.
     *
     * @param elementIndex element global index
     * @param currentQuadPointGlobalIndex global index of the quad point
     */
    virtual void SetupChangeOfBasisMatrix(unsigned elementIndex, unsigned currentQuadPointGlobalIndex)
    {
    }

    /**
     * Compute the term from the surface integral of s*phi, where s is
     * a specified non-zero surface traction (ie Neumann boundary condition)
     * to be added to the residual vector.
     *
     * Calls AssembleOnBoundaryElementForPressureOnDeformedBc() if appropriate
     *
     * @param rBoundaryElement the boundary element to be integrated on
     * @param rAelem The element's contribution to the LHS matrix is returned. There is no
     *     need to zero this matrix before calling.
     * @param rBelem The element's contribution to the RHS vector is returned. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     * @param boundaryConditionIndex index of this boundary (in the vectors
     *     in the problem definition object, in which the boundary conditions are
     *     stored
     */
    void AssembleOnBoundaryElement(BoundaryElement<DIM-1, DIM>& rBoundaryElement,
                                   c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE>& rAelem,
                                   c_vector<double, BOUNDARY_STENCIL_SIZE>& rBelem,
                                   bool assembleResidual,
                                   bool assembleJacobian,
                                   unsigned boundaryConditionIndex);


    /**
     * When pressure-on-deformed-body boundary conditions are used:
     *
     * For some reason the use of the Jacobian corresponding to the pressure term doesn't help
     * (makes things worse!) if the current guess is not close enough to the solution, and can
     * lead to immediate divergence. It will lead to faster convergence once close enough to the
     * solution. This method contains the logic used to decide whether to include the jacobian
     * for this term or not.
     *
     * We only assemble the contribution to the matrix if the last damping value is
     * close enough to 1 (non-snes solver). This will current always return false if the snes solver
     * is being used
     * @return true if we are assembling pressure-on-deformed-body
     */
    bool ShouldAssembleMatrixTermForPressureOnDeformedBc();

    /**
     *  Alternative version of AssembleOnBoundaryElement which is used for problems where a normal pressure
     *  is applied to the deformed body. The traction then depends on the current deformation, specifically
     *  s = -det(F)*P*F^{-T}N.
     *
     *  To compute F we have to find the volume element containing this surface element and use this in the
     *  computation. See comments in implementation for more details.
     *
     * @param rBoundaryElement the boundary element to be integrated on
     * @param rAelem The element's contribution to the LHS matrix is returned. There is no
     *     need to zero this matrix before calling.
     * @param rBelem The element's contribution to the RHS vector is returned. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     * @param boundaryConditionIndex index of this boundary (in the vectors
     *     in the problem definition object, in which the boundary conditions are
     *     stored
     */
    void AssembleOnBoundaryElementForPressureOnDeformedBc(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                                          c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
                                                          c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
                                                          bool assembleResidual,
                                                          bool assembleJacobian,
                                                          unsigned boundaryConditionIndex);

    /////////////////////////////////////////////////////////////
    //
    //    These methods form the non-SNES nonlinear solver
    //
    /////////////////////////////////////////////////////////////

    /**
     * Set up the residual vector (using the current solution), and get its
     * scaled norm (Calculate |r|_2 / length(r), where r is residual vector).
     *
     * @param allowException Sometimes the current solution solution will be such
     *   that the residual vector cannot be computed, as (say) the material law
     *   will throw an exception as the strains are too large. If this bool is
     *   set to true, the exception will be caught, and DBL_MAX returned as the
     *   residual norm
     * @return residual norm
     */
    double ComputeResidualAndGetNorm(bool allowException);

    /** @return  |r|_2 / length(r), where r is the current residual vector. */
    double CalculateResidualNorm();

    /**
     * Simple helper function, computes Z = X + aY, where X and Z are std::vectors and
     * Y a ReplicatableVector.
     *
     * @param rX X
     * @param rY Y (replicatable vector)
     * @param a a
     * @param rZ Z the returned vector
     */
    void VectorSum(std::vector<double>& rX, ReplicatableVector& rY, double a, std::vector<double>& rZ);

    /**
     * Print to std::cout the residual norm for this s, ie ||f(x+su)|| where f is the residual vector,
     * x the current solution and u the update vector.
     *
     * @param s s
     * @param residNorm residual norm.
     */
    void PrintLineSearchResult(double s, double residNorm);

    /**
     * Take one newton step, by solving the linear system -Ju=f, (J the jacobian, f
     * the residual, u the update), and picking s such that a_new = a_old + su (a
     * the current solution) such |f(a)| is the smallest.
     *
     * @return The current norm of the residual after the newton step.
     */
    double TakeNewtonStep();

    /**
     * Using the update vector (of Newton's method), choose s such that ||f(x+su)|| is most decreased,
     * where f is the residual vector, x the current solution (mCurrentSolution) and u the update vector.
     * This checks s=1 first (most likely to be the current solution, then 0.9, 0.8.. until ||f|| starts
     * increasing.
     *
     * @param solution The solution of the linear solve in newton's method, ie the update vector u.
     * @return norm of residual
     */
    double UpdateSolutionUsingLineSearch(Vec solution);

    /**
     * This function may be overloaded by subclasses. It is called after each Newton
     * iteration.
     *
     * @param counter current newton iteration number
     * @param normResidual norm of the residual
     */
    virtual void PostNewtonStep(unsigned counter, double normResidual);

    /**
     *  Solve method which uses a nonlinear solver coded in this class (as opposed
     *  to the SNES solver. Private, user should call Solve()
     *  @param tol absolute solver used in nonlinear solve
     */
    void SolveNonSnes(double tol=-1.0);


    /////////////////////////////////////////////////////////////
    //
    //    These methods form the SNES nonlinear solver
    //
    /////////////////////////////////////////////////////////////
public: // need to be public as are called by global functions
    /**
     * Public method for computing the residual that will be called, effectively,
     * by the SNES solver
     * @param currentGuess Input, the current guess for the solution
     * @param residualVector Output, the residual vector given this guess.
     */
    void ComputeResidual(Vec currentGuess, Vec residualVector);

    /**
     * Public method for computing the jacobian that will be called, effectively,
     * by the SNES solver
     * @param currentGuess Input, the current guess for the solution
     * @param pJacobian Output, the jacobian matrix at this guess
     * @param pPreconditioner Output, the preconditioner matrix
     */
    void ComputeJacobian(Vec currentGuess, Mat* pJacobian, Mat* pPreconditioner);

private:
    /**
     * Alternative solve method which uses a Petsc SNES solver.
     * Private, user should call Solve()
     */
    void SolveSnes();

public:

    /**
     * Constructor.
     *
     * @param rQuadMesh  the quadratic mesh
     * @param rProblemDefinition an object defining in particular the body force and boundary conditions
     * @param outputDirectory output directory
     * @param compressibilityType Should be equal to COMPRESSIBLE or INCOMPRESSIBLE (see enumeration defined at top of file)
     *   (depending on which concrete class is inheriting from this) and is only used in computing mNumDofs and allocating
     *   matrix memory.
     */
    AbstractNonlinearElasticitySolver(AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
                                      SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                      std::string outputDirectory,
                                      CompressibilityType compressibilityType);

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearElasticitySolver();

    /**
     * Solve the problem.
     *
     * @param tol tolerance used in Newton solve (defaults to -1.0). Not used in SNES solves.
     */
    void Solve(double tol=-1.0);

    /**
     * Whether to call AddActiveStressAndStressDerivative() when computing stresses or not.
     *
     * Subclasses, such as the cardiac mechanics solvers, may implement the above method to add
     * active contributions to the stress. However, sometimes we might want to switch this off,
     * which is what this function is for - will generally be called with includeActiveTension=false
     *
     * @param includeActiveTension whether to include active tension
     */
    void SetIncludeActiveTension(bool includeActiveTension = true)
    {
        mIncludeActiveTension = includeActiveTension;
    }

    /**
     * @return number of Newton iterations taken in last solve.
     */
    unsigned GetNumNewtonIterations();


    /**
     * By default only the original and converged solutions are written. Call this
     * by get node positions written after every Newton step (mostly for
     * debugging).
     *
     * @param writeOutputEachNewtonIteration whether to write each iteration
     */
    void SetWriteOutputEachNewtonIteration(bool writeOutputEachNewtonIteration=true)
    {
        mWriteOutputEachNewtonIteration = writeOutputEachNewtonIteration;
    }

    /**
     * Set the absolute tolerance to be used when solving the linear system.
     * If this is not called a relative tolerance is used.
     *
     * @param kspAbsoluteTolerance the tolerance
     */
    void SetKspAbsoluteTolerance(double kspAbsoluteTolerance)
    {
        assert(kspAbsoluteTolerance > 0);
        mKspAbsoluteTol = kspAbsoluteTolerance;
    }

    /**
     *  The following odd behaviour has been observed: for some problems the solver
     *  will fail in the first Newton iteration, with the residual not decreasing
     *  in the direction of the Newton update, but: if you take a full Newton
     *  step anyway (increasing the residual-norm), the solver then converges
     *  perfectly. This method allows the user to choose this.
     *
     *  Does nothing if the SNES solver is used.
     *
     *  See ticket #2304
     *
     *  @param takeFullFirstStep Whether to take a full first Newton step or not.
     */
    void SetTakeFullFirstNewtonStep(bool takeFullFirstStep = true)
    {
        mTakeFullFirstNewtonStep = takeFullFirstStep;
    }

    /**
     *  Get Petsc to do a direct solve on the linear system (instead of using
     *  an iterative solve). This is equivalent to passing in command line
     *  arguments -ksp_type pre_only -pc_type lu through to Petsc, but in the incompressible
     *  case the preconditioner is set equal to the Jacobian with a mass matrix in the
     *  pressure-pressure block (to avoid zeros on the diagonal. Hence a few linear solve
     *  iterations are required for this case. Using a direct solve can lead to huge
     *  computation time savings if
     *  there is enough memory for it: the linear solve may be faster and
     *  nonlinear convergence likely to be much better, as the linear solve is exact.
     *
     *  @param usePetscDirectSolve Whether to use the Petsc direct solver or not
     */
    void SetUsePetscDirectSolve(bool usePetscDirectSolve = true)
    {
        mPetscDirectSolve = usePetscDirectSolve;
    }


    /**
     * This solver is for static problems, however the body force or surface tractions
     * could be a function of time. This method is for setting the time.
     *
     * @param time current time
     */
    void SetCurrentTime(double time)
    {
        mCurrentTime = time;
    }

    /**
     * Convert the output to Cmgui format (placed in a folder called cmgui in the output directory).
     * Writes the original mesh as solution_0.exnode and the (current) solution as solution_1.exnode.
     */
    void CreateCmguiOutput();


    /**
     * Write the strain for each element (evaluated at the centroids of each element). Which strain
     * to compute is determined by the first input parameter, and will be either F, C or E.
     * Each line of the output file corresponds to one element: the DIM*DIM matrix will be written
     * as one line, using the following ordering (assuming F is written).
     * F00 F01 F02 F10 F11 F12 F20 F21 F22.
     *
     * @param strainType Which strain to write, should be one of: DEFORMATION_GRADIENT_F, DEFORMATION_TENSOR_C, or LAGRANGE_STRAIN_E
     * @param fileName The file name stem
     * @param counterToAppend (Optional) number to append in the filename.
     *
     * The final file is [fileName]_[counterToAppend].strain
     */
    void WriteCurrentStrains(StrainType strainType, std::string fileName, int counterToAppend = -1);

    /**
     * The user may request that the stress for each element (averaged over quadrature point stresses)
     * are saved during the Solve(), by calling this.
     *
     * @param setComputeAverageStressPerElement whether to compute stresses (defaults to true)
     */
    void SetComputeAverageStressPerElementDuringSolve(bool setComputeAverageStressPerElement = true);

    /**
     * If SetComputeAverageStressPerElementDuringSolve() was called before the Solve(), then
     * this method can be used to print the average stresses to file)
     *
     * Each line of the output file corresponds to one element: the DIM*DIM matrix will be written
     * as one line, using the ordering:
     * T00 T01 T02 T10 T11 T12 T20 T21 T22.
     *
     * @param fileName The file name stem
     * @param counterToAppend (Optional) number to append in the filename.
     *
     * The final file is [fileName]_[counterToAppend].stress
     */
    void WriteCurrentAverageElementStresses(std::string fileName, int counterToAppend = -1);

    /**
     * @return the deformed position.
     * Note: return_value[i](j) = x_j for node i.
     */
    std::vector<c_vector<double,DIM> >& rGetSpatialSolution();

    /**
     * @return the deformed position. Note: return_value[i](j) = x_j for node i. Just
     * calls rGetSpatialSolution().
     */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();

    /**
     * If SetComputeAverageStressPerElementDuringSolve() was called before the Solve(), then
     * this method can be used to get the average stress for a particular
     * element.
     *
     * @param elementIndex elementIndex
     * @return stress tensor
     */
    c_matrix<double,DIM,DIM> GetAverageStressPerElement(unsigned elementIndex);
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation: first, the non-nonlinear-solve methods
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
AbstractNonlinearElasticitySolver<DIM>::AbstractNonlinearElasticitySolver(AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
                                                                          SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                          std::string outputDirectory,
                                                                          CompressibilityType compressibilityType)
    : AbstractContinuumMechanicsSolver<DIM>(rQuadMesh, rProblemDefinition, outputDirectory, compressibilityType),
      mrProblemDefinition(rProblemDefinition),
      mrJacobianMatrix(this->mSystemLhsMatrix),
      mKspAbsoluteTol(-1),
      mWriteOutputEachNewtonIteration(false),
      mNumNewtonIterations(0),
      mCurrentTime(0.0),
      mCheckedOutwardNormals(false),
      mLastDampingValue(0.0),
      mIncludeActiveTension(true),
      mSetComputeAverageStressPerElement(false)
{
    mUseSnesSolver = (mrProblemDefinition.GetSolveUsingSnes() ||
                      CommandLineArguments::Instance()->OptionExists("-mech_use_snes") );

    mChangeOfBasisMatrix = identity_matrix<double>(DIM,DIM);

    mTakeFullFirstNewtonStep = CommandLineArguments::Instance()->OptionExists("-mech_full_first_newton_step");
    mPetscDirectSolve = CommandLineArguments::Instance()->OptionExists("-mech_petsc_direct_solve");
}

template<unsigned DIM>
AbstractNonlinearElasticitySolver<DIM>::~AbstractNonlinearElasticitySolver()
{
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::FinishAssembleSystem(bool assembleResidual, bool assembleJacobian)
{
    PetscVecTools::Finalise(this->mResidualVector);

    if (assembleJacobian)
    {
        PetscMatTools::SwitchWriteMode(mrJacobianMatrix);
        PetscMatTools::SwitchWriteMode(this->mPreconditionMatrix);

        VecCopy(this->mResidualVector, this->mLinearSystemRhsVector);
    }

    // Apply Dirichlet boundary conditions
    if (assembleJacobian)
    {
        this->ApplyDirichletBoundaryConditions(NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING, this->mCompressibilityType==COMPRESSIBLE);
    }
    else if (assembleResidual)
    {
        this->ApplyDirichletBoundaryConditions(NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY, this->mCompressibilityType==COMPRESSIBLE);
    }

    if (assembleResidual)
    {
        PetscVecTools::Finalise(this->mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::Finalise(mrJacobianMatrix);
        PetscMatTools::Finalise(this->mPreconditionMatrix);
        PetscVecTools::Finalise(this->mLinearSystemRhsVector);
    }
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& AbstractNonlinearElasticitySolver<DIM>::rGetSpatialSolution()
{
    this->mSpatialSolution.clear();
    this->mSpatialSolution.resize(this->mrQuadMesh.GetNumNodes(), zero_vector<double>(DIM));
    for (unsigned i=0; i<this->mrQuadMesh.GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            this->mSpatialSolution[i](j) = this->mrQuadMesh.GetNode(i)->rGetLocation()[j] + this->mCurrentSolution[this->mProblemDimension*i+j];
        }
    }
    return this->mSpatialSolution;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& AbstractNonlinearElasticitySolver<DIM>::rGetDeformedPosition()
{
    return rGetSpatialSolution();
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::WriteCurrentStrains(StrainType strainType, std::string fileName, int counterToAppend)
{
    if (!this->mWriteOutput)
    {
        return;
    }

    std::stringstream file_name;
    file_name << fileName;
    if (counterToAppend >= 0)
    {
        file_name << "_" << counterToAppend;
    }
    file_name << ".strain";

    out_stream p_file = this->mpOutputFileHandler->OpenOutputFile(file_name.str());

    c_matrix<double,DIM,DIM> strain;

    for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = this->mrQuadMesh.GetElementIteratorBegin();
         iter != this->mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        GetElementCentroidStrain(strainType, *iter, strain);
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                *p_file << strain(i,j) << " ";
            }
        }
        *p_file << "\n";
    }
    p_file->close();
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::WriteCurrentAverageElementStresses(std::string fileName, int counterToAppend)
{
    if (!this->mWriteOutput)
    {
        return;
    }

    if (!mSetComputeAverageStressPerElement)
    {
        EXCEPTION("Call SetComputeAverageStressPerElementDuringSolve() before solve if calling WriteCurrentAverageElementStresses()");
    }

    std::stringstream file_name;
    file_name << fileName;
    if (counterToAppend >= 0)
    {
        file_name << "_" << counterToAppend;
    }
    file_name << ".stress";
    assert(mAverageStressesPerElement.size()==this->mrQuadMesh.GetNumElements());

    StressPerElementWriter<DIM> stress_writer(&(this->mrQuadMesh),this);
    stress_writer.WriteData(*(this->mpOutputFileHandler), file_name.str());
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::CreateCmguiOutput()
{
    if (this->mOutputDirectory == "")
    {
        EXCEPTION("No output directory was given so no output was written, cannot convert to cmgui format");
    }

    CmguiDeformedSolutionsWriter<DIM> writer(this->mOutputDirectory + "/cmgui",
                                             "solution",
                                             this->mrQuadMesh,
                                             WRITE_QUADRATIC_MESH);

    std::vector<c_vector<double,DIM> >& r_deformed_positions = this->rGetDeformedPosition();
    writer.WriteInitialMesh(); // this writes solution_0.exnode and .exelem
    writer.WriteDeformationPositions(r_deformed_positions, 1); // this writes the final solution as solution_1.exnode
    writer.WriteCmguiScript(); // writes LoadSolutions.com
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::SetComputeAverageStressPerElementDuringSolve(bool setComputeAverageStressPerElement)
{
    mSetComputeAverageStressPerElement = setComputeAverageStressPerElement;
    if (setComputeAverageStressPerElement && mAverageStressesPerElement.size()==0)
    {
        mAverageStressesPerElement.resize(this->mrQuadMesh.GetNumElements(), zero_vector<double>(DIM*(DIM+1)/2));
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::AddStressToAverageStressPerElement(c_matrix<double,DIM,DIM>& rT, unsigned elemIndex)
{
    assert(mSetComputeAverageStressPerElement);
    assert(elemIndex<this->mrQuadMesh.GetNumElements());

    // In 2d the matrix is
    // [T00 T01]
    // [T10 T11]
    // where T01 = T10. We store this as a vector
    // [T00 T01 T11]
    //
    // Similarly, for 3d we store
    // [T00 T01 T02 T11 T12 T22]
    for (unsigned i=0; i<DIM*(DIM+1)/2; i++)
    {
        unsigned row;
        unsigned col;
        if (DIM == 2)
        {
            row = i<=1 ? 0 : 1;
            col = i==0 ? 0 : 1;
        }
        else // DIM == 3
        {
            row = i<=2 ? 0 : (i<=4? 1 : 2);
            col = i==0 ? 0 : (i==1 || i==3? 1 : 2);
        }

        this->mAverageStressesPerElement[elemIndex](i) += rT(row,col);
    }
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> AbstractNonlinearElasticitySolver<DIM>::GetAverageStressPerElement(unsigned elementIndex)
{
    if (!mSetComputeAverageStressPerElement)
    {
        EXCEPTION("Call SetComputeAverageStressPerElementDuringSolve() before solve if calling GetAverageStressesPerElement()");
    }
    assert(elementIndex<this->mrQuadMesh.GetNumElements());

    c_matrix<double,DIM,DIM> stress;

    // In 2d the matrix is
    // [T00 T01]
    // [T10 T11]
    // where T01 = T10, and was stored as
    // [T00 T01 T11]
    //
    // Similarly, for 3d the matrix was stored as
    // [T00 T01 T02 T11 T12 T22]
    if (DIM == 2)
    {
        stress(0,0) = mAverageStressesPerElement[elementIndex](0);
        stress(1,0) = stress(0,1) = mAverageStressesPerElement[elementIndex](1);
        stress(1,1) = mAverageStressesPerElement[elementIndex](2);
    }
    else
    {
        stress(0,0) = mAverageStressesPerElement[elementIndex](0);
        stress(1,0) = stress(0,1) = mAverageStressesPerElement[elementIndex](1);
        stress(2,0) = stress(0,2) = mAverageStressesPerElement[elementIndex](2);
        stress(1,1) = mAverageStressesPerElement[elementIndex](3);
        stress(2,1) = stress(1,2) = mAverageStressesPerElement[elementIndex](4);
        stress(2,2) = mAverageStressesPerElement[elementIndex](5);
    }

    return stress;
}

///////////////////////////////////////////////////////////////////////////////////
// Methods at the 'element level'.
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::GetElementCentroidStrain(StrainType strainType,
                                                                      Element<DIM,DIM>& rElement,
                                                                      c_matrix<double,DIM,DIM>& rStrain)
{
    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;

    this->mrQuadMesh.GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    // Get the current displacement at the nodes
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[this->mProblemDimension*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    // Allocate memory for the basis functions values and derivative values
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)

    // we need the point in the canonical element which corresponds to the centroid of the
    // version of the element in physical space. This point can be shown to be (1/3,1/3).
    ChastePoint<DIM> quadrature_point;
    if (DIM == 2)
    {
        quadrature_point.rGetLocation()(0) = 1.0/3.0;
        quadrature_point.rGetLocation()(1) = 1.0/3.0;
    }
    else
    {
        assert(DIM==3);
        quadrature_point.rGetLocation()(0) = 1.0/4.0;
        quadrature_point.rGetLocation()(1) = 1.0/4.0;
        quadrature_point.rGetLocation()(2) = 1.0/4.0;
    }

    QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);

    // Interpolate grad_u
    grad_u = zero_matrix<double>(DIM,DIM);
    for (unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
    {
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned M=0; M<DIM; M++)
            {
                grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
            }
        }
    }

    c_matrix<double,DIM,DIM> deformation_gradient;

    for (unsigned i=0; i<DIM; i++)
    {
        for (unsigned M=0; M<DIM; M++)
        {
            deformation_gradient(i,M) = (i==M?1:0) + grad_u(i,M);
        }
    }

    switch(strainType)
    {
        case DEFORMATION_GRADIENT_F:
        {
            rStrain = deformation_gradient;
            break;
        }
        case DEFORMATION_TENSOR_C:
        {
            rStrain = prod(trans(deformation_gradient),deformation_gradient);
            break;
        }
        case LAGRANGE_STRAIN_E:
        {
            c_matrix<double,DIM,DIM> C = prod(trans(deformation_gradient),deformation_gradient);
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    rStrain(M,N) = 0.5* ( C(M,N)-(M==N?1:0) );
                }
            }
            break;
        }
        default:
        {
            NEVER_REACHED;
            break;
        }
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::AssembleOnBoundaryElement(
            BoundaryElement<DIM-1,DIM>& rBoundaryElement,
            c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
            c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
            bool assembleResidual,
            bool assembleJacobian,
            unsigned boundaryConditionIndex)
{
    if (this->mrProblemDefinition.GetTractionBoundaryConditionType() == PRESSURE_ON_DEFORMED
        || this->mrProblemDefinition.GetTractionBoundaryConditionType() == FUNCTIONAL_PRESSURE_ON_DEFORMED)
    {
        AssembleOnBoundaryElementForPressureOnDeformedBc(rBoundaryElement, rAelem, rBelem,
                                                         assembleResidual, assembleJacobian, boundaryConditionIndex);
        return;
    }

    rAelem.clear();
    rBelem.clear();

    if (assembleJacobian && !assembleResidual)
    {
        // Nothing to do
        return;
    }

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;
    this->mrQuadMesh.GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);

    c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> phi;

    for (unsigned quad_index=0; quad_index<this->mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * this->mpBoundaryQuadratureRule->GetWeight(quad_index);

        const ChastePoint<DIM-1>& quad_point = this->mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);

        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

        // Get the required traction, interpolating X (slightly inefficiently,
        // as interpolating using quad bases) if necessary
        c_vector<double,DIM> traction = zero_vector<double>(DIM);
        switch (this->mrProblemDefinition.GetTractionBoundaryConditionType())
        {
            case FUNCTIONAL_TRACTION:
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                for (unsigned node_index=0; node_index<NUM_NODES_PER_BOUNDARY_ELEMENT; node_index++)
                {
                    X += phi(node_index)*this->mrQuadMesh.GetNode( rBoundaryElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                traction = this->mrProblemDefinition.EvaluateTractionFunction(X, this->mCurrentTime);
                break;
            }
            case ELEMENTWISE_TRACTION:
            {
                traction = this->mrProblemDefinition.rGetElementwiseTractions()[boundaryConditionIndex];
                break;
            }
            default:
                NEVER_REACHED;
        }


        for (unsigned index=0; index<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            assert(node_index < NUM_NODES_PER_BOUNDARY_ELEMENT);

            rBelem(index) -=   traction(spatial_dim)
                             * phi(node_index)
                             * wJ;
        }
    }
}

template<unsigned DIM>
bool AbstractNonlinearElasticitySolver<DIM>::ShouldAssembleMatrixTermForPressureOnDeformedBc()
{
    if (mUseSnesSolver)
    {
        // although not using this in the first few steps might be useful when the deformation
        // is large, the snes solver is more robust, so we have this on all the time. (Also because
        // for cardiac problems and in timesteps after the initial large deformation, we want this on
        // in the first step
        return true;

        // could do something like this, if make the snes a member variable
        //PetscInt iteration_number;
        //SNESGetIterationNumber(mSnes,&iteration_number);
        //return (iteration_number >= 3);
    }
    else
    {
        return (mLastDampingValue >= 0.5);
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::AssembleOnBoundaryElementForPressureOnDeformedBc(
            BoundaryElement<DIM-1,DIM>& rBoundaryElement,
            c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
            c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
            bool assembleResidual,
            bool assembleJacobian,
            unsigned boundaryConditionIndex)
{
    assert(   this->mrProblemDefinition.GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED
           || this->mrProblemDefinition.GetTractionBoundaryConditionType()==FUNCTIONAL_PRESSURE_ON_DEFORMED);

    rAelem.clear();
    rBelem.clear();

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;
    // note: jacobian determinant may be over-written below
    this->mrQuadMesh.GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);

    ///////////////////////////////////////////////////////
    // Find the volume element of the mesh which
    // contains this boundary element
    ///////////////////////////////////////////////////////

    Element<DIM,DIM>* p_containing_vol_element = nullptr;

    std::set<unsigned> potential_elements = rBoundaryElement.GetNode(0)->rGetContainingElementIndices();
    for (std::set<unsigned>::iterator iter = potential_elements.begin();
        iter != potential_elements.end();
        iter++)
    {
        p_containing_vol_element = this->mrQuadMesh.GetElement(*iter);

        bool this_vol_ele_contains_surf_ele = true;
        // loop over the nodes of boundary element and see if they are in the volume element
        for (unsigned i=1; i<NUM_NODES_PER_BOUNDARY_ELEMENT; i++) // don't need to start at 0, given looping over contain elems of node 0
        {
            unsigned surf_element_node_index = rBoundaryElement.GetNodeGlobalIndex(i);
            bool found_this_node = false;
            for (unsigned j=0; j<p_containing_vol_element->GetNumNodes(); j++)
            {
                unsigned vol_element_node_index = p_containing_vol_element->GetNodeGlobalIndex(j);
                if (surf_element_node_index == vol_element_node_index)
                {
                    found_this_node = true;
                    break;
                }
            }
            if (!found_this_node)
            {
                this_vol_ele_contains_surf_ele = false;
                break;
            }
        }
        if (this_vol_ele_contains_surf_ele)
        {
            break;
        }
    }

    // Find the local node index in the volume element for each node in the boundary element
    std::vector<unsigned> surf_to_vol_map(NUM_NODES_PER_BOUNDARY_ELEMENT);
    for (unsigned i=0; i<NUM_NODES_PER_BOUNDARY_ELEMENT; i++)
    {
        unsigned index = rBoundaryElement.GetNodeGlobalIndex(i);
        for (unsigned j=0; j<NUM_NODES_PER_ELEMENT; j++)
        {
            if (p_containing_vol_element->GetNodeGlobalIndex(j)==index)
            {
                surf_to_vol_map[i] = j;
                break;
            }
        }
    }


    // We require the volume element to compute F, which requires grad_phi on the volume element. For this we will
    // need the inverse jacobian for the volume element
    static c_matrix<double,DIM,DIM> jacobian_vol_element;
    static c_matrix<double,DIM,DIM> inverse_jacobian_vol_element;
    double jacobian_determinant_vol_element;
    this->mrQuadMesh.GetInverseJacobianForElement(p_containing_vol_element->GetIndex(), jacobian_vol_element, jacobian_determinant_vol_element, inverse_jacobian_vol_element);

    // Get the current displacements at each node of the volume element, to be used in computing F
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[this->mProblemDimension*p_containing_vol_element->GetNodeGlobalIndex(II) + JJ];
        }
    }


    // We will need both {grad phi_i} for the quadratic bases of the volume element, for computing F..
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi_vol_element;
    // ..the phi_i for each of the quadratic bases of the surface element, for the standard FE assembly part.
    c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> quad_phi_surf_element;
    // We need this too, which is obtained by taking a subset of grad_quad_phi_vol_element
    static c_matrix<double, DIM, NUM_NODES_PER_BOUNDARY_ELEMENT> grad_quad_phi_surf_element;

    c_matrix<double,DIM,DIM> F;
    c_matrix<double,DIM,DIM> invF;

    c_vector<double,DIM> normal = rBoundaryElement.CalculateNormal();
    c_matrix<double,1,DIM> normal_as_mat;
    for (unsigned i=0; i<DIM; i++)
    {
        normal_as_mat(0,i) = normal(i);
    }

    double normal_pressure;
    switch (this->mrProblemDefinition.GetTractionBoundaryConditionType())
    {
        case PRESSURE_ON_DEFORMED:
            normal_pressure = this->mrProblemDefinition.GetNormalPressure();
            break;
        case FUNCTIONAL_PRESSURE_ON_DEFORMED:
            normal_pressure = this->mrProblemDefinition.EvaluateNormalPressureFunction(this->mCurrentTime);
            break;
        default:
            NEVER_REACHED;
    }

    for (unsigned quad_index=0; quad_index<this->mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * this->mpBoundaryQuadratureRule->GetWeight(quad_index);

        // Get the quadrature point on this surface element (in canonical space) - so eg, for a 2D problem,
        // the quad point is in 1D space
        const ChastePoint<DIM-1>& quadrature_point = this->mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);
        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quadrature_point, quad_phi_surf_element);

        // We will need the xi coordinates of this quad point in the volume element. We could do this by figuring
        // out how the nodes of the surface element are ordered in the list of nodes in the volume element,
        // however it is less fiddly to compute directly. Firstly, compute the corresponding physical location
        // of the quad point, by interpolating
        c_vector<double,DIM> X = zero_vector<double>(DIM);
        for (unsigned node_index=0; node_index<NUM_NODES_PER_BOUNDARY_ELEMENT; node_index++)
        {
            X += quad_phi_surf_element(node_index)*rBoundaryElement.GetNode(node_index)->rGetLocation();
        }


        // Now compute the xi coordinates of the quad point in the volume element
        c_vector<double,DIM+1> weight = p_containing_vol_element->CalculateInterpolationWeights(X);
        c_vector<double,DIM> xi;
        for (unsigned i=0; i<DIM; i++)
        {
            xi(i) = weight(i+1); // Note, in 2d say, weights = [1-xi(0)-xi(1), xi(0), xi(1)]
        }

        // Check one of the weights was zero, as the quad point is on the boundary of the volume element
        if (DIM == 2)
        {
            assert( DIM!=2 || (fabs(weight(0))<1e-6) || (fabs(weight(1))<1e-6) || (fabs(weight(2))<1e-6) );
        }
        else
        {
            assert( DIM!=3 || (fabs(weight(0))<1e-6) || (fabs(weight(1))<1e-6) || (fabs(weight(2))<1e-6)  || (fabs(weight(3))<1e-6) ); // LCOV_EXCL_LINE
        }

        // Now we can compute the grad_phi and then interpolate F
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(xi, inverse_jacobian_vol_element, grad_quad_phi_vol_element);

        F = identity_matrix<double>(DIM,DIM);
        for (unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned M=0; M<DIM; M++)
                {
                    F(i,M) += grad_quad_phi_vol_element(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }

        double detF = Determinant(F);
        invF = Inverse(F);

        if (assembleResidual)
        {
            c_vector<double,DIM> traction = detF*normal_pressure*prod(trans(invF),normal);

            // assemble
            for (unsigned index=0; index<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                assert(node_index < NUM_NODES_PER_BOUNDARY_ELEMENT);

                rBelem(index) -=   traction(spatial_dim)
                                 * quad_phi_surf_element(node_index)
                                 * wJ;
            }
        }

        // Sometimes we don't include the analytic jacobian for this integral
        // in the jacobian that is assembled - the ShouldAssembleMatrixTermForPressureOnDeformedBc()
        // bit below - see the documentation for this method to see why.
        if (assembleJacobian && ShouldAssembleMatrixTermForPressureOnDeformedBc())
        {
            for (unsigned II=0; II<NUM_NODES_PER_BOUNDARY_ELEMENT; II++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    grad_quad_phi_surf_element(N,II) = grad_quad_phi_vol_element(N,surf_to_vol_map[II]);
                }
            }

            static FourthOrderTensor<DIM,DIM,DIM,DIM> tensor1;
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned e=0; e<DIM; e++)
                {
                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned d=0; d<DIM; d++)
                        {
                            tensor1(N,e,M,d) = invF(N,e)*invF(M,d) - invF(M,e)*invF(N,d);
                        }
                    }
                }
            }

            // tensor2(II,e,M,d) = tensor1(N,e,M,d)*grad_quad_phi_surf_element(N,II)
            static FourthOrderTensor<NUM_NODES_PER_BOUNDARY_ELEMENT,DIM,DIM,DIM> tensor2;
            tensor2.template SetAsContractionOnFirstDimension<DIM>( trans(grad_quad_phi_surf_element), tensor1);

            // tensor3 is really a third-order tensor
            // tensor3(II,e,0,d) = tensor2(II,e,M,d)*normal(M)
            static FourthOrderTensor<NUM_NODES_PER_BOUNDARY_ELEMENT,DIM,1,DIM> tensor3;
            tensor3.template SetAsContractionOnThirdDimension<DIM>( normal_as_mat, tensor2);

            for (unsigned index1=0; index1<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index1++)
            {
                unsigned spatial_dim1 = index1%DIM;
                unsigned node_index1 = (index1-spatial_dim1)/DIM;

                for (unsigned index2=0; index2<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    rAelem(index1,index2) -=    normal_pressure
                                              * detF
                                              * tensor3(node_index2,spatial_dim2,0,spatial_dim1)
                                              * quad_phi_surf_element(node_index1)
                                              * wJ;
                }
            }
        }
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::Solve(double tol)
{
    // Check the problem definition is set up correctly (and fully).
    mrProblemDefinition.Validate();

    // If the problem includes specified pressures on deformed surfaces (as opposed
    // to specified tractions), the code needs to compute normals, and they need
    // to be consistently all facing outward (or all facing inward). Check the undeformed
    // mesh boundary elements has nodes that are ordered so that all normals are
    // outward-facing
    if (mrProblemDefinition.GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED && mCheckedOutwardNormals==false)
    {
        this->mrQuadMesh.CheckOutwardNormals();
        mCheckedOutwardNormals = true;
    }

    // Write the initial solution
    this->WriteCurrentSpatialSolution("initial", "nodes");

    if (mUseSnesSolver)
    {
        SolveSnes();
    }
    else
    {
        SolveNonSnes(tol);
    }

    // Remove pressure dummy values (P=0 at internal nodes, which should have been
    // been the result of the solve above), by linear interpolating from vertices of
    // edges to the internal node of the edge
    if (this->mCompressibilityType==INCOMPRESSIBLE)
    {
        this->RemovePressureDummyValuesThroughLinearInterpolation();
    }

    // Write the final solution
    this->WriteCurrentSpatialSolution("solution", "nodes");
}

///////////////////////////////////////
///\todo #2057 Make better choices..
///////////////////////////////////////
template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::SetKspSolverAndPcType(KSP solver)
{
    // Four alternatives
    //   (a) Petsc direct solve
    //   Otherwise iterative solve with:
    //   (b) Incompressible: GMRES with ILU preconditioner (or bjacobi=ILU on each process) [default]. Very poor on large problems.
    //   (c) Incompressible: GMRES with AMG preconditioner. Uncomment #define MECH_USE_HYPRE above. Requires Petsc3 with HYPRE installed.
    //   (d) Compressible: CG with ICC

    PC pc;
    KSPGetPC(solver, &pc);

    if (mPetscDirectSolve)
    {
        if (this->mCompressibilityType==COMPRESSIBLE)
        {
            KSPSetType(solver,KSPPREONLY);

        }
        PCSetType(pc, PCLU);

        // See #2057
        // PCFactorSetMatSolverPackage(pc,"mumps");
    }
    else
    {
        if (this->mCompressibilityType==COMPRESSIBLE)
        {
            KSPSetType(solver,KSPCG);
            if (PetscTools::IsSequential())
            {
                PCSetType(pc, PCICC);
                //Note that PetscOptionsSetValue is dangerous because we can't easily do
                //regression testing.  If a name changes, then the behaviour of the code changes
                //because it won't recognise the old name.  However, it won't fail to compile/run.
                #if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 1) //PETSc 3.1 or later
                    PetscTools::SetOption("-pc_factor_shift_type", "positive_definite");
                #else
                    PetscTools::SetOption("-pc_factor_shift_positive_definite", "");
                #endif
            }
            else
            {
                PCSetType(pc, PCBJACOBI);
            }
        }
        else
        {
            unsigned num_restarts = 100;
            KSPSetType(solver,KSPGMRES);
            KSPGMRESSetRestart(solver,num_restarts);

            #ifndef MECH_USE_HYPRE
                PCSetType(pc, PCBJACOBI); // BJACOBI = ILU on each block (block = part of matrix on each process)
            #else
                /////////////////////////////////////////////////////////////////////////////////////////////////////
                // Speed up linear solve time massively for larger simulations (in fact GMRES may stagnate without
                // this for larger problems), by using a AMG preconditioner -- needs HYPRE installed
                /////////////////////////////////////////////////////////////////////////////////////////////////////
                PetscTools::SetOption("-pc_hypre_type", "boomeramg");
                // PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
                // PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");

                PCSetType(pc, PCHYPRE);
                #if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >=2) //PETSc 3.2 or later
                    KSPSetPCSide(solver, PC_RIGHT);
                #else
                    KSPSetPreconditionerSide(solver, PC_RIGHT);
                #endif

                // other possible preconditioners..
                //PCBlockDiagonalMechanics* p_custom_pc = new PCBlockDiagonalMechanics(solver, this->mPreconditionMatrix, mBlock1Size, mBlock2Size);
                //PCLDUFactorisationMechanics* p_custom_pc = new PCLDUFactorisationMechanics(solver, this->mPreconditionMatrix, mBlock1Size, mBlock2Size);
                //remember to delete memory..
                //KSPSetPreconditionerSide(solver, PC_RIGHT);
            #endif
        }
    }
}

////////////////////////////////////////////////////////////////////
//  The code for the non-SNES solver - maybe remove all this
//  as SNES solver appears better
////////////////////////////////////////////////////////////////////

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::ComputeResidualAndGetNorm(bool allowException)
{
    if (!allowException)
    {
        // Assemble the residual
        AssembleSystem(true, false);
    }
    else
    {
        try
        {
            // Try to assemble the residual using this current solution
            AssembleSystem(true, false);
        }
        catch(Exception&)
        {
            // If fail (because e.g. ODEs fail to solve, or strains are too large for material law), return infinity
            return DBL_MAX;
        }
    }

    // Return the scaled norm of the residual
    return CalculateResidualNorm();
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::CalculateResidualNorm()
{
    double norm;

    //\todo Change to NORM_1 and remove the division by mNumDofs...
    VecNorm(this->mResidualVector, NORM_2, &norm);
    return norm/this->mNumDofs;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::VectorSum(std::vector<double>& rX,
                                                       ReplicatableVector& rY,
                                                       double a,
                                                       std::vector<double>& rZ)
{
    assert(rX.size()==rY.GetSize());
    assert(rY.GetSize()==rZ.size());
    for (unsigned i=0; i<rX.size(); i++)
    {
        rZ[i] = rX[i] + a*rY[i];
    }
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::TakeNewtonStep()
{
    if (this->mVerbose)
    {
        Timer::Reset();
    }

    /////////////////////////////////////////////////////////////
    // Assemble Jacobian (and preconditioner)
    /////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    AssembleSystem(true, true);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);
    if (this->mVerbose)
    {
        Timer::PrintAndReset("AssembleSystem");
    }

    ///////////////////////////////////////////////////////////////////
    // Solve the linear system.
    ///////////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::SOLVE);

    Vec solution;
    VecDuplicate(this->mResidualVector,&solution);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD,&solver);

#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5))
    KSPSetOperators(solver, mrJacobianMatrix, this->mPreconditionMatrix);
#else
    KSPSetOperators(solver, mrJacobianMatrix, this->mPreconditionMatrix, DIFFERENT_NONZERO_PATTERN /*in precond between successive solves*/);
#endif

    // Set the type of KSP solver (CG, GMRES etc) and preconditioner (ILU, HYPRE, etc)
    SetKspSolverAndPcType(solver);

    //PetscTools::SetOption("-ksp_monitor","");
    //PetscTools::SetOption("-ksp_norm_type","natural");

    KSPSetFromOptions(solver);
    KSPSetUp(solver);


    // Set the linear system absolute tolerance.
    // This is either the user provided value, or set to
    // max {1e-6 * initial_residual, 1e-12}
    if (mKspAbsoluteTol < 0)
    {
        Vec temp;
        VecDuplicate(this->mResidualVector, &temp);
        Vec temp2;
        VecDuplicate(this->mResidualVector, &temp2);
        Vec linsys_residual;
        VecDuplicate(this->mResidualVector, &linsys_residual);

        KSPInitialResidual(solver, solution, temp, temp2, linsys_residual, this->mLinearSystemRhsVector);
        double initial_resid_norm;
        VecNorm(linsys_residual, NORM_2, &initial_resid_norm);

        PetscTools::Destroy(temp);
        PetscTools::Destroy(temp2);
        PetscTools::Destroy(linsys_residual);

        double ksp_rel_tol = 1e-6;
        double absolute_tol = ksp_rel_tol * initial_resid_norm;
        if (absolute_tol < 1e-12)
        {
            absolute_tol = 1e-12;
        }
        KSPSetTolerances(solver, 1e-16, absolute_tol, PETSC_DEFAULT, 1000 /* max iters */); // Note: some machines - max iters seems to be 1000 whatever we give here
    }
    else
    {
        KSPSetTolerances(solver, 1e-16, mKspAbsoluteTol, PETSC_DEFAULT, 1000 /* max iters */); // Note: some machines - max iters seems to be 1000 whatever we give here
    }

    if (this->mVerbose)
    {
        Timer::PrintAndReset("KSP Setup");
    }

    KSPSolve(solver,this->mLinearSystemRhsVector,solution);

//    ///// For printing matrix when debugging
//    OutputFileHandler handler("TEMP",false);
//    std::stringstream ss;
//    static unsigned counter = 0;
//    ss << "all_" << counter++ << ".txt";
//    out_stream p_file = handler.OpenOutputFile(ss.str());
//    *p_file << std::setprecision(10);
//    for (unsigned i=0; i<this->mNumDofs; i++)
//    {
//        for (unsigned j=0; j<this->mNumDofs; j++)
//        {
//            *p_file << PetscMatTools::GetElement(mrJacobianMatrix, i, j) << " ";
//        }
//        *p_file << PetscVecTools::GetElement(this->mLinearSystemRhsVector, i) << " ";
//        *p_file << PetscVecTools::GetElement(solution, i) << "\n";
//    }
//    p_file->close();


    /////////////////////////////////////////////
    // Error checking for linear solve
    /////////////////////////////////////////////

    // warn if ksp reports failure
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver,&reason);

    if (reason != KSP_DIVERGED_ITS)
    {
        // Throw an exception if the solver failed for any reason other than DIVERGED_ITS.
        // This is not covered as would be difficult to cover - requires a bad matrix to
        // assembled, for example.
        // LCOV_EXCL_START
        KSPEXCEPT(reason);
        // LCOV_EXCL_STOP
    }
    else
    {
        // DIVERGED_ITS just means it didn't converge in the given maximum number of iterations,
        // which is potentially not a problem, as the nonlinear solver may (and often will) still converge.
        // Just warn once.
        // (Very difficult to cover in normal tests, requires relative and absolute ksp tols to be very small, there
        // is no interface for setting both of these. Could be covered by setting up a problem the solver
        // finds difficult to solve, but this would be overkill.)
        // LCOV_EXCL_START
        WARN_ONCE_ONLY("Linear solve (within a mechanics solve) didn't converge, but this may not stop nonlinear solve converging")
        // LCOV_EXCL_STOP
    }

    // quit if no ksp iterations were done
    int num_iters;
    KSPGetIterationNumber(solver, &num_iters);
    if (num_iters==0)
    {
        PetscTools::Destroy(solution);
        KSPDestroy(PETSC_DESTROY_PARAM(solver));
        EXCEPTION("KSP Absolute tolerance was too high, linear system wasn't solved - there will be no decrease in Newton residual. Decrease KspAbsoluteTolerance");
    }


    if (this->mVerbose)
    {
        Timer::PrintAndReset("KSP Solve");
        std::cout << "[" << PetscTools::GetMyRank() << "]: Num iterations = " << num_iters << "\n" << std::flush;
    }

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::SOLVE);

    ///////////////////////////////////////////////////////////////////////////
    // Update the solution
    //  Newton method:       sol = sol - update, where update=Jac^{-1}*residual
    //  Newton with damping: sol = sol - s*update, where s is chosen
    //   such that |residual(sol)| is minimised. Damping is important to
    //   avoid initial divergence.
    //
    // Normally, finding the best s from say 0.05,0.1,0.2,..,1.0 is cheap,
    // but this is not the case in cardiac electromechanics calculations.
    // Therefore, we initially check s=1 (expected to be the best most of the
    // time, then s=0.9. If the norm of the residual increases, we assume
    // s=1 is the best. Otherwise, check s=0.8 to see if s=0.9 is a local min.
    ///////////////////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::UPDATE);
    double new_norm_resid = UpdateSolutionUsingLineSearch(solution);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::UPDATE);

    PetscTools::Destroy(solution);
    KSPDestroy(PETSC_DESTROY_PARAM(solver));

    return new_norm_resid;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::PrintLineSearchResult(double s, double residNorm)
{
    if (this->mVerbose)
    {
        std::cout << "\tTesting s = " << s << ", |f| = " << residNorm << "\n" << std::flush;
    }
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::UpdateSolutionUsingLineSearch(Vec solution)
{
    double initial_norm_resid = CalculateResidualNorm();
    if (this->mVerbose)
    {
        std::cout << "\tInitial |f| [corresponding to s=0] is " << initial_norm_resid << "\n"  << std::flush;
    }

    ReplicatableVector update(solution);

    std::vector<double> old_solution = this->mCurrentSolution;

    std::vector<double> damping_values; // = {1.0, 0.9, .., 0.2, 0.1, 0.05} ie size 11
    for (unsigned i=10; i>=1; i--)
    {
        damping_values.push_back((double)i/10.0);
    }
    damping_values.push_back(0.05);
    assert(damping_values.size()==11);

    //// Try s=1 and see what the residual-norm is
    // let mCurrentSolution = old_solution - damping_val[0]*update; and compute residual
    unsigned index = 0;
    VectorSum(old_solution, update, -damping_values[index], this->mCurrentSolution);
    double current_resid_norm = ComputeResidualAndGetNorm(true);
    PrintLineSearchResult(damping_values[index], current_resid_norm);

    //// Try s=0.9 and see what the residual-norm is
    // let mCurrentSolution = old_solution - damping_val[1]*update; and compute residual
    index = 1;
    VectorSum(old_solution, update, -damping_values[index], this->mCurrentSolution);
    double next_resid_norm = ComputeResidualAndGetNorm(true);
    PrintLineSearchResult(damping_values[index], next_resid_norm);

    index = 2;
    // While f(s_next) < f(s_current), [f = residnorm], keep trying new damping values,
    // ie exit thus loop when next norm of the residual first increases
    while (    (next_resid_norm==DBL_MAX) // the residual is returned as infinity if the deformation is so large to cause exceptions in the material law/EM contraction model
            || ( (next_resid_norm < current_resid_norm) && (index<damping_values.size()) ) )
    {
        current_resid_norm = next_resid_norm;

        // let mCurrentSolution = old_solution - damping_val*update; and compute residual
        VectorSum(old_solution, update, -damping_values[index], this->mCurrentSolution);
        next_resid_norm = ComputeResidualAndGetNorm(true);
        PrintLineSearchResult(damping_values[index], next_resid_norm);

        index++;
    }

    unsigned best_index;

    if (index==damping_values.size() && (next_resid_norm < current_resid_norm))
    {
        // Difficult to come up with large forces/tractions such that it had to
        // test right down to s=0.05, but overall doesn't fail.
        // The possible damping values have been manually temporarily altered to
        // get this code to be called, it appears to work correctly. Even if it
        // didn't tests wouldn't fail, they would just be v. slightly less efficient.
        // LCOV_EXCL_START
        // if we exited because we got to the end of the possible damping values, the
        // best one was the last one (excl the final index++ at the end)
        current_resid_norm = next_resid_norm;
        best_index = index-1;
        // LCOV_EXCL_STOP
    }
    else
    {
        // else the best one must have been the second last one (excl the final index++ at the end)
        // (as we would have exited when the resid norm first increased)
        best_index = index-2;
    }

    // See documentation for SetTakeFullFirstNewtonStep()
    bool full_first_step = mTakeFullFirstNewtonStep && mFirstStep;


    // Check out best was better than the original residual-norm
    if (initial_norm_resid < current_resid_norm && !full_first_step)
    {
        // LCOV_EXCL_START
        EXCEPTION("Residual does not appear to decrease in newton direction, quitting");
        // LCOV_EXCL_STOP
    }

    // See documentation for SetTakeFullFirstNewtonStep()
    if (full_first_step)
    {
        if (this->mVerbose)
        {
            std::cout << "\tTaking full first Newton step...\n";
        }
        best_index = 0;
    }

    if (this->mVerbose)
    {
        std::cout << "\tChoosing s = " << damping_values[best_index] << "\n"  << std::flush;
    }


///// Print the maximum change in the displacement and pressure (for debugging robustness issues). Assumes incompressible
//
//    double l_inf_disp = 0.0;
//    double l_inf_pressure = 0.0;
//
//    if (this->mCompressibilityType==INCOMPRESSIBLE)
//    {
//        for (unsigned i=0; i<this->mrQuadMesh.GetNumNodes(); i++)
//        {
//            for (unsigned j=0; j<DIM; j++)
//            {
//                double value = update[(DIM+1)*i + j]*damping_values[best_index];
//                l_inf_disp = std::max(l_inf_disp, fabs(value));
//            }
//            l_inf_pressure = std::max(l_inf_pressure, fabs(update[(DIM+1)*i + DIM]*damping_values[best_index]));
//        }
//        std::cout << "l_inf_disp, l_inf_pressure = " << l_inf_disp << " " << l_inf_pressure << "\n";
//    }
//    else
//    {
//        for (unsigned i=0; i<this->mrQuadMesh.GetNumNodes(); i++)
//        {
//            for (unsigned j=0; j<DIM; j++)
//            {
//                double value = update[DIM*i + j]*damping_values[best_index];
//                l_inf_disp = std::max(l_inf_disp, fabs(value));
//            }
//        }
//        std::cout << "l_inf_disp = " << l_inf_disp << "\n";
//    }

    VectorSum(old_solution, update, -damping_values[best_index], this->mCurrentSolution);

    mLastDampingValue = damping_values[best_index];
    return current_resid_norm;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::PostNewtonStep(unsigned counter, double normResidual)
{
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::SolveNonSnes(double tol)
{
    mLastDampingValue = 0;

    if (mWriteOutputEachNewtonIteration)
    {
        this->WriteCurrentSpatialSolution("newton_iteration", "nodes", 0);
    }

    // Compute residual
    double norm_resid = ComputeResidualAndGetNorm(false);
    if (this->mVerbose)
    {
        std::cout << "\nNorm of residual is " << norm_resid << "\n";
    }

    mNumNewtonIterations = 0;
    unsigned iteration_number = 1;

    if (tol < 0) // i.e. if wasn't passed in as a parameter
    {
        tol = NEWTON_REL_TOL*norm_resid;

        // LCOV_EXCL_START // not going to have tests in cts for everything
        if (tol > MAX_NEWTON_ABS_TOL)
        {
            tol = MAX_NEWTON_ABS_TOL;
        }
        if (tol < MIN_NEWTON_ABS_TOL)
        {
            tol = MIN_NEWTON_ABS_TOL;
        }
        // LCOV_EXCL_STOP
    }

    if (this->mVerbose)
    {
        std::cout << "Solving with tolerance " << tol << "\n";
    }

    while (norm_resid > tol)
    {
        if (this->mVerbose)
        {
            std::cout <<  "\n-------------------\n"
                      <<   "Newton iteration " << iteration_number
                      << ":\n-------------------\n";
        }

        // take newton step (and get returned residual)
        mFirstStep = (iteration_number==1);
        norm_resid = TakeNewtonStep();

        if (this->mVerbose)
        {
            std::cout << "Norm of residual is " << norm_resid << "\n";
        }

        if (mWriteOutputEachNewtonIteration)
        {
            this->WriteCurrentSpatialSolution("newton_iteration", "nodes", iteration_number);
        }

        mNumNewtonIterations = iteration_number;

        PostNewtonStep(iteration_number,norm_resid);

        iteration_number++;
        if (iteration_number==20)
        {
            // LCOV_EXCL_START
            EXCEPTION("Not converged after 20 newton iterations, quitting");
            // LCOV_EXCL_STOP
        }
    }

    if (norm_resid > tol)
    {
        // LCOV_EXCL_START
        EXCEPTION("Failed to converge");
        // LCOV_EXCL_STOP
    }
}

template<unsigned DIM>
unsigned AbstractNonlinearElasticitySolver<DIM>::GetNumNewtonIterations()
{
    return mNumNewtonIterations;
}

//////////////////////////////////////////////////////////////
//  SNES version of the nonlinear solver
//////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::SolveSnes()
{
    // Set up solution guess for residuals
    Vec initial_guess;
    VecDuplicate(this->mResidualVector, &initial_guess);
    double* p_initial_guess;
    VecGetArray(initial_guess, &p_initial_guess);
    int lo, hi;
    VecGetOwnershipRange(initial_guess, &lo, &hi);
    for (int global_index=lo; global_index<hi; global_index++)
    {
        int local_index = global_index - lo;
        p_initial_guess[local_index] = this->mCurrentSolution[global_index];
    }
    VecRestoreArray(initial_guess, &p_initial_guess);
    PetscVecTools::Finalise(initial_guess);

    Vec snes_residual_vec;
    VecDuplicate(this->mResidualVector, &snes_residual_vec);

    SNES snes;

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, snes_residual_vec, &AbstractNonlinearElasticitySolver_ComputeResidual<DIM>, this);
    SNESSetJacobian(snes, mrJacobianMatrix, this->mPreconditionMatrix, &AbstractNonlinearElasticitySolver_ComputeJacobian<DIM>, this);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 4) //PETSc 3.4 or later
    SNESSetType(snes, SNESNEWTONLS);
#else
    SNESSetType(snes, SNESLS);
#endif
    SNESSetTolerances(snes, 1e-5, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 3) //PETSc 3.3
    SNESLineSearch linesearch;
    SNESGetSNESLineSearch(snes, &linesearch);
    // Use 'critical point' line search algorithm.  This was changed from 'backtracking'; see #2916
    SNESLineSearchSetType(linesearch, "cp");
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 4) //PETSc 3.4 or later
    SNESLineSearch linesearch;
    SNESGetLineSearch(snes, &linesearch);
    // Use 'critical point' line search algorithm.  This was changed from 'backtracking'; see #2916
    SNESLineSearchSetType(linesearch, "cp");
#endif

    SNESSetMaxLinearSolveFailures(snes,100);

    KSP ksp;
    SNESGetKSP(snes, &ksp);

    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1000 /* max iters */); // Note: some machines - max iters seems to be 1000 whatever we give here

    // Set the type of KSP solver (CG, GMRES etc) and preconditioner (ILU, HYPRE, etc)
    SetKspSolverAndPcType(ksp);

    if (this->mVerbose)
    {
        PetscTools::SetOption("-snes_monitor","");
    }
    SNESSetFromOptions(snes);

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscErrorCode err = SNESSolve(snes, initial_guess);
#else
    PetscErrorCode err = SNESSolve(snes, PETSC_NULL, initial_guess);
#endif

    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);

// LCOV_EXCL_START
    if (err != 0)
    {
        std::stringstream err_stream;
        err_stream << err;
        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(snes_residual_vec);
        SNESDestroy(PETSC_DESTROY_PARAM(snes));
        EXCEPTION("Nonlinear Solver failed. PETSc error code: "+err_stream.str()+" .");
    }

    if (reason < 0)
    {
        std::stringstream reason_stream;
        reason_stream << reason;
        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(snes_residual_vec);
        SNESDestroy(PETSC_DESTROY_PARAM(snes));
        EXCEPTION("Nonlinear Solver did not converge. PETSc reason code: "+reason_stream.str()+" .");
    }
// LCOV_EXCL_STOP

    PetscInt num_iters;
    SNESGetIterationNumber(snes,&num_iters);
    mNumNewtonIterations = num_iters;

    PetscTools::Destroy(initial_guess);
    PetscTools::Destroy(snes_residual_vec);
    SNESDestroy(PETSC_DESTROY_PARAM(snes));
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::ComputeResidual(Vec currentGuess, Vec residualVector)
{
    // Note: AssembleSystem() assumes the current solution is in this->mCurrentSolution and assembles
    // this->mResiduaVector and/or this->mrJacobianMatrix. Since PETSc wants us to use the input
    // currentGuess, and write the output to residualVector, we have to copy do some copies below.

    ReplicatableVector guess_repl(currentGuess);
    for (unsigned i=0; i<guess_repl.GetSize(); i++)
    {
        this->mCurrentSolution[i] = guess_repl[i];
    }
    AssembleSystem(true,false);
    VecCopy(this->mResidualVector, residualVector);
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::ComputeJacobian(Vec currentGuess, Mat* pJacobian, Mat* pPreconditioner)
{
    // Note: AssembleSystem() assumes the current solution is in this->mCurrentSolution and assembles
    // this->mResiduaVector and/or this->mrJacobianMatrix.
    // We need to copy the input currentGuess into the local mCurrentGuess.
    // We don't have to copy mrJacobianMatrix to pJacobian, which would be expensive, as they will
    // point to the same memory.

    // check Petsc data corresponds to internal Mats
    assert(mrJacobianMatrix==*pJacobian);
    assert(this->mPreconditionMatrix==*pPreconditioner);

    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    ReplicatableVector guess_repl(currentGuess);
    for (unsigned i=0; i<guess_repl.GetSize(); i++)
    {
        this->mCurrentSolution[i] = guess_repl[i];
    }

    AssembleSystem(false,true);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);
}

template<unsigned DIM>
PetscErrorCode AbstractNonlinearElasticitySolver_ComputeResidual(SNES snes,
                                                                 Vec currentGuess,
                                                                 Vec residualVector,
                                                                 void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearElasticitySolver<DIM>* p_solver = (AbstractNonlinearElasticitySolver<DIM>*)pContext;
    p_solver->ComputeResidual(currentGuess, residualVector);
    return 0;
}

template<unsigned DIM>
#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5))
    PetscErrorCode AbstractNonlinearElasticitySolver_ComputeJacobian(SNES snes,
                                                                     Vec currentGuess,
                                                                     Mat globalJacobian,
                                                                     Mat preconditioner,
                                                                     void* pContext)
    {
        // Extract the solver from the void*
        AbstractNonlinearElasticitySolver<DIM>* p_solver = (AbstractNonlinearElasticitySolver<DIM>*) pContext;
        p_solver->ComputeJacobian(currentGuess, &globalJacobian, &preconditioner);
        return 0;
    }
#else
    PetscErrorCode AbstractNonlinearElasticitySolver_ComputeJacobian(SNES snes,
                                                                     Vec currentGuess,
                                                                     Mat* pGlobalJacobian,
                                                                     Mat* pPreconditioner,
                                                                     MatStructure* pMatStructure,
                                                                     void* pContext)
    {
        // Extract the solver from the void*
        AbstractNonlinearElasticitySolver<DIM>* p_solver = (AbstractNonlinearElasticitySolver<DIM>*) pContext;
        p_solver->ComputeJacobian(currentGuess, pGlobalJacobian, pPreconditioner);
        return 0;
    }
#endif


// Constant setting definitions

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::MAX_NEWTON_ABS_TOL = 1e-7;

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::MIN_NEWTON_ABS_TOL = 1e-10;

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::NEWTON_REL_TOL = 1e-4;

#endif /*ABSTRACTNONLINEARELASTICITYSOLVER_HPP_*/
