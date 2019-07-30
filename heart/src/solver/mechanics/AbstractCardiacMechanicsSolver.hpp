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

#ifndef ABSTRACTCARDIACMECHANICSSOLVER_HPP_
#define ABSTRACTCARDIACMECHANICSSOLVER_HPP_

#include <map>
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "QuadraticBasisFunction.hpp"
#include "LinearBasisFunction.hpp"
#include "AbstractContractionModel.hpp"
#include "FibreReader.hpp"
#include "FineCoarseMeshPair.hpp"
#include "AbstractCardiacMechanicsSolverInterface.hpp"
#include "HeartRegionCodes.hpp"
#include "ElectroMechanicsProblemDefinition.hpp"


/**
 *  This struct is used to collect the three things that are stored for
 *  each physical quadrature point: a contraction model, the stretch at that
 *  point, and the stretch at the last time-step (used to compute
 *  stretch rate).
 */
typedef struct DataAtQuadraturePoint_
{
    AbstractContractionModel* ContractionModel; /**< Pointer to contraction model at this quadrature point */
//    bool Active;/**<whether this quad point is active or not*/
    double Stretch; /**< Stretch (in fibre direction) at this quadrature point */
    double StretchLastTimeStep; /**< Stretch (in fibre direction) at the previous timestep, at this quadrature point */

} DataAtQuadraturePoint;


/**
 *  AbstractCardiacMechanicsSolver
 *
 *  Base class to implicit and explicit cardiac mechanics solvers. Inherits from IncompressibleNonlinearElasticitySolver
 *  or CompressibleNonlinearElasticityAssembler (depending on what the template parameter ELASTICITY_SOLVER
 *  is), and also from AbstractCardiacMechanicsSolverInterface which just declares this classes
 *  main public methods.
 *
 *  Overloads AddActiveStressAndStressDerivative() which adds on the active tension term to the stress. The
 *  child classes hold the contraction models and need to implement a method for getting the active tension from
 *  the model.
 */
template<class ELASTICITY_SOLVER, unsigned DIM>
class AbstractCardiacMechanicsSolver : public ELASTICITY_SOLVER, public AbstractCardiacMechanicsSolverInterface<DIM>
{
protected:
    static const unsigned NUM_VERTICES_PER_ELEMENT = ELASTICITY_SOLVER::NUM_VERTICES_PER_ELEMENT; /**< Useful const from base class */

    /**
     *  A map from the index of a quadrature point to the data (contraction
     *  model, stretch, stretch at the last time-step) at that quad point.
     *  Note that there is no vector of all the quadrature points of the mesh;
     *  the quad point index is the index that would be obtained by looping over
     *  elements and then looping over quad points.
     *
     *  DISTRIBUTED - only holds data for the quad points within elements
     *  owned by this process.
     */
    std::map<unsigned,DataAtQuadraturePoint> mQuadPointToDataAtQuadPointMap;

    /**
     *  An iterator to the map, used to avoid having to repeatedly search the map.
     */
    std::map<unsigned,DataAtQuadraturePoint>::iterator mMapIterator;

    /** A mesh pair object that can be set by the user to inform the solver about the electrics mesh. */
    FineCoarseMeshPair<DIM>* mpMeshPair;

    /** Total number of quad points in the (mechanics) mesh */
    unsigned mTotalQuadPoints;

    /** Current time */
    double mCurrentTime;

    /** Time to which the solver has been asked to solve to */
    double mNextTime;

    /** Time used to integrate the contraction model */
    double mOdeTimestep;

    /**
     * The fibre-sheet-normal directions (in a matrix), if constant
     * (defaults to the identity, ie fibres in the X-direction, sheet in the XY plane)
     */
    c_matrix<double,DIM,DIM> mConstantFibreSheetDirections;

    /**
     * The fibre-sheet-normal directions (matrices), one for each element. Only non-NULL if SetVariableFibreSheetDirections()
     * is called, if not mConstantFibreSheetDirections is used instead
     */
    std::vector<c_matrix<double,DIM,DIM> >* mpVariableFibreSheetDirections;

    /**
     *  Whether the fibre-sheet directions that where read in where define per element or per quadrature point.
     *  Only valid if mpVariableFibreSheetDirections!=NULL
     */
    bool mFibreSheetDirectionsDefinedByQuadraturePoint;

    /** The fibre direction for the current element being assembled on */
    c_vector<double,DIM> mCurrentElementFibreDirection;

    /** The sheet direction for the current element being assembled on */
    c_vector<double,DIM> mCurrentElementSheetDirection;

    /** The sheet normal direction for the current element being assembled on */
    c_vector<double,DIM> mCurrentElementSheetNormalDirection;


    /**
     *  This class contains all the information about the electro mechanics problem (except the material law)
     */
    ElectroMechanicsProblemDefinition<DIM>& mrElectroMechanicsProblemDefinition;

    /**
     *  @return whether the solver is implicit or not (i.e. whether the contraction model depends on lambda (and depends on
     *  lambda at the current time)). For whether dTa_dLam dependent terms need to be added to the Jacbobian
     */
    virtual bool IsImplicitSolver()=0;


    /**
     *  Overloaded AddActiveStressAndStressDerivative(), which calls on the contraction model to get
     *  the active stress and add it on to the stress tensor
     *
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param elementIndex Index of the current element
     *  @param currentQuadPointGlobalIndex The index (assuming an outer loop over elements and an inner
     *    loop over quadrature points), of the current quadrature point.
     *  @param rT The stress to be added to
     *  @param rDTdE the stress derivative to be added to, assuming
     *    the final parameter is true
     *  @param addToDTdE A boolean flag saying whether the stress derivative is
     *   required or not.
     */
    void AddActiveStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                            unsigned elementIndex,
                                            unsigned currentQuadPointGlobalIndex,
                                            c_matrix<double,DIM,DIM>& rT,
                                            FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                            bool addToDTdE);


    /**
     *  Over-ridden method which sets up an internal variable in the parent class, using
     *  the provided fibre-sheet direction information.
     *
     *  @param elementIndex element global index
     *  @param currentQuadPointGlobalIndex quad point global index
     */
    void SetupChangeOfBasisMatrix(unsigned elementIndex, unsigned currentQuadPointGlobalIndex);

    /**
     * Sets relevant data at all quadrature points, including whether it is an active region or not.
     *
     * Calls a contraction cell factory to assign a model to each (quadrature point in each) element.
     */
    void Initialise();

    /**
     *  Pure method called in AbstractCardiacMechanicsSolver::AddActiveStressAndStressDerivative(), which needs to provide
     *  the active tension (and other info if implicit (if the contraction model depends on stretch
     *  or stretch rate)) at a particular quadrature point. Takes in the current fibre stretch.
     *
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex quadrature point the integrand is currently being evaluated at in AssembleOnElement
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch. Only should be set in implicit solvers
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate.  Only should be set in implicit solver
     */
    virtual void GetActiveTensionAndTensionDerivs(double currentFibreStretch,
                                                  unsigned currentQuadPointGlobalIndex,
                                                  bool assembleJacobian,
                                                  double& rActiveTension,
                                                  double& rDerivActiveTensionWrtLambda,
                                                  double& rDerivActiveTensionWrtDLambdaDt)=0;
public:
    /**
     * Constructor
     *
     * @param rQuadMesh A reference to the mesh.
     * @param rProblemDefinition Object defining body force and boundary conditions
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     */
    AbstractCardiacMechanicsSolver(QuadraticMesh<DIM>& rQuadMesh,
                                   ElectroMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                   std::string outputDirectory);

    /**
     *  Destructor
     */
    virtual ~AbstractCardiacMechanicsSolver();

    /**
     * Sets the fine-coarse mesh pair object so that the solver knows about electrics too.
     * It checks that the coarse mesh of the fine-mesh pair has the same number of elements as
     * the quad mesh of this object and throws an exception otherwise.
     *
     * @param pMeshPair the FineCoarseMeshPair object to be set
     */
    void SetFineCoarseMeshPair(FineCoarseMeshPair<DIM>* pMeshPair);

    /** @return the total number of quad points in the mesh. Pure, implemented in concrete solver */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }

    /** @return the quadrature rule used in the elements. */
    virtual GaussianQuadratureRule<DIM>* GetQuadratureRule()
    {
        return this->mpQuadratureRule;
    }

    /**
     * @return access mQuadPointToDataAtQuadPointMap. See doxygen for this variable
     */
    std::map<unsigned,DataAtQuadraturePoint>& rGetQuadPointToDataAtQuadPointMap()
    {
        return mQuadPointToDataAtQuadPointMap;
    }


    /**
     *  Set a constant fibre-sheet-normal direction (a matrix) to something other than the default (fibres in X-direction,
     *  sheet in the XY plane)
     *  @param rFibreSheetMatrix The fibre-sheet-normal matrix (fibre dir the first column, normal-to-fibre-in sheet in second
     *  column, sheet-normal in third column).
     */
    void SetConstantFibreSheetDirections(const c_matrix<double,DIM,DIM>& rFibreSheetMatrix);

    /**
     *  Set a variable fibre-sheet-normal direction (matrices), from file.
     *  If the second parameter is false, there should be one fibre-sheet definition for each element; otherwise
     *  there should be one fibre-sheet definition for each *quadrature point* in the mesh.
     *  In the first case, the file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir
     *  for that element), in the second it should have .orthoquad as the format.
     *
     *  @param rOrthoFile the file containing the fibre/sheet directions
     *  @param definedPerQuadraturePoint whether the fibre-sheet definitions are for each quadrature point in the mesh
     *   (if not, one for each element is assumed).
     */
    void SetVariableFibreSheetDirections(const FileFinder& rOrthoFile, bool definedPerQuadraturePoint);

    /**
     *  Set the intracellular Calcium concentrations and voltages at each quad point. Pure.
     *
     *  Implicit solvers (for contraction models which are functions of stretch (and maybe
     *  stretch rate) would integrate the contraction model with this Ca/V/t using the current
     *  stretch (ie inside AssembleOnElement, ie inside GetActiveTensionAndTensionDerivs).
     *  Explicit solvers (for contraction models which are NOT functions of stretch can immediately
     *  integrate the contraction models to get the active tension.
     *
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     */
    void SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations,
                              std::vector<double>& rVoltages);

    /**
     *  Solve for the deformation, integrating the contraction model ODEs.
     *
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    virtual void Solve(double time, double nextTime, double odeTimestep)=0;

    /**
     *  Compute the deformation gradient, and stretch in the fibre direction, for each element in the mesh.
     *  Note: using quadratic interpolation for position, the deformation gradients and stretches
     *  actually vary linearly in each element. However, for computational efficiency reasons, when computing
     *  deformation gradients and stretches to pass back to the electrophysiology solver, we just assume
     *  they are constant in each element (ie ignoring the quadratic correction to the displacement). This means
     *  that  the (const) deformation gradient and stretch for each element can be computed in advance and
     *  stored, and we don't have to worry about interpolation onto the precise location of the cell-model (electrics-mesh)
     *  node, just which element it is in, and ditto the electric mesh element centroid.
     *
     *  To compute this (elementwise-)constant F (and from it the constant stretch), we just have to compute
     *  F using the deformed positions at the vertices only, with linear bases, rather than all the
     *  nodes and quadratic bases.
     *
     *  @param rDeformationGradients A reference of a std::vector in which the deformation gradient
     *  in each element will be returned. Must be allocated prior to being passed in.
     *  @param rStretches A reference of a std::vector in which the stretch in each element will be returned.
     *  Must be allocated prior to being passed in.
     */
    void ComputeDeformationGradientAndStretchInEachElement(std::vector<c_matrix<double,DIM,DIM> >& rDeformationGradients,
                                                           std::vector<double>& rStretches);
};


#endif /*ABSTRACTCARDIACMECHANICSSOLVER_HPP_*/
