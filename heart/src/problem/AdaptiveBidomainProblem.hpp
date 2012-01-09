/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

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



#ifndef ADAPTIVEBIDOMAINPROBLEM_HPP_
#define ADAPTIVEBIDOMAINPROBLEM_HPP_

#define CXXBLAS_H    // This stops multiple definition of blas headers (one via Chaste, one via libadaptivity/include/cxxblas.h)
                    // Might break things, no idea, will keep it here until problems appear....

#ifdef CHASTE_ADAPTIVITY

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>

// Include libadaptivity header files.
// Need to add $LIBADAPTIVITY_DIR/include/ to the include path.
#define HAVE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtk.h"
#include "../metric_field/include/DiscreteGeometryConstraints.h"
#include "../metric_field/include/ErrorMeasure.h"
#include "../adapt3d/include/Adaptivity.h"

#include "AbstractCardiacProblem.hpp"
#include "BidomainProblem.hpp"
#include "BidomainSolver.hpp"
#include "BidomainTissue.hpp"
#include "AdaptiveTetrahedralMesh.hpp"
#include "AbstractStimulusFunction.hpp"
#include "StimulusBoundaryCondition.hpp"
#include "VtkMeshReader.hpp"

#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

/**
 * Class which specifies and solves a bidomain problem using adaptivity.
 *
 * The solution vector is of the form: (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the extracellular potential at
 * node j
 *
 * Adaptivity occurs at each printing time step (to avoid having to add a variable to
 * the XML schema at this stage): this will need tidying up at some point in the future
 */
class AdaptiveBidomainProblem : public BidomainProblem<3>
{
    friend class TestAdaptiveBidomainProblem;
    friend class TestAdaptiveBidomainProblemNightly;

private:

    /**
     *  We need to save the solver that is being used to switch off the
     *  electrodes (by adding default boundary conditions to the solver)
     */
    BidomainSolver<3,3>* mpSolver;

    /** Adaptive tetrahedral mesh: used to interface with adaptivity library */
    AdaptiveTetrahedralMesh *mpAdaptiveMesh;

    bool mIsMeshAdapting;     /**< Whether the mesh is to be adapted or not */

    bool mInitializeFromVtu;    /**< Whether or not we are initializing from a .vtu file */

    /** Whether we are applying Neumann bounadry conditions or not (default = false) */
    bool mUseNeumannBoundaryConditions;
    /** The value i when applying Neumann boundary conditions to x_i=a and x_i=b (a<b) */
    double mNeumannStimulusIndex;
    /** The value a when applying Neumann boundary conditions to x_i=a and x_i=b (a<b) */
    double mNeumannStimulusLowerValue;
    /** The value b when applying Neumann boundary conditions to x_i=a and x_i=b (a<b) */
    double mNeumannStimulusUpperValue;
    /** Magnitude of any Neumann stimulus */
    double mNeumannStimulusMagnitude;
    /** Duration of any Neumann stimulus */
    double mNeumannStimulusDuration;
    /** Period of any repeating Neumannn stimulus */
    double mNeumannStimulusPeriod;
    /** Neumann stimulus function */
    AbstractStimulusFunction* mpNeumannStimulus;
    /** Neumann boundary condition */
    StimulusBoundaryCondition<3>* mpNeumannStimulusBoundaryCondition;

//    /**
//     * Determine whether or not an edge of the mesh is of sufficient quality. Edge is "good" if the error metric
//     * associated  with it is in the range [ 1 - mGoodEdgeRange , 1 + mGoodEdgeRange ]
//     */
//    double mGoodEdgeRange;
//
//    /** Proportion of edges that must be deemed "bad" (i.e. not good) before an adapt takes place */
//    double mBadEdgeCriterion;

    /** Directory in which output files should be saved. */
    std::string mOutputDirectory;

    /** Filename prefix for output files */
    std::string mOutputFilenamePrefix;

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the tissue should
     * create cells.
     *
     * @param hasBath Whether the simulation has a bath (if this is true, all elements with
     * attribute = 1 will be set to be bath elements (the rest should have
     * attribute = 0)).
     */
    AdaptiveBidomainProblem(AbstractCardiacCellFactory<3>* pCellFactory, bool hasBath=false);

    /**
     * Destructor
     */
    ~AdaptiveBidomainProblem();

    /**
     * Method to stop any further adapts taking place
     */
    void DoNotAdaptMesh();

//
//     * Set the values of mGoodEdgeRange and mBadEdgeCriterion to be used by the AdaptiveTetrahedralMesh
//     * in determining whether an adapt is necessary or if the current mesh is of sufficient quality.
//     *
//     * @param range Value of mGoodEdgeRange to be used
//     * @param criterion Value of mBadEdgeCriterion to be used
//
//    void SetAdaptCriterion(double range, double criterion);

    /**
     * Add the values of all state variables into the adaptive mesh (as VTK node point data)
     * so that they can (a) be used to drive mesh adaption; and/or (b) be interpolated onto
     * the new mesh nodes.
     *
     * @param solution The current solution vector
     */
    void AddCurrentSolutionToAdaptiveMesh( Vec solution );

    /**
     * Extract the values of state variables from VTK format following adaption and initialize
     * the values of mSolution and the cardiac cell models at each mesh node for simulation on
     * the new mesh
     *
     * @param reader
     */
    void InitializeSolutionOnAdaptedMesh( VtkMeshReader<3,3>* reader );

    /**
     * Adapt the mesh, recontruct the mesh, cell factory, boundary conditions container,
     * bidomain tissue and solver and call InitializeSolutionOnAdaptedMesh()
     *
     * If the adapt fails for any reason then nothing is done (i.e. we don't recreate
     * objects that remain unchanged)
     */
    void AdaptMesh();

    /**
     * Return value of mError
     */
    double GetTargetError();

    /**
     * Return value of mSigma
     */
    double GetSigma();

    /**
     * Return value of mMaxEdgeLength
     */
    double GetMaxEdgeLength();

    /**
     * Return value of mMinEdgeLength
     */
    double GetMinEdgeLength();

    /**
     * Return value of mGradation
     */
    double GetGradation();

    /**
     * Return value of mMaxMeshNodes
     */
    unsigned GetMaxMeshNodes();

    /**
     * Return value of mNumAdaptSweeps
     */
    unsigned GetNumAdaptSweeps();

    /**
     * Specify that Neumann boundary conditions will be used. Use with caution: not all Neumann boundary conditions
     * lead to a well-posed problem. Best used only on cuboid geometries, when a stimulus is applied to the face
     * x_i=a and the face x_i=b is grounded.
     *
     * @param index The value i when applying the Neumann BCs to x_i=a and x_i=b (a<b) (default = 0)
     */
    void UseNeumannBoundaryCondition( unsigned index = 0 );

    /**
     * Initialize magnitude and duration of a Neumann stimulus (must have called
     * UseNeumannBoundaryConditions() or UseNeumannBoundaryConditions( true ) for this to
     * have any effect)
     *
     * @param magnitude Value of mNeumannStimulusMagnitude
     * @param duration Value of mNeumannStimulusDuration
     * @param period Optional parameter used to specify after how long a regular stimulus is repeated
     */
    void SetNeumannStimulusMagnitudeAndDuration(double magnitude, double duration, double period = DBL_MAX);

    /**
     * Create a new boundary conditions container to handle the specified Neumann boundary
     * conditions on the mesh. Needs to be called every time a new mesh is constructed (i.e.
     * after every successful adapt)
     */
    void SetupNeumannBoundaryConditionOnMesh();

    /**
     * Restore a simulation previously saved as a .vtu file. Generates mesh and initial conditions (on V, phi
     * and ODE state variables).
     *
     * Location of .vtu file is assumed to have been set by HeartConfig::SetMesh(). The alternative (constructing
     * its own mesh) doesn't make much sense since the .vtu file must contain a mesh in order to store the
     * variable values at each point of it.
     */
    void LoadSimulationFromVtuFile();

    /**
     * Performs some checks by calling the PreSolveChecks method, creates a solver to which
     * it passes the boundary conditions specified by the user (otherwise it passes the defauls bcc),
     * then at each printing step calls the AdaptMesh method and the Solve method in the solver
     * class to solve the system up until the next printing time step (for the moment, printing time
     * steps and adapt time steps are synonymous)
     *
     * Also outputs VTK files that record the solution (and the mesh) at each printing time step
     */
    void Solve();
};

#endif // CHASTE_ADAPTIVITY

#endif /*ADAPTIVEBIDOMAINPROBLEM_HPP_*/
