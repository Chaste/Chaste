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


#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include <vector>
#include <string>
#include "UblasIncludes.hpp"

#include "AbstractCardiacCellFactory.hpp"
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "AbstractOdeBasedContractionModel.hpp"
#include "AbstractCardiacMechanicsSolver.hpp"
#include "AbstractCardiacMechanicsSolverInterface.hpp"
#include "FineCoarseMeshPair.hpp"
#include "AbstractConductivityModifier.hpp"
#include "ElectroMechanicsProblemDefinition.hpp"

/**
 * Enumeration of the possible electrics problem types
 * to be used in an EM problem.
 */
typedef enum ElectricsProblemType_
{
    MONODOMAIN,//!< MONODOMAIN
    BIDOMAIN,   //!< BIDOMAIN
    BIDOMAIN_WITH_BATH
    //EXTENDED_BIDOMAIN
} ElectricsProblemType;

/**
 *  CardiacElectroMechanicsProblem
 *
 *  For solving full electro-mechanical problems.
 *
 *  Solves a monodomain problem (diffusion plus cell models) on a (fine) electrics
 *  mesh, and a mechanics problem (finite elasticity plus contraction model) on a coarse
 *  mesh. An implicit scheme (Jon Whiteley's algorithm) be be used.
 *
 *  For solving problems on regular grids use CardiacElectroMechProbRegularGeom
 *
 *  The implicit algorithm:
 *
 *  Store the position in the electrics mesh of each quad point in the mechanics mesh
 *  For every time:
 *    Solve the monodomain problem (ie integrate ODEs, solve PDE)
 *    Get intracellular [Ca] at each electrics node and interpolate on each mechanics quad point
 *    Set [Ca] on each contraction model (one for each point)
 *    Solve static finite elasticity problem implicity
 *       - guess solution
 *       - this gives the fibre stretch and stretch rate to be set on contraction models
 *       - integrate contraction models (implicitly if NHS) to get for active tension
 *       - use this active tension in computing the stress for that guess of the deformation
 *  end
 */

template<unsigned DIM, unsigned ELEC_PROB_DIM=1>
class CardiacElectroMechanicsProblem
    : public AbstractConductivityModifier<DIM,DIM> // this only inherits from this class so it can be passed to the tissue to
                                                   // allow deformation-based altering of the conductivity
{
friend class TestAbstractContractionCellFactory;
friend class TestCardiacElectroMechanicsProblem;
friend class TestCardiacElectroMechanicsFurtherFunctionality;
friend class TestElectroMechanicsExactSolution;

protected :
    /** Either COMPRESSIBLE or INCOMPRESSIBLE */
    CompressibilityType mCompressibilityType;

    /** The cardiac problem class */
    AbstractCardiacProblem<DIM,DIM,ELEC_PROB_DIM>* mpElectricsProblem;

    /** The mechanics solver - a pointer to the part that sees the cardiac mechanics interface bit.
     * (Object pointed to is the same as with mpMechanicsSolver) */
    AbstractCardiacMechanicsSolverInterface<DIM>* mpCardiacMechSolver;
    /** The mechanics solver - a pointer to the part that sees the solid mechanics solver
     * (Object pointed to is the same as with mpCardiacMechSolver) */
    AbstractNonlinearElasticitySolver<DIM>* mpMechanicsSolver;

    /** The number of electrics timesteps per mechanics timestep */
    unsigned mNumElecTimestepsPerMechTimestep;

    /**
     * A cache for the interpolated calcium concentrations from electrics to mechanics mesh.
     * Memory is allocated within Initialise(). Filled in during Solve() and passed on to the mechanics solver
     */
    std::vector<double> mInterpolatedCalciumConcs;

    /**
     * A cache for the interpolated voltages from electrics to mechanics mesh.
     * Memory is allocated within Initialise(). Filled in during Solve() and passed on to the mechanics solver
     */
    std::vector<double> mInterpolatedVoltages;

    /** The mesh for the electrics */
    TetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /** The mesh for the mechanics */
    QuadraticMesh<DIM>* mpMechanicsMesh;

    /** Object containing information about the problem to be solved */
    ElectroMechanicsProblemDefinition<DIM>* mpProblemDefinition;

    /**Whether the mesh has a bath (non-active) region or not. False by default. */
    bool mHasBath;

    /** Class wrapping both meshes, useful for transferring information */
    FineCoarseMeshPair<DIM>* mpMeshPair;

    /** Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /** Deformation output-sub-directory */
    std::string mDeformationOutputDirectory;
    /** Whether to write any output */
    bool mWriteOutput;
    /** Whether to not write out voltages */
    bool mNoElectricsOutput;

    /** when to write output */
    static const int WRITE_EVERY_NTH_TIME = 1; //hardcoded for the time being ///\todo, allow user to set this

    /** Whether any location has been set to be watched (lots of output for that location */
    bool mIsWatchedLocation;
    /** The watched location if there is one */
    c_vector<double,DIM> mWatchedLocation;
    /** The node in the electrics mesh corresponding to the watched location */
    unsigned mWatchedElectricsNodeIndex;
    /** The node in the mechanics mesh corresponding to the watched location */
    unsigned mWatchedMechanicsNodeIndex;
    /** File where watched location info is written */
    out_stream mpWatchedLocationFile;

    /** How often to print the deformation gradients and stress to file (if ever) */
    unsigned mNumTimestepsToOutputDeformationGradientsAndStress;

    /** A vector of stretches (in the fibre direction), one for each element in the mechanics mesh */
    std::vector<double> mStretchesForEachMechanicsElement;
    /** A vector of deformation gradients (each entry a matrix), one for each element in the mechanics mesh */
    std::vector<c_matrix<double,DIM,DIM> > mDeformationGradientsForEachMechanicsElement;

    /** Somewhere to store the modified conductivity tensor */
    c_matrix<double,DIM,DIM> mModifiedConductivityTensor;

    /**
     *  Determine which node is closest to the watched location
     */
    void DetermineWatchedNodes();


    /**
     *  Write info (x, y, [z], V) for the watched node.
     *
     * @param time  Time-step now, to write out
     * @param voltage  Vm vector (this is Monodomain)
     */
    void WriteWatchedLocationData(double time, Vec voltage);


public :

    /**
     * Constructor.
     * @param compressibilityType Should be either INCOMPRESSIBLE or COMPRESSIBLE
     * @param electricsProblemType the type of electrics problem (MONODOMAIN or BIDOMAIN)
     * @param pElectricsMesh  Mesh on which to solve electrics (Monodomain)
     * @param pMechanicsMesh  Mesh (2nd order) on which to solve mechanics
     * @param pCellFactory factory to use to create cells
     * @param pProblemDefinition electro-mechanics problem definition
     * @param outputDirectory the output directory
     */
    CardiacElectroMechanicsProblem(CompressibilityType compressibilityType,
                                   ElectricsProblemType electricsProblemType,
                                   TetrahedralMesh<DIM,DIM>* pElectricsMesh,
                                   QuadraticMesh<DIM>* pMechanicsMesh,
                                   AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   ElectroMechanicsProblemDefinition<DIM>* pProblemDefinition,
                                   std::string outputDirectory);

    /**
     *  Delete allocated memory and close the watched location file
     *
     *  NOTE if SetWatchedLocation but not Initialise has been
     *  called, mpWatchedLocationFile will be uninitialised and
     *  using it will cause a seg fault. Hence the mpMechanicsMesh!=NULL
     *  it is true if Initialise has been called.
     */
    virtual ~CardiacElectroMechanicsProblem();

    /**
     *  Initialise the class. Initialises the MonodomainProblem
     *  and sets up the electrics mesh to mechanics mesh data.
     */
    void Initialise();

    /**
     *  Solve the electromechanics problem
     */
    void Solve();

    /**
     *  Short helper function
     *  @return the max of a std::vector
     *  @param vec a vector of doubles
     */
    double Max(std::vector<double>& vec);

    /** Call to not write out voltages */
    void SetNoElectricsOutput();

    /**
     *  Set a location to be watched - for which lots of output
     *  is given. Should correspond to nodes in both meshes.
     *
     *  The watched file will have rows that look like:
     *  time x_pos y_pos [z_pos] voltage Ca_i_conc.
     *
     *  NOTE: for the Calcium - assumes LUO_RUDY IS USED
     *  @param watchedLocation  location (x,y,z) in space.  Watched node is the closest to this point.
     */
    void SetWatchedPosition(c_vector<double,DIM> watchedLocation);

    /**
     *  Call this for a files containing the deformation gradients (F), evaluated at the
     *  centroids of element, and the 2nd PK stress for each element (averaged over the values
     *  at the quadrature points of that element), to be written,
     *
     *  @param timestep how often to write this. Must be a multiple of the mechanics timestep.
     */
    void SetOutputDeformationGradientsAndStress(double timestep);

    /** @return the current deformed position of the nodes */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();


    /**
     *  The implementation of the pure method defined in the base class AbstractConductivityModifier. The tissue class will
     *  call this method.
     *  @param elementIndex Index of current element
     *  @param rOriginalConductivity Reference to the original (for example, undeformed) conductivity tensor
     *  @param domainIndex Used to tailor modification to the domain. 0 = intracellular, 1 = extracellular, 2 = second intracellular (tridomain)
     *  @return Reference to a modified conductivity tensor.
     */
    c_matrix<double,DIM,DIM>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,DIM,DIM>& rOriginalConductivity, unsigned domainIndex);


    /**
     * Called in Solve() before the time loop
     */
    virtual void PrepareForSolve(){}

    /**
     * Called in Solve() at the end of every time step
     * @param counter time step
     */
    virtual void OnEndOfTimeStep(unsigned counter)
    {
    }

    /**
     * Get the mechanics solver used in the solve. Only call after a solve.
     *
     * @return The mechanics (nonlinear elasticity) PDE solver.
     */
    AbstractNonlinearElasticitySolver<DIM>* GetMechanicsSolver()
    {
        assert(mpMechanicsSolver);
        return mpMechanicsSolver;
    }
};

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
