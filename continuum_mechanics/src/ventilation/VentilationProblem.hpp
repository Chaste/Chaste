/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef VENTILATIONPROBLEM_HPP_
#define VENTILATIONPROBLEM_HPP_

#include <map>
#include "TetrahedralMesh.hpp"
#include "LinearSystem.hpp"
#include "TimeStepper.hpp"
#include "VtkMeshWriter.hpp"
#include "Swan2012AcinarUnit.hpp"

/**
 * A class for solving one-dimensional flow in pipe problems on branching trees.
 *
 * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
 * At graph nodes: the flux is balanced so that mass is conserved.
 *
 * Works in 3D <1,3>
 * Current functionality: pressure boundary conditions are set on each of the boundary nodes
 * Solves for pressure at internal nodes and flux on edges
 */
class VentilationProblem
{
private:
    TetrahedralMesh<1,3> mMesh; /**< The 1d in 3d branching tree mesh */
    unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */
    bool mDynamicResistance; /**< Use dynamic (flux related) resistance and a nonlinear solver */
    bool mRadiusOnEdge; /**< False by default (conical pipes with radius defined at nodes).  When true pipes are cylindrical.*/

    /** (Dynamic) viscosity in kg/(mm*second).
     * Default to value from Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
     *  mu = 1.92e-5 Pa*s <-- USED HERE
           = 1.92e-5 kg/(m*s)
           = 1.92e-8 kg/(mm*s) or kPa.s
     * Consider http://en.wikipedia.org/wiki/Viscosity#Air which gives
     * mu = 1.81x10^(-5) kg/(m*s) -- 1.86x10^(-5) kg/(m*s)
     */
    double mViscosity;

    /** Density in kg/(mm^3).
     *  rho (for dry air) ~ 1.2041 kg/m^3 = 1.2e-9 kg/mm^3
     *  Default to Swan (page 224)
     *  rho = 1.51e-6 g/mm^3 <-- USED HERE
     *      = 1.51e-9 kg/mm^3
     *      = 1.51e-6 kg/(m s^2) * s^2 / (mm)^2
     *      - 1.51e-6 Pa s^2 /mm^2
     *  This is used in the dynamic (Pedley) resistance calculation
     */
    double mDensity;

    std::map<unsigned, AbstractAcinarUnit*> mAcinarUnits; /**< One acinar unit for each terminal node. */
    std::vector<double> mFlux; /**< Used to hold the flux solution (and boundary conditions) in edge index ordering */
    std::vector<double> mPressure; /**< Used to hold the pressure solution (and outlet boundary pressure) in node index ordering). */
    std::map<unsigned, double> mPressureCondition; /**< Pressure boundary conditions at terminal nodes. \todo This could be a vector and/or share a map with the acinar units. */
    bool mFluxGivenAtInflow; /**< Used to switch solution methods.  If the flux is given at the boundary then the entire system can be solved directly by back substitution. */

    /**
     * The symmetric matrix is an estimate of the dense matrix system which determines how flux changes at terminal
     * nodes are reflected in pressure changes at terminals in the form
     *  P_{terminal} = A Q_{terminal}.
     * Each entry in the dense matrix (for Poiseuille flow) is the sum of the resistances of all pipes which are
     * common ancestors of a pair of terminals.
     * In order to iteratively match pressure conditions we must invert this equation.
     */
    Mat mTerminalInteractionMatrix;
    std::map<unsigned, unsigned> mTerminalToNodeIndex; /**< A mapping from the indexing scheme used in the mTerminalInteractionMatrix to the full mesh node indexing */
    std::map<unsigned, unsigned> mTerminalToEdgeIndex; /**< A mapping from the indexing scheme used in the mTerminalInteractionMatrix to the full mesh element indexing */

    Vec mTerminalFluxChangeVector; /**< Used as the output of the mTerminalInteractionMatrix terminal pressure to flux solver*/
    Vec mTerminalPressureChangeVector; /**< Used as the input of the mTerminalInteractionMatrix terminal pressure to flux solver*/
    KSP mTerminalKspSolver; /**< The linear solver for the mTerminalInteractionMatrix terminal pressure to flux solver*/

    /**
     * Use flux boundary conditions at leaves (and pressure condition at root) to perform a direct solve.
     * This involves
     *  * solving directly for parent flux up the tree (using flux balance at each node)
     *  * solving for child pressure (using Poiseuille or Pedley resistance) down the tree
     */
    void SolveDirectFromFlux();


    /**
     * Set up the PETSc machinery for solving the iterative problem: given pressure conditions at the terminals
     * guess and refine flux conditions which match them.
     *
     * This creates and fills a PETSc Mat.  It also creates two PETSc Vecs and a KSP solver.
     */
    void SetupIterativeSolver();
    /**
     * Use pressure boundary conditions at leaves (and pressure condition at root) to perform
     * convert to flux boundary conditions (assuming Poiseuille flow) and then perform a direct solve
     * with SolveDirectFromFlux().  Note that pressure to flux conversion requires accumulation of the
     * resistance down the tree since weighted_sum(resistance)*flux = delta pressure.
     *
     * Note in the mDynamicResistance case we ignore dynamic resistance when setting up the weighted sum
     * of resistances.  However, since we use dynamic resistance in SolveDirectFromFlux() the solution will
     * converge to the one which uses dynamic resistance, despite the flux corrections being slightly off.
     */
    void SolveIterativelyFromPressure();

    /**
     * Get the resistance of an edge.  This defaults to Poiseuille resistance (in which only the geometry is used.
     * Otherwise, Pedley's correction is calculated, which requires a flux to be given.
     *
     * @param rElement  The edge on which to perform this calculation
     * @param usePedley  Whether to add Pedley's increasing correction term.  Here the resistance increases
     * which the sqrt of Reynold's number (dependent on flux).
     * @param flux  The flux in the edge (used for Pedley correction).
     * @return the resistance of this element/edge
     */
    double CalculateResistance(Element<1,3>& rElement, bool usePedley=false, double flux=DBL_MAX);


public:
    /** Default constructor
     * Attempts to read all parameters from a hard-coded file
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.  We could also check that
     *   on trees with more than one bifurcation there are no boundary nodes in
     *   its 2nd generation successors.
     * Creates a linear system of the appropriate size to match the mesh
     */
    VentilationProblem();

    /** Constructor
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.  We could also check that
     *   on trees with more than one bifurcation there are no boundary nodes in
     *   its 2nd generation successors.
     * Creates a linear system of the appropriate size to match the mesh
     * @param rMeshDirFilePath  the path and root name of the .node and .edge files for the mesh
     * @param rootIndex  the global index of the root/outlet node in the mesh (defaults to node zero).
     */
    VentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex=0u);
    /** Destructor
     *  destroys the linear system
     */
    ~VentilationProblem();

    /** Sets the pressure at the outflow/root of the tree
     *
     * @param pressure  The pressure value in Pascals
     */
    void SetOutflowPressure(double pressure);

    /** Sets the pressure at each inflow/leaf of the tree
     * @param pressure  The pressure value in Pascals
     */
    void SetConstantInflowPressures(double pressure);

    /** Sets the flux at each inflow/leaf-edge of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (mm^3)/s
     */
    void SetConstantInflowFluxes(double flux);


    /**
     * Sets a Dirichlet pressure boundary condition for a given node.
     *
     * The given boundary condition will be applied at the next time step and persist through
     * time unless overwritten.
     *
     * @param rNode The node to set the boundary condition for
     * @param pressure The pressure boundary condition in Pascals
     */
    void SetPressureAtBoundaryNode(const Node<3>& rNode, double pressure);

    /**
     * Gets the most recent pressure at a boundary node
     *
     * @param rNode The node to get the pressure for.
     * @return The pressure at the node.
     */
    double GetPressureAtBoundaryNode(const Node<3>& rNode);

    /**
     * Gets the most recent flux at the outflow/outlet (mouth)
     *  @return The flux at outflow.
     */
    double GetFluxAtOutflow();
    /**
     * Sets a Dirichlet flux boundary condition for a given node.
     *
     * The given boundary condition will be applied at the next time step and persist through
     * time unless overwritten.
     *
     * @param rNode The node to set the boundary condition for
     * @param flux The flux boundary condition in (mm^3)/s
     */
    void SetFluxAtBoundaryNode(const Node<3>& rNode, double flux);

    /**
     *  Solve the linear system either
     *   * directly from fluxes
     *   * iteratively from pressures
     */
    void Solve();


    /**
     * @param rFluxesOnEdges The fluxes ordered by edge index (this vector is resized)
     * @param rPressuresOnNodes The pressures ordered by node index  (this vector is resized)
     */
    void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes);

#ifdef CHASTE_VTK

    /**
     * Add flux and pressure data to a VtkMeshWriter.
     * @param rVtkWriter  the mesh writer ready for the data
     * @param rSuffix  Suffix with which to annotate e.g. pressure_001
     */
    void AddDataToVtk(VtkMeshWriter<1, 3>& rVtkWriter, const std::string& rSuffix);

    /**
     * Output the solution to a Vtk file
     * @param rDirName A directory name relative to CHASTE_TEST_OUTPUT.
     * @param rFileBaseName The base name of the new VTK file.
     */
    void WriteVtk(const std::string& rDirName, const std::string& rFileBaseName);

#endif // CHASTE_VTK

    /**
     * Used to set mRadiusOnEdge flag.
     * This is false by default in the constructor (conic pipes with radius defined at nodes).  When true pipes are cylindrical.
     * @param isOnEdges  The new value of mRadiusOnEdge
     *
     */
    void SetRadiusOnEdge(bool isOnEdges=true);

    /**
     * @return  reference to the mesh
     */
    TetrahedralMesh<1,3>& rGetMesh();

    /** Assemble the linear system by writing in
     *  * flux balance at the nodes
     *  * Poiseuille flow in the edges
     *
     *  Solve the linear system repeatedly
     *  @param rTimeStepper  The start, end and time-step
     *  @param pBoundaryConditionFunction setter
     *  @param rDirName A directory name relative to CHASTE_TEST_OUTPUT.
     *  @param rFileBaseName The base name of the new VTK file.
     */
    void Solve(TimeStepper& rTimeStepper, void (*pBoundaryConditionFunction)(VentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&), const std::string& rDirName, const std::string& rFileBaseName);

    /**
     * Read a problem definition from a file and use then solve that problem
     *
     * @param rInFilePath  Path to file which contains the problem definition
     * @param rOutFileDir  Path to folder for output (relative to CHASTE_TEST_OUTPUT)
     * @param rOutFileName  Name for VTK output
     */
    void SolveProblemFromFile(const std::string& rInFilePath, const std::string& rOutFileDir,const std::string& rOutFileName);

    /**
     * @return the viscosity in kg/(mm*sec)
     */
    double GetViscosity() const
    {
        return mViscosity;
    }

    /**
     * @param viscosity  the viscosity in kg/(mm*sec)
     */
    void SetViscosity(double viscosity)
    {
        mViscosity = viscosity;
    }

    /**
     * @return the density in kg/(m^3)
     */
    double GetDensity() const
    {
        return mDensity;
    }

    /**
     * @param density  the density in kg/(m^3)
     */
    void SetDensity(double density)
    {
        mDensity = density;
    }

    /**
     * @param dynamicResistance
     **/
    void SetDynamicResistance(bool dynamicResistance = true)
    {
        mDynamicResistance = dynamicResistance;
    }

    /**
     * @param rNode The node to get the acinus at. Must be a boundary node!
     * @return The acinar unit at the given node.
     */
    AbstractAcinarUnit* GetAcinus(const Node<3>& rNode)
    {
        assert(mAcinarUnits.count(rNode.GetIndex()));
        return mAcinarUnits[rNode.GetIndex()];
    }

};

#endif /* VENTILATIONPROBLEM_HPP_ */
