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

#ifndef VENTILATIONPROBLEM_HPP_
#define VENTILATIONPROBLEM_HPP_

#include <map>
#include <petscsnes.h>
#include "AbstractVentilationProblem.hpp"
#include "LinearSystem.hpp"
#include "TimeStepper.hpp"
#include "VtkMeshWriter.hpp"
/**
 * A class for solving one-dimensional flow in pipe problems on branching trees.
 *
 * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
 * At graph nodes: the flux is balanced so that mass is conserved.
 *
 * Works in 3D <1,3>
 * Current functionality: pressure boundary conditions are set on each of the boundary nodes
 * Solves for pressure at internal nodes and flux on edges
 *
 * In this subclass fluxes are propagated up the tree and pressures are propagated down the tree.
 * If pressure boundary conditions are given at the bottom of the tree then this is done iteratively
 * (a KSP matrix solution is used to estimate corrections to the fluxes on the terminal edges).
 */
class VentilationProblem : public AbstractVentilationProblem
{
private:
    friend class TestVentilationProblem;
    friend PetscErrorCode ComputeSnesResidual(SNES snes, Vec terminal_flux_solution, Vec terminal_pressure_difference, void* pContext);

    /**< Used to hold the flux solution (and boundary conditions) in edge index ordering */
    std::vector<double> mFlux;

    /**< Used to hold the pressure solution (and outlet boundary pressure) in node index ordering). */
    std::vector<double> mPressure;

    /**< Pressure boundary conditions at terminal nodes.*/
    std::map<unsigned, double> mPressureCondition;

    /**< Used to switch solution methods.  If the flux is given at the boundary then the entire system can be solved directly by back substitution. */
    bool mFluxGivenAtInflow;

    /**< \todo #2300 Used to switch solution methods.  If the flux is given at the top of the tree then we can no longer reflect off a known pressure condition and need to adjust the top pressure to match flux. */
    bool mFluxGivenAtOutflow;

    /**
     * The symmetric matrix is an estimate of the dense matrix system which determines how flux changes at terminal
     * nodes are reflected in pressure changes at terminals in the form
     *  P_{terminal} = A Q_{terminal}.
     * Each entry in the dense matrix (for Poiseuille flow) is the sum of the resistances of all pipes which are
     * common ancestors of a pair of terminals.
     * In order to iteratively match pressure conditions we must invert this equation.
     */
    Mat mTerminalInteractionMatrix;

    /**
     * Number of non-zeroes per row in #mTerminalInteractionMatrix (which is a sparse approximation to a fully dense matrix).
     * This is used to cut down the fill of the matrix.   The time spent constructing the approximation goes up
     * quadratically in this parameter.
     *
     * A sensible number for this parameter is 25 (or somewhere between 10 and 50).
     *
     * Why? Because when running profiled time-varying Poiseuille and Pedley simulations on large airway trees:
     *  * If mNumNonZeroesPerRow is very large (~500) then the time to construct the becomes infeasible.
     *  * If mNumNonZeroesPerRow is large (250) then the iteration count is quite small but the cost of doing in
     * PETSc solve dominates the run time
     *  * If mNumNonZeroesPerRow is between 10 and 50 then the cost of doing the linear solve and the cost of doing the
     * direct solve on each iteration are balanced.  The total run time is mimimised.
     *  * If mNumNonZeroesPerRow is 1 (the terminal unit are modelled as long pipes with no interaction between them) then the
     * matrix solve is trivial but the iteration count goes too high.
     *
     */
    unsigned mNumNonZeroesPerRow;
    /**
     * The set of all terminals which are the descendants of each particular edge.
     * The root edge has all terminals as its descendants.  Maximum depth edges have only one descendant.
     *
     * This structure is used in the computation of mTerminalInteractionMatrix
     */
    std::vector<std::set<unsigned> > mEdgeDescendantNodes;

    /** A mapping from the indexing scheme used in the mTerminalInteractionMatrix to the full mesh node indexing */
    std::map<unsigned, unsigned> mTerminalToNodeIndex;

    /** A mapping from the indexing scheme used in the mTerminalInteractionMatrix to the full mesh element indexing */
    std::map<unsigned, unsigned> mTerminalToEdgeIndex;

    /** Used as the output of the mTerminalInteractionMatrix terminal pressure to flux solver*/
    Vec mTerminalFluxChangeVector;

    /** Used as the input of the mTerminalInteractionMatrix terminal pressure to flux solver*/
    Vec mTerminalPressureChangeVector;

    /** The linear solver for the mTerminalInteractionMatrix terminal pressure to flux solver*/
    KSP mTerminalKspSolver;

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
     * This creates but doesn't fill a PETSc Mat.  It also creates two PETSc Vecs and a KSP solver.
     */
    void SetupIterativeSolver();
    /**
     * Set up the PETSc machinery for solving the iterative problem: given pressure conditions at the terminals
     * guess and refine flux conditions which match them.
     *
     * This fills the PETSc Mat created by #SetupIterativeSolver.
     *
     * @param  redoExisting  Indicates that existing resistances are changing (due to changed flux).
     */
    void FillInteractionMatrix(bool redoExisting);


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
     * Common code used by constructors
     *
     */
    void Initialise();


public: ///\todo #2300
    /** Experimental code
     *  See SolveIterativelyFromPressure() documentation
     */
    void SolveFromPressureWithSnes();

public:
    /** Default constructor
     * Attempts to read all parameters from a hard-coded file
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.
     * Creates a linear system of the appropriate size to match the mesh
     */
    VentilationProblem();

    /** Constructor
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.
     * Creates a linear system of the appropriate size to match the mesh
     *
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

    /**
     * Sets the flux at outflow/trachea/top of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (m^3)/s
     */
    void SetOutflowFlux(double flux);

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
     * @param flux The flux boundary condition in (m^3)/s
     */
    void SetFluxAtBoundaryNode(const Node<3>& rNode, double flux);

    /**
     *  Solve the system either
     *   * directly from fluxes
     *   * iteratively from pressures
     */
    void Solve();


    /**
     * @param rFluxesOnEdges The fluxes ordered by edge index (this vector is resized)
     * @param rPressuresOnNodes The pressures ordered by node index  (this vector is resized)
     */
    void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes);
};


#endif /* VENTILATIONPROBLEM_HPP_ */

