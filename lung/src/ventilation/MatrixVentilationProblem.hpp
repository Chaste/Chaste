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

#ifndef MATRIXVENTILATIONPROBLEM_HPP_
#define MATRIXVENTILATIONPROBLEM_HPP_

//#define LUNG_USE_KLU 1 //Uncomment to use a direct solver

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
 * In this subclass all node pressures and edge fluxes are solved simultaneously using a direct matrix solution.
 */
class MatrixVentilationProblem : public AbstractVentilationProblem
{
private:
    LinearSystem* mpLinearSystem; /**< Linear system for pressure (at nodes) and flux (in edges).  Allocated by constructor */

    double mFluxScaling;  /**< In order to keep the pressure and flux solution at a comparable magnitude, so solve for mFluxScaling * flux.  This should be the same scale as Poiseuille resistance (comparable to viscosity).*/
    Vec mSolution; /**< Allow access to the solution of the linear system and use as a guess later */


    /** Assemble the linear system by writing in
     *  * flux balance at the nodes
     *  * Poiseuille flow in the edges
     *
     * Does not set any boundary conditions - these are assumed to be set elsewhere
     * Used by solvers
     *
     * @param dynamicReassemble  reassemble the Poiseuille component of the system using dynamic resistance
     * (assuming that the system has been solved once already).
     */
    void Assemble(bool dynamicReassemble=false);

public:
    /** Default constructor
     * Attempts to read all parameters from a hard-coded file
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.
     * Creates a linear system of the appropriate size to match the mesh
     */
    MatrixVentilationProblem();

    /** Constructor
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.
     * Creates a linear system of the appropriate size to match the mesh
     *
     * @param rMeshDirFilePath  the path and root name of the .node and .edge files for the mesh
     * @param rootIndex  the global index of the root/outlet node in the mesh (defaults to node zero).
     */
    MatrixVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex=0u);
    /** Destructor
     *  destroys the linear system
     */
    ~MatrixVentilationProblem();

    /**
     * Tells the solver that the supplied mesh has units in milli metres rather than metres.
     * Overridden so we can scale other quantities in subclasses when working in SI units.
     */
    void SetMeshInMilliMetres();

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

    /** Assemble the linear system by writing in
     *  * flux balance at the nodes
     *  * Poiseuille flow in the edges
     *
     *  Solve the linear system
     */
    void Solve();


    /**
     * The mSolution Vec is a mixture of flux and pressure solutions and, in parallel, it is distributed across
     * processors.  This method replicates the solution across all processors and then splits it into its flux
     * and pressure components.  Because of the replication it makes sense to get both solutions in a single call.
     *
     * @param rFluxesOnEdges The component of the mSolution Vec which represents fluxes (this vector is resized)
     * @param rPressuresOnNodes The component of the mSolution Vec which represents pressures (this vector is resized)
     */
    void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes);
};

#endif /* MATRIXVENTILATIONPROBLEM_HPP_ */
