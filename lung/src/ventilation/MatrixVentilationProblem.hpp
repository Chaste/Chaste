/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef MatrixVentilationProblem_HPP_
#define MatrixVentilationProblem_HPP_

//#define LUNG_USE_UMFPACK 1 //Uncomment to use a direct solver

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
class MatrixVentilationProblem
{
private:
    TetrahedralMesh<1,3> mMesh; /**< The 1d in 3d branching tree mesh */
    unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */
    LinearSystem* mpLinearSystem; /**< Linear system for pressure (at nodes) and flux (in edges).  Allocated by constructor */
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
     *
     *  Swan has a typo in the paper, actual value is 1.15e-6 g/mm^3 not 1.51e-6 g/mm^3!
     *  rho = 1.15    kg/m^3 <-- USED HERE
     *      = 1.15e-6 g/mm^3
     *      = 1.15e-9 kg/mm^3
     *      = 1.15e-6 kg/(m s^2) * s^2 / (mm)^2
     *      - 1.15e-6 Pa s^2 /mm^2
     *  This is used in the dynamic (Pedley) resistance calculation
     */
    double mDensity;
    double mFluxScaling;  /**< In order to keep the pressure and flux solution at a comparable magnitude, so solve for mFluxScaling * flux.  This should be the same scale as Poiseuille resistance (comparable to viscosity).*/
    Vec mSolution; /**< Allow access to the solution of the linear system and use as a guess later */

    std::vector<Swan2012AcinarUnit*> mAcinarUnits; /**< One acinar unit for each terminal node. \todo These will be abstract*/

    double mLengthScaling; /**< This solver is designed to be used with SI units, but meshes in mm are common. This scaling allows this to be handled.*/

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
     *   A check is made that it is a boundary node.  We could also check that
     *   on trees with more than one bifurcation there are no boundary nodes in
     *   its 2nd generation successors.
     * Creates a linear system of the appropriate size to match the mesh
     */
    MatrixVentilationProblem();

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
    MatrixVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex=0u);
    /** Destructor
     *  destroys the linear system
     */
    ~MatrixVentilationProblem();

    /**
     * Tells the solver that the supplied mesh has units in milli metres rather than metres.
     */
    void SetMeshInMilliMetres();

    /** Sets the pressure at the outflow/root of the tree
     *
     * @param pressure  The pressure value in Pascals
     */
    void SetOutflowPressure(double pressure);

    /**
     * Sets the flux at each inflow/leaf-edge of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (m^3)/s
     */
    void SetOutflowFlux(double flux);

    /** Sets the pressure at each inflow/leaf of the tree
     * @param pressure  The pressure value in Pascals
     */
    void SetConstantInflowPressures(double pressure);

    /** Sets the flux at each inflow/leaf-edge of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (m^3)/s
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
     * Sets a Dirichlet flux boundary condition for a given node.
     *
     * The given boundary condition will be applied at the next time step and persist through
     * time unless overwritten.
     *
     * @param rNode The node to set the boundary condition for
     * @param flux The flux boundary condition in (mm^3)/s
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
     * @return the PETSc solution vector (for both node pressures and edge fluxes)
     */
    Vec GetSolution();

    /**
     * The mSolution Vec is a mixture of flux and pressure solutions and, in parallel, it is distributed across
     * processors.  This method replicates the solution across all processors and then splits it into its flux
     * and pressure components.  Because of the replication it makes sense to get both solutions in a single call.
     *
     * @param rFluxesOnEdges The component of the mSolution Vec which represents fluxes (this vector is resized)
     * @param rPressuresOnNodes The component of the mSolution Vec which represents pressures (this vector is resized)
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
    void Solve(TimeStepper& rTimeStepper, void (*pBoundaryConditionFunction)(MatrixVentilationProblem*, double, const Node<3>&), const std::string& rDirName, const std::string& rFileBaseName);

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

};

#endif /* MatrixVentilationProblem_HPP_ */
