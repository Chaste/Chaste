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
#ifndef ABSTRACTVENTILATIONPROBLEM_HPP_
#define ABSTRACTVENTILATIONPROBLEM_HPP_

#include "TetrahedralMesh.hpp"
#include "TimeStepper.hpp"
#include "VtkMeshWriter.hpp"

/**
 * A class for solving one-dimensional flow in pipe problems on branching trees.
 *
 * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
 * At graph nodes: the flux is balanced so that mass is conserved.
 *
 * Works in 3D <1,3>
 *
 * Current functionality: pressure boundary conditions are set on each of the boundary nodes
 * Solves for pressure at internal nodes and flux on edges
 *
 * In subclasses EITHER:
 *  * All node pressures and edge fluxes are solved simultaneously using a direct matrix solution
 *  * Fluxes are propagated up the tree and pressures are propagated down the tree.  If pressure BCs are
 *    given then this is done iteratively until the terminal pressure are matched.
 */
class AbstractVentilationProblem
{

protected:
    /**
     * A helper enumeration used to keep track of the order in which things get added to edge attributes
     */
    enum
    {
        RADIUS,         //!< Radius of edge (where given in file)
        SEGMENT_LENGTH,  //!< Length of segment when several edges are linked by non-bifurcating intermediate nodes.  This is a derived quantity and not (necessarily) the same as edge length.
        PEDLEY_CORRECTION
    };

    TetrahedralMesh<1,3> mMesh; /**< The 1d in 3d branching tree mesh */
    unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */
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
     *  Default to Swan (page 224)
     *  rho = 1.15    kg/m^3 <-- USED HERE
     *      = 1.15e-6 g/mm^3
     *      = 1.15e-9 kg/mm^3
     *      = 1.15e-6 kg/(m s^2) * s^2 / (mm)^2
     *      - 1.15e-6 Pa s^2 /mm^2
     *  This is used in the dynamic (Pedley) resistance calculation
     */
    double mDensity;

    double mLengthScaling; /**< This solver is designed to be used with SI units, but meshes in mm are common. This scaling allows this to be handled.*/

    /** Use dynamic (flux related) resistance and a nonlinear solver */
    bool mDynamicResistance;
    /**
     * When using dynamic (flux related) resistance.  It's possible to set a different Pedley resistance factor on
     * each individual edge.  This can be
     *  * Set on a generational basis (van Ertbruggen 2005)
     *  * Set in the mesh file \todo This is not yet implemented
     */
    bool mPerGenerationDynamicResistance;

    /** False by default (conical pipes with radius defined at nodes).  When true pipes are cylindrical.*/
    bool mRadiusOnEdge;

    /** True by default.  Nodes should be ordered so that parents have a lower index than there children.
     * When false this will trigger warnings in the direct-solving subclass.
     */
    bool mNodesInGraphOrder;

    /**
     * Get the resistance of an edge.  This defaults to Poiseuille resistance (in which only the geometry is used.
     * Otherwise, Pedley's correction is calculated, which requires a flux to be given.
     *
     * @param rElement  The edge on which to perform this calculation
     * @param usePedley  Whether to add Pedley's increasing correction term.  Here the resistance increases
     * with the sqrt of Reynold's number (dependent on flux).
     * @param flux  The flux in the edge (used for Pedley correction).
     * @return the resistance of this element/edge
     */
    double CalculateResistance(Element<1,3>& rElement, bool usePedley=false, double flux=DBL_MAX);
    /**
     * Common code used in all constructors.  Over-ridden in direct solver
     *
     */
    virtual void Initialise();


public:
    /** Main constructor
     * Loads a mesh from file(s)
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.  We could also check that
     *   on trees with more than one bifurcation there are no boundary nodes in
     *   its 2nd generation successors.
     * Creates a linear system of the appropriate size to match the mesh
     * @param rMeshDirFilePath  the path and root name of the .node and .edge files for the mesh
     * @param rootIndex  the global index of the root/outlet node in the mesh (defaults to node zero).
     */
    AbstractVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex=0u);

    /** Virtual destructor is empty. */
    virtual ~AbstractVentilationProblem()
    {

    }
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
     * Tells the solver that the supplied mesh has units in milli metres rather than metres.
     * Overridden so we can scale other quantities in subclasses when working in SI units.
     */
    virtual void SetMeshInMilliMetres();

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

    /**
     * @param dynamicResistance
     **/
    void SetDynamicResistance(bool dynamicResistance = true)
    {
        mDynamicResistance = dynamicResistance;
    }

    /**
     * Sets a Dirichlet flux boundary condition for a given node.
     *
     * The given boundary condition will be applied at the next time step and persist through
     * time unless overwritten.
     *
     * @param rNode The node to set the boundary condition for
     * @param flux The flux boundary condition in (m^3)/s
     */
    virtual void SetFluxAtBoundaryNode(const Node<3>& rNode, double flux)=0;

    /**
     * Sets a Dirichlet pressure boundary condition for a given node.
     *
     * The given boundary condition will be applied at the next time step and persist through
     * time unless overwritten.
     *
     * @param rNode The node to set the boundary condition for
     * @param pressure The pressure boundary condition in Pascals
     */
    virtual void SetPressureAtBoundaryNode(const Node<3>& rNode, double pressure)=0;

    /** Sets the pressure at each inflow/leaf of the tree
     * @param pressure  The pressure value in Pascals
     */
    void SetConstantInflowPressures(double pressure);

    /** Sets the flux at each inflow/leaf-edge of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (m^3)/s
     */
    void SetConstantInflowFluxes(double flux);

    /** Sets the pressure at the outflow/root of the tree
     *
     * @param pressure  The pressure value in Pascals
     */
    virtual void SetOutflowPressure(double pressure);

    /**
     * Sets the flux at outflow/trachea/top of the tree
     * Flux is "volumetric flow rate"
     * @param flux  The flux value in (m^3)/s
     */
    virtual void SetOutflowFlux(double flux);


    /**
     * @param rFluxesOnEdges The fluxes ordered by edge index (this vector is resized)
     * @param rPressuresOnNodes The pressures ordered by node index  (this vector is resized)
     */
    virtual void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)=0;

    /**
     * Read a problem definition from a file and use then solve that problem
     *
     * @param rInFilePath  Path to file which contains the problem definition
     * @param rOutFileDir  Path to folder for output (relative to CHASTE_TEST_OUTPUT)
     * @param rOutFileName  Name for VTK output
     */
    void SolveProblemFromFile(const std::string& rInFilePath, const std::string& rOutFileDir,const std::string& rOutFileName);

    /**
     *  Solve linear system using
     *  * flux balance at the nodes
     *  * Poiseuille flow in the edges
     *
     *  Solve the system repeatedly
     *  @param rTimeStepper  The start, end and time-step
     *  @param pBoundaryConditionFunction setter
     *  @param rDirName A directory name relative to CHASTE_TEST_OUTPUT.
     *  @param rFileBaseName The base name of the new VTK file.
     */
    void SolveOverTime(TimeStepper& rTimeStepper, void (*pBoundaryConditionFunction)(AbstractVentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&), const std::string& rDirName, const std::string& rFileBaseName);
    /**
     *  Solve the system at one time-point.  (Implementation dependent on the concrete class.)
     */
    virtual void Solve()=0;

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
     * Set up per element dynamic resistance in the manner prescribed by
     * van Ertbruggen et al. 2005.
     * This method walks the tree and sets the attributes for the Pedley correction factor C
     * which is related to van Ertbruggen's gamma: C = 4*sqrt(2)*gamma.
     * Later use of the Solve method will pick up these factors in altering the local resistance.
     */
    void SetPerGenerationDynamicResistance();
};

#endif /* ABSTRACTVENTILATIONPROBLEM_HPP_ */
