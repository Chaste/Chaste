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
#ifndef ABSTRACTVENTILATIONPROBLEM_HPP_
#define ABSTRACTVENTILATIONPROBLEM_HPP_

#include "TetrahedralMesh.hpp"

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

    /**< Use dynamic (flux related) resistance and a nonlinear solver */
    bool mDynamicResistance;

    /**< False by default (conical pipes with radius defined at nodes).  When true pipes are cylindrical.*/
    bool mRadiusOnEdge;

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



 };

#endif /* ABSTRACTVENTILATIONPROBLEM_HPP_ */
