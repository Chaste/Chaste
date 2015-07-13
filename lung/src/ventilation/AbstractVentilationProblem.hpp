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

 };

#endif /* ABSTRACTVENTILATIONPROBLEM_HPP_ */
