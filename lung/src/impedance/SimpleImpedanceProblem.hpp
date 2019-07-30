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

#ifndef SIMPLEIMPEDANCEPROBLEM_HPP_
#define SIMPLEIMPEDANCEPROBLEM_HPP_

#include <complex>

#include "TetrahedralMesh.hpp"
#include "LinearSystem.hpp"
#include "TimeStepper.hpp"
#include "VtkMeshWriter.hpp"
#include "AirwayTreeWalker.hpp"

/**
 * A class for solving a one-dimensional acoustic impedance problem on a branched 1D tree.
 *
 * This is perhaps the simplest possible impedance calculation that can be performed on an
 * airway tree. Airway wall compliance, alveolar gas shunting and many other effects are neglected.
 *
 * The impedance of a single airway is given by
 *
 * Z = R + i*I*omega
 *
 * with R & I calculated assuming Pouseille’s flow and omega is the angular frequency.
 * The impedance of an acinus is given by
 *
 * Z = -i*E/omega
 *
 * where E is the elastance of the acinus.
 */
class SimpleImpedanceProblem
{
    friend class TestSimpleImpedanceProblem;

public:
    /**
     * Constructor
     *
     * Identifies the outlet node (a.k.a root of tree or the mouth end)
     *   A check is made that it is a boundary node.  We could also check that
     *   on trees with more than one bifurcation there are no boundary nodes in
     *   its 2nd generation successors.
     *
     * @param rAirwaysMesh  The mesh to calculate on
     * @param rootIndex  the global index of the root/outlet node in the mesh (defaults to node zero).
     */
    SimpleImpedanceProblem(TetrahedralMesh<1,3>& rAirwaysMesh, unsigned rootIndex=0u);

    /**
     * Destructor
     */
    ~SimpleImpedanceProblem();

    /**
     * Sets the oscillation frequency applied at the outflow/root of the tree
     *
     * This is a special case of SetFrequencies & is largely used for testing.
     *
     * @param frequency The frequency
     */
    void SetFrequency(double frequency);

    /**
     * Sets the frequencies to test impedance for.
     *
     * Note that frequencies MUST be monotonically increasing. Defaults to F = 1,2,...,30
     *
     * @param frequencies The frequencies to test impedance for.
     */
    void SetFrequencies(std::vector<double> frequencies);

    /**
     * Returns the impedance test frequencies
     *
     * @return vector of test frequencies.
     */
    std::vector<double>& rGetFrequencies();

    /**
     * Sets the elastance of the whole lung.
     *
     * The elastance will be divided by the number of acini in the airway tree
     * to provide a homogeneous acinar elastance. Should be set in kPa/mm^3
     *
     * @param elastance The elastance of the whole lung
     */
    void SetElastance(double elastance);

    /**
     *  Performs a depth first iteration over the tree to
     *  calculate total impedance
     */
    void Solve();

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
     * Recursively calculates the impedance of an element (and all its children)
     *
     * @param pElement The root element to calculate impedance for
     * @param frequency The input frequency in Hz
     * @return The impedance
     */
    std::complex<double> CalculateElementImpedance(Element<1,3>* pElement, double frequency);

    /**
     * Calculates the impedance of a single acinar unit.
     *
     * @param pNode The boundary node representing the acinus
     * @param frequency The input frequency in Hz
     * @return The impedance
     */
    std::complex<double> CalculateAcinusImpedance(Node<3>* pNode, double frequency);

    /**
     * Set the density of air.
     *
     * Should be specified in Kg/mm^3
     *
     * @param rho The density of the fluid
     */
    void SetRho(double rho);

    /**
     * Set the dynamic viscosity of air
     *
     * Should be specific in kPa s
     *
     * @param mu The dynamic viscosity of the fluid
     */
    void SetMu(double mu);

    /**
     * Returns the impedance of the airway tree in kPa.s.mm^-3.
     *
     * This is a special case of rGetImpedances() that returns the first impedance. Largely used for testing.
     * Solve must be called first!
     *
     * @return The impedance of the airway tree
     */
    std::complex<double> GetImpedance();

    /**
     * Returns the impedances of the airway tree in kPa.s.mm^-3.
     *
     * @return The impedance of the airway tree
     */
    std::vector<std::complex<double> >& rGetImpedances();

    /**
     * Tells the solver that the supplied mesh has units in milli metres rather than metres.
     */
    void SetMeshInMilliMetres();

private:
     TetrahedralMesh<1,3>& mrMesh; /**<The 1d in 3d branching tree mesh */
     unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */

     /** Allows easy traversal of the airway tree */
     AirwayTreeWalker mWalker;

     double mRho;   /**< Density of the fluid (Defaults to (air 20 degrees C) */
     double mMu;    /**< The dynamic viscosity of the fluid */
     double mH;     /**< The elastance of the whole lung */
     double mAcinarH; /**< The elastance of a single acinus */

     double mLengthScaling; /**< Scaling to allow for a mesh specified in non-SI units */

     std::vector<double> mFrequencies; /**<The applied frequency in Hz */
     std::vector<std::complex<double> > mImpedances; /**< The calculated impedance for the network */

    /**
     * Calculate the Poiseille flow resistance of an element
     *
     * @param radius Radius of the tube
     * @param length The length of the tube
     * @return The element resistance
     */
    double CalculateElementResistance(double radius, double length);

    /**
     * Calculate the Poiseille flow inertance of an element
     *
     * @param radius Radius of the tube
     * @param length The length of the tube
     * @return The element inertance
     */
    double CalculateElementInertance(double radius, double length);
};

#endif /* SIMPLEIMPEDANCEPROBLEM_HPP_ */
