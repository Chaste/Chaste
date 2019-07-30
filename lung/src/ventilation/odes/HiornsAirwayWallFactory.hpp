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


#ifndef HIORNSAIRWAYWALLFACTORY_HPP_
#define HIORNSAIRWAYWALLFACTORY_HPP_

#include "AbstractAirwayWallFactory.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "AirwayTreeWalker.hpp"
#include "HiornsAirwayWall.hpp"

/**
 * A factory to ease creating Hiorns based airway wall models.
 *
 * This creates Hiorns models based on Horsfield or Strahler order.
 */
class HiornsAirwayWallFactory : public AbstractAirwayWallFactory
{
public:

    /**
     * Default constructor.
     *
     * By default airway walls models will be defined by interpolating Hiorns's parameters
     * onto element Horsfield order.
     *
     * @param useStrahlerOrder Choose to use Strahler order instead of Horsfield order
     */
    HiornsAirwayWallFactory(bool useStrahlerOrder = false);

    /**
     * Destructor
     */
    virtual ~HiornsAirwayWallFactory();


    /**
     * @return a newly airway wall object for the given node.
     *
     * @param pNode  Pointer to node object.
     * @return An airway wall object parameterised for this airway
     */
    virtual HiornsAirwayWall* CreateAirwayWallForElement(Element<1,3>* pElement);

    /**
     * Creates a new Hiorns airway wall.
     *
     * Can be overridden to allow specialised subclasses of the model to be used.
     *
     * @return An airway wall object
     */
    virtual HiornsAirwayWall* CreateBasicAirwayWall();

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter mu for the given generation
     */
    double GetmuForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter phi1 for the given generation
     */
    double Getphi1ForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter phi2 for the given generation
     */
    double Getphi2ForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter C1 for the given generation
     */
    double GetC1ForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter C2 for the given generation
     */
    double GetC2ForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter A for the given generation
     */
    double GetAForGeneration(unsigned generation);

    /**
     * @param generation The generation number to get the parameter for
     * @return The parameter alpha0 for the given generation
     */
    double GetAlpha0ForGeneration(unsigned generation);


    /**
     *
     * @param time The current time in seconds
     * @param pNode Pointer to node object.
     * @return The pleural pressure at the given node at the given time
     */
    virtual double GetPleuralPressureForAirway(double time, Element<1,3>* pElement);

    /**
     * @param pMesh  the mesh for which to create acinar units.
     */
    virtual void SetMesh(AbstractTetrahedralMesh<1,3>* pMesh);

private:
    /** Walker to determine order of airway elements */
    AirwayTreeWalker* mpWalker;

    /** By default use Horsfield order, this flag allows use of Strahler order */
    bool mUseStrahlerOrder;

    /** The largest order in the airway tree*/
    double mMaxOrder;

    /** alpha_0 values by generation dimensionless */
    static const double mAlpha0[17];

    /** alpha_0_prime values by generation in Pascals^-1 */
    static const double mAlpha0Prime[17];

    /** n1 values by generation dimensionless */
    static const double mN1[17];

    /** n2 values by generation dimensionless  */
    static const double mN2[17];

    /** The number of generations. Nb generation number starts at 0! */
    static const double mMaxGeneration;

    /** mu values by generation dimensionless  */
    static const double mmu[17];

    /** phi1 values by generation dimensionless  */
    static const double mphi1[17];

    /** phi2 values by generation dimensionless  */
    static const double mphi2[17];

    /** C1 values by generation dimensionless  */
    static const double mC1[17];

    /** C2 values by generation dimensionless  */
    static const double mC2[17];
};

#endif /*HIORNSAIRWAYWALLFACTORY_HPP_*/
