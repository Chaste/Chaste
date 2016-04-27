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

#ifndef IMMERSEDBOUNDARYCELLCELLINTERACTIONFORCE_HPP_
#define IMMERSEDBOUNDARYCELLCELLINTERACTIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <iostream>

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * Energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104.
 */


template<unsigned DIM>
class ImmersedBoundaryCellCellInteractionForce : public AbstractImmersedBoundaryForce<DIM>
{
//friend class TestForces;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
//        archive & mRestLength;
//        archive & mSpringConstant;
    }

protected:

    /** The immersed boundary cell population */
    ImmersedBoundaryCellPopulation<DIM>* mpCellPopulation;

    /** The immersed boundary mesh */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** The cell-cell spring constant */
    double mSpringConst;

    /** The cell-cell rest length */
    double mRestLength;

    /** The number of transmembrane proteins represented in this force class */
    unsigned mNumProteins;

    /** Whether to use linear spring forces */
    bool mLinearSpring;

    /** Whether to use a force derived from the Morse potential */
    bool mMorse;

    /** A vector storing in which position of the node attributes vector each protein is represented */
    std::vector<unsigned> mProteinNodeAttributeLocations;

public:

    /**
     * Constructor.
     */
    ImmersedBoundaryCellCellInteractionForce(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * For serialization.
     */
    ImmersedBoundaryCellCellInteractionForce();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryCellCellInteractionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the immersed boundary cell population as a result of cell-cell interactions.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     * @return mProteinNodeAttributeLocations
     */
    const std::vector<unsigned>& rGetProteinNodeAttributeLocations() const;

    /**
     * Helper method for the constructor.
     *
     * Initializes the levels of each protein.
     */
    void InitializeProteinLevels();

    /**
     * Helper method for AddForceContribution().
     *
     * Updates the levels of each protein at each timestep.
     */
    void UpdateProteinLevels();

    /**
     * Set the spring constant.
     */
    void SetSpringConstant(double springConst);

    /**
     * Set the rest length.
     */
    void SetRestLength(double restLength);

    /**
     * Set the force law to linear spring (default)
     */
    void UseLinearSpringLaw();

    /**
     * Set the force law to be based on the Morse potential
     */
    void UseMorsePotential();

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryCellCellInteractionForce)

#endif /*IMMERSEDBOUNDARYCELLCELLINTERACTIONFORCE_HPP_*/
