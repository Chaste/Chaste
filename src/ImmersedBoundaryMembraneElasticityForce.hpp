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

#ifndef IMMERSEDBOUNDARYMEMBRANEELASTICITYFORCE_HPP_
#define IMMERSEDBOUNDARYMEMBRANEELASTICITYFORCE_HPP_

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
class ImmersedBoundaryMembraneElasticityForce : public AbstractImmersedBoundaryForce<DIM>
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
    }

protected:

    /** The cell population */
    ImmersedBoundaryCellPopulation<DIM>* mpCellPopulation;

    /** The immersed boundary mesh */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** The membrane spring constant associated with each element */
    double mSpringConst;

    /** The membrane rest length associated with each element */
    double mRestLength;

    /** Vector containing locations of corner-node-indices in the element attribute vectors */
    std::vector<unsigned> mCornerLocationsInAttributeVector;

    /** Vector containing locations of apical and basal rest-lengths in the element attribute vectors */
    std::vector<unsigned> mRestLengthLocationsInAttributeVector;

public:

    /**
     * Constructor.
     */
    ImmersedBoundaryMembraneElasticityForce(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * For serialization.
     */
    ImmersedBoundaryMembraneElasticityForce();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryMembraneElasticityForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     * Splits the nodes into three categories: basal, apical, and lateral.  We keep this information in the node
     * attribute called region, with 0, 1, and 2 representing basal, apical, and lateral respectively.
     */
    void TagNodeRegions();

    /*
     * We calculate the 'corners' of each element, in order to alter behaviour on apical, lateral, and basal regions
     * separately. Corners are stored as four consecutive element attributes, numbered sequentially clockwise from the
     * apical left-hand corner.
     */
    void TagElementCorners();

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMembraneElasticityForce)

#endif /*IMMERSEDBOUNDARYMEMBRANEELASTICITYFORCE_HPP_*/
