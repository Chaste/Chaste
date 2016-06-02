/*

Copyright (c) 2005-2016, University of Oxford.
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

    /** The cell population. */
    ImmersedBoundaryCellPopulation<DIM>* mpCellPopulation;

    /** The immersed boundary mesh. */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** How far through the element attributes vector we are when this constructor is called. */
    unsigned mReferenceLocationInAttributesVector;

    /** Node region code for basal, used only by this class. */
    const static unsigned msBas = 1;

    /** Node region code for apical, used only by this class. */
    const static unsigned msApi = 2;

    /** Node region code for lateral, used only by this class. */
    const static unsigned msLat = 3;

    /**
     * The membrane spring constant associated with each element.
     *
     * Initialised to 1e6 in the constructor.
     */
    double mSpringConst;

    /**
     * The membrane rest length associated with each element.
     *
     * Initialised to 0.5 in the constructor.
     */
    double mRestLengthMultiplier;

    /**
     * The multiplicative quantity by which we alter the spring constant of the basement lamina, if present.
     *
     * Initialised to 5 in the constructor.
     */
    double mBasementSpringConstantModifier;

    /**
     * The multiplicative quantity by which we alter the rest length of the basement lamina, if present.
     *
     * Initialised to 0.5 in the constructor.
     */
    double mBasementRestLengthModifier;

    /** Whether the elements have corners tagged. */
    bool mElementsHaveCorners;

    /** Vector containing locations of apical and basal rest-lengths in the element attribute vectors. */
    std::vector<unsigned> mRestLengthLocationsInAttributeVector;

    /**
     * @param elemIndex index of the element to retrieve apical length of
     * @return apical length of the specified element
     */
    double GetApicalLengthForElement(unsigned elemIndex);

    /**
     * @param elemIndex index of the element to retrieve basal length of
     * @return basal length of the specified element
     */
    double GetBasalLengthForElement(unsigned elemIndex);

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
    void TagApicalAndBasalLengths();

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
     * Overridden AddImmersedBoundaryForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     * Set mSpringConst
     *
     * @param new value of the spring constant
     */
    void SetSpringConstant(double springConstant);

    /**
     * Set mRestLength
     *
     * @param new value of the rest length
     */
    void SetRestLengthMultiplier(double restLengthMultiplier);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMembraneElasticityForce)

#endif /*IMMERSEDBOUNDARYMEMBRANEELASTICITYFORCE_HPP_*/
