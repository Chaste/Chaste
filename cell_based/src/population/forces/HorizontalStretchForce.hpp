/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef HORIZONTALSTRETCHFORCE_HPP_
#define HORIZONTALSTRETCHFORCE_HPP_

#include "AbstractForce.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement;

///\todo document class (#2850)
template <unsigned DIM>
class HorizontalStretchForce : public AbstractForce<DIM>
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
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive& mForceMagnitude;
        archive& mRelativeWidth;
    }

protected:
    /**
     * Internal variable to save the magnitude.
     * Initialised to 1.0 in the constructor.
     */
    double mForceMagnitude;

    /**
     * The relative width from the boundary upon with the force is acting.
     * Initialised to 0.1 in the constructor.
     */
    double mRelativeWidth;

    /**
     * The force vector acting on mMovingPinnedElements.
     * This can be changed by SetForceVector for more general cases (like shearing).
     */
    c_vector<double, DIM> mForceVector;

    /**
     * Set of fixed elements.
     */
    std::set<VertexElement<DIM, DIM>*> mNonMovingPinnedElements;

    /**
     * Set of moving elements upon which mForceVector acts.
     */
    std::set<VertexElement<DIM, DIM>*> mMovingPinnedElements;

    /**
     * Helper function to set mForceVector when mForceMagnitude is modified.
     */
    void SetUpForceVectorUsingMagnitude();

public:
    /**
     * Constructor.
     *
     * @param ForceMagnitude initial value of mForceMagnitude
     * @param RelativeWidth initial value of mRelativeWidth
     */
    HorizontalStretchForce(const double ForceMagnitude = 1,
                           const double RelativeWidth = 0.1);

    /**
     * Destructor.
     */
    virtual ~HorizontalStretchForce();

    /**
     * Set up mMovingPinnedElements and mNonMovingPinnedElements according to mRelativeWidth.
     * Besides, it also changes all pinned elements to NoCellCycleModel.
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetUpPinnedElements(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Clear mMovingPinnedElements and mNonMovingPinnedElements.
     */
    void ClearPinnedElements();

    /**
     * Get the set of pinned elements.
     * @param isMovingElements default or true for mMovingPinnedElements, false for mNonMovingPinnedElements
     * @return return either mMovingPinnedElements and mNonMovingPinnedElements
     */
    std::set<VertexElement<DIM, DIM>*>& GetPinnedElements(
        bool isMovingElements = true);

    /**
     * Set mForceVector for more general cases (like shearing).
     * @param ForceVector new value for mForceVector
     */
    void SetForceVector(const c_vector<double, DIM>& ForceVector);

    /**
     * @return reference to mForceVector
     */
    c_vector<double, DIM>& rGetForceVector();

    /**
     * Set mForceMagnitude.
     *
     * @param forceMagnitude the new value of mForceMagnitude
     */
    void SetForceMagnitude(double forceMagnitude);

    /**
     * Set mRelativeWidth.
     *
     * @param relativeWidth the new value of mRelativeWidth
     */
    void SetRelativeWidth(double relativeWidth);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(
        AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HorizontalStretchForce)

#endif /*HORIZONTALSTRETCHFORCE_HPP_*/
