/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef SIMPLETARGETVOLUMEMODIFIER_HPP_
#define SIMPLETARGETVOLUMEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractTargetVolumeModifier.hpp"

/**
 * A target volume modifier class in which the target volume of a cell grows linearly, up to
 * mReferenceTargetVolume, over a prescribed duration.
 *
 * If used with a phase-based cell-cycle model (such as FixedG1GenerationalCellCycleModel),
 * the target volume of a cell increases linearly from the value 0.5*mReferenceTargetVolume
 * up to mReferenceTargetVolume over the course of the cell's G1 phase.
 *
 * If used with a non-phase-based cell-cycle model, the target volume of a cell increases
 * linearly from the value 0.5*mReferenceTargetVolume up to mReferenceTargetVolume while the
 * cell's age is less than mGrowthDuration.
 *
 * Here mReferenceTargetVolume and mGrowthDuration are settable member variables. The default
 * value of mReferenceTargetVolume is 1.0 and the default value of mGrowthDuration is DOUBLE_UNSET.
 *
 * Note that if mGrowthDuration is set by the user, then this value is used to prescribe target
 * volume growth as described earlier, regardless of whether a phase-based cell-cycle model is
 * present.
 *
 * This is the volume analogue of SimpleTargetAreaModifier, intended for the volume energy term of
 * the 3D monolayer forces (GeneralMonolayerVertexMeshForce and MisraForce).
 */
template<unsigned DIM>
class SimpleTargetVolumeModifier : public AbstractTargetVolumeModifier<DIM>
{
    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractTargetVolumeModifier<DIM> >(*this);
        archive & mGrowthDuration;
    }

    /**
     * The duration over which a cell's target volume increases from 0.5*mReferenceTargetVolume
     * to mReferenceTargetVolume at the start of its cell cycle. Defaults to DOUBLE_UNSET.
     * If this variable is set using SetGrowthDuration(), then it is used regardless of
     * whether a phase-based cell-cycle model is used.
     */
    double mGrowthDuration;

public:

    /**
     * Default constructor.
     */
    SimpleTargetVolumeModifier();

    /**
     * Destructor.
     */
    virtual ~SimpleTargetVolumeModifier();

    /**
     * Overridden UpdateTargetVolumeOfCell() method.
     *
     * @param pCell pointer to the cell
     */
    virtual void UpdateTargetVolumeOfCell(const CellPtr pCell);

    /**
     * @return #mGrowthDuration
     */
    double GetGrowthDuration();

    /**
     * Set #mGrowthDuration.
     *
     * @param growthDuration the new value of #mGrowthDuration
     */
    void SetGrowthDuration(double growthDuration);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleTargetVolumeModifier)

#endif /*SIMPLETARGETVOLUMEMODIFIER_HPP_*/
