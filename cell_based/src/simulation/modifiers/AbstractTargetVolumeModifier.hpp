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

#ifndef ABSTRACTTARGETVOLUMEMODIFIER_HPP_
#define ABSTRACTTARGETVOLUMEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class in which the "target volume" CellData property of each cell is updated. It is
 * used to implement growth in 3D (monolayer) vertex-based simulations, where the volume energy term
 * of the monolayer forces pulls each cell's actual volume towards its target/equilibrium volume: as
 * the target grows over the cell cycle the cell grows, and dividing daughters restart from a smaller
 * target. This is the volume analogue of AbstractTargetAreaModifier (the "target area" used by 2D
 * vertex forces), following the equilibrium-volume growth of Okuda et al. (2013), where the actual
 * cell volume is emergent from relaxation rather than prescribed directly.
 */
template<unsigned DIM>
class AbstractTargetVolumeModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mReferenceTargetVolume;
    }

protected:

    /**
     * A reference volume. Usually this is the target volume of mature cells in the population.
     */
    double mReferenceTargetVolume;

public:

    /**
     * Default constructor.
     */
    AbstractTargetVolumeModifier();

    /**
     * Destructor.
     */
    virtual ~AbstractTargetVolumeModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Get the reference target volume.
     *
     * @return the reference target volume.
     */
    double GetReferenceTargetVolume();

    /**
     * Set the reference target volume. The standard value is 1.0.
     *
     * @param referenceTargetVolume the new value of #mReferenceTargetVolume
     */
    void SetReferenceTargetVolume(double referenceTargetVolume);

    /**
     * Helper method to update the target volume property of all cells in the population.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateTargetVolumes(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to update the target volume property of an individual cell.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell pointer to the cell
     */
    virtual void UpdateTargetVolumeOfCell(const CellPtr pCell)=0;

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractTargetVolumeModifier)

#endif /*ABSTRACTTARGETVOLUMEMODIFIER_HPP_*/
