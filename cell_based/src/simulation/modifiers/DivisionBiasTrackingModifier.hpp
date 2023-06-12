/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef DIVISIONBIASTRACKINGMODIFIER_HPP_
#define DIVISIONBIASTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class which at each simulation time step calculates a real number between 
 * zero and 1 for each cell, which represents the normalised distance of the cell along 
 * a given axis of the cell population through its centroid, and stores this number in 
 * the CellData property as "bias". To be used in conjunction with the 
 * BiasedBernoulliTrialCellCycleModel class.
 */
template<unsigned DIM>
class DivisionBiasTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    /**
     * The specified axis along which division probability is biased.
     * Initialized in the constructor.
     */
    c_vector<double, DIM> mDivisionBiasVector;

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
    }

public:

    /**
     * Constructor.
     * 
     * @param centre the specified axis along which division probability is biased
     */
    DivisionBiasTrackingModifier(c_vector<double, DIM> divisionBiasVector);

    /**
     * Destructor.
     */
    virtual ~DivisionBiasTrackingModifier();
  
    /**
     * @return #mDivisionBiasVector.
     */
    const c_vector<double, DIM>& rGetDivisionBiasVector() const;
    
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
     * Helper method to calculate the normalised distance of each cell in the population 
     * along the division bias axis and store these in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DivisionBiasTrackingModifier)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a DivisionBiasTrackingModifier.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const DivisionBiasTrackingModifier<DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    c_vector<double, DIM> point = t->rGetDivisionBiasVector();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << point[i];
    }
}

/**
 * De-serialize constructor parameters and initialize a DivisionBiasTrackingModifier.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, DivisionBiasTrackingModifier<DIM>* t, const unsigned int file_version)
{
    // Retrieve c_vector one component at a time
    c_vector<double, DIM> point;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> point[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)DivisionBiasTrackingModifier<DIM>(point);
}
}
} // namespace ...

#endif /*DIVISIONBIASTRACKINGMODIFIER_HPP_*/
