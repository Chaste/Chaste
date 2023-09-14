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

#ifndef EXTRINSICPULLMODIFIER_HPP_
#define EXTRINSICPULLMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class which at each simulation time step applies an extrinsic 
 * pull to the right at a specified speed, optionally applying this to all 
 * nodes or else just the right-most nodes.
 */
template<unsigned DIM>
class ExtrinsicPullModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

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
        archive & mApplyExtrinsicPullToAllNodes;
        archive & mSpeed;
    }

    /**
     * Whether to apply extrinsic pull on all nodes (or just the right-most nodes).
     */
    bool mApplyExtrinsicPullToAllNodes;

    /**
     * The speed of the extrinsic pull.
     */
    double mSpeed;

public:

    /**
     * Default constructor.
     */
    ExtrinsicPullModifier();

    /**
     * Destructor.
     */
    virtual ~ExtrinsicPullModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * @param applyExtrinsicPullToAllNodes whether to apply the extrinsic pull to all nodes in the tissue
     */
    void SetApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes);

    /**
     * @param speed the speed of the extrinsic pull
     */
    void SetSpeed(double speed);
    
    /**
     * @return #mApplyExtrinsicPullToAllNodes
     */
    bool GetApplyExtrinsicPullToAllNodes();

    /**
     * @return #mSpeed
     */
    double GetSpeed();

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtrinsicPullModifier)

#endif /*EXTRINSICPULLMODIFIER_HPP_*/
