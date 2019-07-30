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

#ifndef GENERALPLANESTIMULUSCELLFACTORY_HPP_
#define GENERALPLANESTIMULUSCELLFACTORY_HPP_

#include "PlaneStimulusCellFactory.hpp"

/**
 * GeneralPlaneStimulusCellFactory
 *
 * Supplies cells with a stimuli that depend upon the number of cells and width of the mesh.
 *
 * Applied to cells within 1e-5 of x=0.
 */
template <class CELL, unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class GeneralPlaneStimulusCellFactory : public PlaneStimulusCellFactory<CELL,ELEMENT_DIM,SPACE_DIM>
{

public:
    /**
     * Constructor
     *
     * @param numEleAcross  Number of elements across which to apply the stimulus
     * @param meshWidth  Width of the mesh (used to calculate magnitude of stimulus)
     * @param stimulusMagnitude  Magnitude of the applied stimulus (defaults to -1e7, modified in the constructor dependent on mesh size).
     * @param stimulusDuration  Duration of the applied stimulus (defaults to 0.5ms).
     */
    GeneralPlaneStimulusCellFactory(unsigned numEleAcross, double meshWidth, double stimulusMagnitude=-1e7, double stimulusDuration=0.5)
        : PlaneStimulusCellFactory<CELL,ELEMENT_DIM, SPACE_DIM>(stimulusMagnitude,stimulusDuration) // These values are overridden below anyway.
    {
        stimulusMagnitude*=numEleAcross/(64.0);
        // ELEMENT_DIM==1 Justification: elements go half size with each refinement
        // ELEMENT_DIM==2 Justification: Triangles go quarter size with each refinement, but there are twice as many nodes on boundary
        // ELEMENT_DIM==3 Hypothesis: Triangles go eighth size with each refinement, but there are four-times as many nodes on boundary
        stimulusMagnitude*=meshWidth/(0.2);

        this->mpStimulus.reset(new SimpleStimulus(stimulusMagnitude, stimulusDuration));
        LOG(1, "Defined a GeneralPlaneStimulusCellFactory<"<<SPACE_DIM<<"> with SimpleStimulus("<<stimulusMagnitude<<","<< stimulusDuration<< ")\n");
    }
};

#endif /*GENERALPLANESTIMULUSCELLFACTORY_HPP_*/
