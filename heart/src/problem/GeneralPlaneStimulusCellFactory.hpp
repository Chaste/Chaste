/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
     * \todo The useMeshWidth is temporary, while we are sorting out
     * 3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
     * scale stimulus depending on space_step of elements
     *
     * \todo It looks like the value of the stimulus is specific to 3D
     *
     * @param numEleAcross  Number of elements across which to apply the stimulus
     * @param meshWidth  Width of the mesh (used to calculate magnitude of stimulus)
     * @param useMeshWidthAsMag  see todo comments above (defaults to false).
     * @param stimulusMagnitude  Magnitude of the applied stimulus (defaults to -1e7, modified in the constructor dependent on mesh size).
     * @param stimulusDuration  Duration of the applied stimulus (defaults to 0.5ms).
     */
    GeneralPlaneStimulusCellFactory(unsigned numEleAcross, double meshWidth, bool useMeshWidthAsMag=false, double stimulusMagnitude=-1e7, double stimulusDuration=0.5)
        : PlaneStimulusCellFactory<CELL,ELEMENT_DIM, SPACE_DIM>(stimulusMagnitude,stimulusDuration) // These values are overridden below anyway.
    {
        if (useMeshWidthAsMag)
        {
            #define COVERAGE_IGNORE
            this->mpStimulus.reset(new SimpleStimulus(meshWidth, 0.5));
            #undef COVERAGE_IGNORE
        }
        else
        {
            stimulusMagnitude*=numEleAcross/(64.0);
            // ELEMENT_DIM==1 Justification: elements go half size with each refinement
            // ELEMENT_DIM==2 Justification: Triangles go quarter size with each refinement, but there are twice as many nodes on boundary
            // ELEMENT_DIM==3 Hypothesis: Triangles go eighth size with each refinement, but there are four-times as many nodes on boundary
            stimulusMagnitude*=meshWidth/(0.2);

            //std::cout<<"Mag is "<<stimulusMagnitude<<"\n";
            this->mpStimulus.reset(new SimpleStimulus(stimulusMagnitude, stimulusDuration));
            LOG(1, "Defined a GeneralPlaneStimulusCellFactory<"<<SPACE_DIM<<"> with SimpleStimulus("<<stimulusMagnitude<<","<< stimulusDuration<< ")\n");
        }
    }

};

#endif /*GENERALPLANESTIMULUSCELLFACTORY_HPP_*/
