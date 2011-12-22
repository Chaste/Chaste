/*

Copyright (C) University of Oxford, 2005-2011

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

#include "DynamicRestitutionStimulus.hpp"

DynamicRestitutionStimulus::DynamicRestitutionStimulus(double magnitude, double stimulusDuration, double startTime, std::vector<double> pacingCycleLengths, unsigned numberOfPulses)
    : MultiStimulus()
{
    for (unsigned stim=0; stim < pacingCycleLengths.size() ; stim++)
    {
        boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(magnitude, stimulusDuration, pacingCycleLengths[stim], startTime, startTime + numberOfPulses*pacingCycleLengths[stim]));
        this->AddStimulus(p_stim);
        // increment start time to the end of this pulse train.
        startTime = startTime + numberOfPulses*pacingCycleLengths[stim];
    }
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DynamicRestitutionStimulus)


