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

#ifndef ZERONETCHARGEELECTRODES_HPP_
#define ZERONETCHARGEELECTRODES_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedVector.hpp"
#include "ElectrodesStimulusFactory.hpp"
#include "AbstractChasteRegion.hpp"
#include "AbstractStimulusFactory.hpp"
#include "LinearBasisFunction.hpp"

/**
 * This class implements the specification of two electrodes with a RegularStimulusZeroNetCharge applied to them.
 * Note that the compatibility conditions are the same as its parent class (ElectrodesStimulusFactory), becuase the
 * magnitudes involved are either + or minus 'Magnitude', hence the stimulus will be compatible.
 *
 * See documentation of ElectrodesStimulusFactory for full functionality.
 *
 */
template<unsigned DIM>
class ZeroNetChargeElectrodes : public ElectrodesStimulusFactory<DIM>
{

public:

    /**
     * Constructor. Electrodes and stimulation parameters need to be passed in.
     *
     * @param rElectrodePairs the pairs of electrodes
     * @param rStimulusMagnitudes the magnitudes of the stimuli (microA / cm^3). First electrode will have magnitude value and second electrode will have -magnitude (before being corrected to ensure equal flux).
     * @param rDurations the duration of each stimulus (ms)
     * @param rPeriods the period of each stimulus (ms)
     * @param rStarts the start time of each stimulus (ms).
     * @param rEnds the end of each stimulation (ms)
     */
    ZeroNetChargeElectrodes(std::vector<std::pair<AbstractChasteRegion<DIM>*, AbstractChasteRegion<DIM>*> >& rElectrodePairs,
                              std::vector<double>& rStimulusMagnitudes,
                              std::vector<double>& rDurations,
                              std::vector<double>& rPeriods,
                              std::vector<double>& rStarts,
                              std::vector<double>& rEnds);

    /**
     * Creates an appropriate stimuus for each node as to abide compatibility conditions.
     *
     * @param nodeIndex the node index for which to create the stimulus
     */
    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(unsigned nodeIndex);

    /**
     * Destructor
     */
    ~ZeroNetChargeElectrodes();

};

#endif /*ZERONETCHARGEELECTRODES_HPP_*/
