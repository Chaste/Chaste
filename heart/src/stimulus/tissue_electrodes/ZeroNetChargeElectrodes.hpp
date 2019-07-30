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
    ZeroNetChargeElectrodes(std::vector<std::pair<AbstractChasteRegion<DIM>*,
                            AbstractChasteRegion<DIM>*> >& rElectrodePairs,
                            std::vector<double>& rStimulusMagnitudes,
                            std::vector<double>& rDurations,
                            std::vector<double>& rPeriods,
                            std::vector<double>& rStarts,
                            std::vector<double>& rEnds);

    /**
     * @return Creates an appropriate stimulus for each node as to abide compatibility conditions.
     *
     * @param pNode the node for which to create the stimulus
     */
    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<DIM>* pNode);

    /**
     * Destructor
     */
    ~ZeroNetChargeElectrodes();
};

#endif /*ZERONETCHARGEELECTRODES_HPP_*/
