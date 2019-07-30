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

#include "ZeroNetChargeElectrodes.hpp"
#include "ElectrodesStimulusFactory.hpp"
#include "RegularStimulusZeroNetCharge.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "IsNan.hpp"
#include "HeartConfig.hpp"

template<unsigned DIM>
ZeroNetChargeElectrodes<DIM>::ZeroNetChargeElectrodes(std::vector<std::pair<AbstractChasteRegion<DIM>*, AbstractChasteRegion<DIM>*> >& rElectrodePairs,
                                                          std::vector<double>& rStimulusMagnitudes,
                                                          std::vector<double>& rDurations,
                                                          std::vector<double>& rPeriods,
                                                          std::vector<double>& rStarts,
                                                          std::vector<double>& rEnds)
    : ElectrodesStimulusFactory<DIM>(rElectrodePairs, rStimulusMagnitudes, rDurations, rPeriods, rStarts, rEnds)
{
}

template<unsigned DIM>
ZeroNetChargeElectrodes<DIM>::~ZeroNetChargeElectrodes()
{
}

template<unsigned DIM>
boost::shared_ptr<AbstractStimulusFunction> ZeroNetChargeElectrodes<DIM>::CreateStimulusForNode(Node<DIM>* pNode)
{
    boost::shared_ptr<RegularStimulusZeroNetCharge> p_stimulus;
    for (unsigned pair_index = 0; pair_index < this->mrElectrodePairs.size(); pair_index++)
    {
        if (this->mrElectrodePairs[pair_index].first->DoesContain(pNode->GetPoint()) )
        {
            p_stimulus.reset ( new RegularStimulusZeroNetCharge(this->mMagnitudesElectrode1[pair_index], this->mrDurations[pair_index], this->mrPeriods[pair_index], this->mrStarts[pair_index], this->mrEnds[pair_index]));
        }
        else if (this->mrElectrodePairs[pair_index].second->DoesContain(pNode->GetPoint()) )
        {
            p_stimulus.reset ( new RegularStimulusZeroNetCharge(this->mMagnitudesElectrode2[pair_index], this->mrDurations[pair_index], this->mrPeriods[pair_index], this->mrStarts[pair_index], this->mrEnds[pair_index]));
        }
        else//no stimulus here
        {
            p_stimulus.reset ( new RegularStimulusZeroNetCharge(0.0, this->mrDurations[pair_index], this->mrPeriods[pair_index], this->mrStarts[pair_index], this->mrEnds[pair_index]) );
        }
    }
    return p_stimulus;
}

// Explicit instantiation
template class ZeroNetChargeElectrodes<1>;
template class ZeroNetChargeElectrodes<2>;
template class ZeroNetChargeElectrodes<3>;
