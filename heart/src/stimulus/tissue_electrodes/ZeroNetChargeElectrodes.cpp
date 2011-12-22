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

#include "ZeroNetChargeElectrodes.hpp"
#include "ElectrodesStimulusFactory.hpp"
#include "RegularStimulusZeroNetCharge.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "IsNan.hpp"
#include "HeartConfig.hpp"
#include "GaussianQuadratureRule.hpp"

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
boost::shared_ptr<AbstractStimulusFunction> ZeroNetChargeElectrodes<DIM>::CreateStimulusForNode(unsigned nodeIndex)
{
    boost::shared_ptr<RegularStimulusZeroNetCharge> p_stimulus;
    for (unsigned pair_index = 0; pair_index < this->mrElectrodePairs.size(); pair_index++)
    {
        if (this->mrElectrodePairs[pair_index].first->DoesContain(this->mpMesh->GetNode(nodeIndex)->GetPoint()) )
        {
            p_stimulus.reset ( new RegularStimulusZeroNetCharge(this->mMagnitudesElectrode1[pair_index], this->mrDurations[pair_index], this->mrPeriods[pair_index], this->mrStarts[pair_index], this->mrEnds[pair_index]));
        }
        else if (this->mrElectrodePairs[pair_index].second->DoesContain(this->mpMesh->GetNode(nodeIndex)->GetPoint()) )
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


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ZeroNetChargeElectrodes<1>;
template class ZeroNetChargeElectrodes<2>;
template class ZeroNetChargeElectrodes<3>;

