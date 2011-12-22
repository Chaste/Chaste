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

#include "AbstractStimulusFactory.hpp"
#include "ZeroStimulus.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<AbstractStimulusFunction>  AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::CreateStimulusForNode(unsigned nodeIndex)
{
    //this is the default implementation
    boost::shared_ptr<ZeroStimulus> p_stim (  new ZeroStimulus() );
    return p_stim;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::GetNumberOfCells()
{
    assert(mpMesh != NULL);
    return mpMesh->GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::AbstractStimulusFactory()
    : mpMesh(NULL)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::~AbstractStimulusFactory()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::SetCompatibleExtracellularStimulus()
{
    //empty in the default case (zero stimulus -> no need to do anything)
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::GetMesh()
{
    if (mpMesh == NULL)
    {
        EXCEPTION("The mesh object has not been set in the stimulus factory");
    }
    return mpMesh;
}

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//void AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::SetRegionToBeGrounded(AbstractChasteRegion<SPACE_DIM>* pRegion)
//{
//    mGroundedRegion = pRegion;
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<AbstractChasteRegion<SPACE_DIM>* > AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::GetRegionsToBeGrounded()
{
    return mGroundedRegions;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractStimulusFactory<1,1>;
template class AbstractStimulusFactory<2,2>;
template class AbstractStimulusFactory<3,3>;
