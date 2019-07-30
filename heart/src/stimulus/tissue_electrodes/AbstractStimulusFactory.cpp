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

#include "AbstractStimulusFactory.hpp"
#include "ZeroStimulus.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<AbstractStimulusFunction>  AbstractStimulusFactory<ELEMENT_DIM,SPACE_DIM>::CreateStimulusForNode(Node<SPACE_DIM>* pNode)
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

// Explicit instantiation
template class AbstractStimulusFactory<1,1>;
template class AbstractStimulusFactory<2,2>;
template class AbstractStimulusFactory<3,3>;
