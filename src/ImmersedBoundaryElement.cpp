/*

Copyright (c) 2005-2015, University of Oxford.
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
#include "ImmersedBoundaryElement.hpp"
#include "Exception.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElement(unsigned index,
                                                                         const std::vector<Node<SPACE_DIM>*>& rNodes)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
          mpFluidSource(NULL),
          mAverageNodeSpacing(DOUBLE_UNSET)
{
    assert(ELEMENT_DIM == SPACE_DIM);

    // Ensure number of nodes is at least 2
    assert(rNodes.size() > 2);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::~ImmersedBoundaryElement()
{
    // Do not delete fluid source - that is taken care of in the mesh (which owns the sources)
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::SetFluidSource(FluidSource<SPACE_DIM>* fluidSource)
{
    mpFluidSource = fluidSource;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FluidSource<SPACE_DIM>* ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::GetFluidSource(void)
{
    return mpFluidSource;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<Node<SPACE_DIM>*>& ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::rGetCornerNodes(void)
{
    return mCornerNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::GetAverageNodeSpacing(void)
{
    return mAverageNodeSpacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>::SetAverageNodeSpacing(double averageNodeSpacing)
{
    mAverageNodeSpacing = averageNodeSpacing;
}


//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
ImmersedBoundaryElement<1, SPACE_DIM>::ImmersedBoundaryElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<1, SPACE_DIM>(index, rNodes)
{
}

template<unsigned SPACE_DIM>
void ImmersedBoundaryElement<1, SPACE_DIM>::SetFluidSource(FluidSource<SPACE_DIM>* fluidSource)
{
}

template<unsigned SPACE_DIM>
FluidSource<SPACE_DIM>* ImmersedBoundaryElement<1, SPACE_DIM>::GetFluidSource(void)
{
    return NULL;
}

template<unsigned SPACE_DIM>
std::vector<Node<SPACE_DIM>*>& ImmersedBoundaryElement<1, SPACE_DIM>::rGetCornerNodes(void)
{
    std::vector<Node<SPACE_DIM>*> dummy_vec;
    return dummy_vec;
}

template<unsigned SPACE_DIM>
double ImmersedBoundaryElement<1, SPACE_DIM>::GetAverageNodeSpacing(void)
{
    return 0.0;
}

template<unsigned SPACE_DIM>
void ImmersedBoundaryElement<1, SPACE_DIM>::SetAverageNodeSpacing(double averageNodeSpacing)
{
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class ImmersedBoundaryElement<1,1>;
template class ImmersedBoundaryElement<1,2>;
template class ImmersedBoundaryElement<1,3>;
template class ImmersedBoundaryElement<2,2>;
template class ImmersedBoundaryElement<2,3>;
template class ImmersedBoundaryElement<3,3>;
