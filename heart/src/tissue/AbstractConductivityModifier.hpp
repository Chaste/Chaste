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

#ifndef ABSTRACTCONDUCTIVITYMODIFIER_HPP_
#define ABSTRACTCONDUCTIVITYMODIFIER_HPP_

#include "UblasCustomFunctions.hpp"

/**
 *  Abstract class which just defines an interface and caching method. The pure method
 *  rCalculateModifiedConductivityTensor() should take in a conductivity and return a modified
 *  conductivity (with some dependence e.g. on tissue deformation in cardiac electromechanics).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractConductivityModifier
{
private:
    /** Cache recently-seen elementindex-tensor pairs, one per "domain" */
    std::vector<std::pair<unsigned, c_matrix<double,SPACE_DIM,SPACE_DIM> > > mCache;

public:
    /** Destructor */
    virtual ~AbstractConductivityModifier()
    {
    }

    /** Method that checks element index and the domain (intra/extra/other) and returns cached value if available.
     *  @param elementIndex Index of current element
     *  @param rOriginalConductivity Reference to the original (for example, undeformed) conductivity tensor
     *  @param domainIndex Used to tailor modification to the domain. 0 = intracellular, 1 = extracellular, 2 = second intracellular (tridomain)
     *  @return Reference to a modified conductivity tensor.
     */
    c_matrix<double,SPACE_DIM,SPACE_DIM>& rGetModifiedConductivityTensor(unsigned elementIndex,
                                                                         const c_matrix<double,SPACE_DIM,SPACE_DIM>& rOriginalConductivity,
                                                                         unsigned domainIndex)
    {
        // Have we got space for this domain?
        if (mCache.size() <= domainIndex)
        {
            // Not a pretty line! Initialises every new entry with an UNSIGNED_UNSET and a zero matrix.
            mCache.resize(domainIndex+1, std::pair<unsigned, c_matrix<double,SPACE_DIM,SPACE_DIM> >
                                             (UNSIGNED_UNSET, zero_matrix<double>(SPACE_DIM,SPACE_DIM)));
        }
        // Is this not the same element as last time?
        if (mCache[domainIndex].first != elementIndex)
        {
            mCache[domainIndex].first = elementIndex;
            mCache[domainIndex].second = rCalculateModifiedConductivityTensor(elementIndex, rOriginalConductivity, domainIndex);
        }
        // Return cached tensor
        return mCache[domainIndex].second;
    }

    /**
     * Pure method that alters the given conductivity tensor on an element-wise basis for a particular domain.
     *
     * @param elementIndex Global index of current element.
     * @param rOriginalConductivity Reference to the original (for example, undeformed) conductivity tensor.
     * @param domainIndex  The index of the domain (0=intracellular, 1=extracellular, [2=second intracellular in extended bidomain]).
     * @return Reference to a modified conductivity tensor,
     *         N.B. the fact this is a reference means the tensor object has to persist,
     *         and should therefore generally be a member variable of your subclass
     *         (it gets copied/cached appropriately by the calling code, so doesn't have to persist for subsequent calls!).
     */
    virtual c_matrix<double,SPACE_DIM,SPACE_DIM>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                                                       const c_matrix<double,SPACE_DIM,SPACE_DIM>& rOriginalConductivity,
                                                                                       unsigned domainIndex)=0;
};


#endif /* ABSTRACTCONDUCTIVITYMODIFIER_HPP_ */
