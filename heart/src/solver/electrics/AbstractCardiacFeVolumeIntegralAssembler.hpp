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


#ifndef ABSTRACTCARDIACFEVOLUMEINTEGRALASSEMBLER_HPP_
#define ABSTRACTCARDIACFEVOLUMEINTEGRALASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "HeartConfig.hpp"
#include "AbstractCardiacTissue.hpp"

/**
 *  Simple implementation of AbstractFeVolumeIntegralAssembler which provides access to a cardiac tissue
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractCardiacFeVolumeIntegralAssembler
   : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>
{
protected:
    /** The Cardiac tissue on which to solve. */
    AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* mpCardiacTissue;

    /** Local cache of the configuration singleton pointer*/
    HeartConfig* mpConfig;

public:

    /**
     *  Constructor
     *  @param pMesh the mesh
     *  @param pTissue  pointer to the tissue used for getting conductivity values
     */
    AbstractCardiacFeVolumeIntegralAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                             AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX,INTERPOLATION_LEVEL>(pMesh),
          mpCardiacTissue(pTissue),
          mpConfig(HeartConfig::Instance())
    {
        assert(pTissue);
    }
};

#endif /*ABSTRACTCARDIACFEVOLUMEINTEGRALASSEMBLER_HPP_*/
