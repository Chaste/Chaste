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

public:
    /**
     *  Constructor
     *  @param pMesh the mesh
     *  @param pTissue  pointer to the tissue used for getting conductivity values
     *  @param numQuadPoints  The number of quadrature points to use (other constructor takes none and uses default specified in AbstractFeVolumeIntegralAssembler).
     */
    AbstractCardiacFeVolumeIntegralAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                             AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                                             unsigned numQuadPoints)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX,INTERPOLATION_LEVEL>(pMesh, numQuadPoints),
          mpCardiacTissue(pTissue)
    {
        assert(pTissue);
    }

    /**
     *  Constructor
     *  @param pMesh the mesh
     *  @param pTissue  pointer to the tissue used for getting conductivity values
     */
    AbstractCardiacFeVolumeIntegralAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                             AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
        : AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX,INTERPOLATION_LEVEL>(pMesh),
          mpCardiacTissue(pTissue)
    {
        assert(pTissue);
    }
};

#endif /*ABSTRACTCARDIACFEVOLUMEINTEGRALASSEMBLER_HPP_*/
