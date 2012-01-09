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

#include "MonodomainTissue.hpp"

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
MonodomainTissue<ELEMENT_DIM,SPACE_DIM>::MonodomainTissue(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory,
            bool exchangeHalos)
    :  AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>(pCellFactory, exchangeHalos)
{
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
MonodomainTissue<ELEMENT_DIM,SPACE_DIM>::MonodomainTissue(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        :  AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>(pMesh)
{
}



/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainTissue<1,1>;
template class MonodomainTissue<1,2>;
template class MonodomainTissue<1,3>;
template class MonodomainTissue<2,2>;
template class MonodomainTissue<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 1)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 3)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 3, 3)
