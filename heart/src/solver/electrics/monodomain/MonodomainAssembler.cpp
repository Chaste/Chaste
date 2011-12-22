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

#ifndef MONODOMAINASSEMBLER_CPP_
#define MONODOMAINASSEMBLER_CPP_

#include "MonodomainAssembler.hpp"
#include "PdeSimulationTime.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    /// Am and Cm are set as scaling factors for the mass matrix in its constructor.
    return (PdeSimulationTime::GetPdeTimeStepInverse())*mMassMatrixAssembler.ComputeMatrixTerm(rPhi,rGradPhi,rX,rU,rGradU,pElement)
            + mStiffnessMatrixAssembler.ComputeMatrixTerm(rPhi,rGradPhi,rX,rU,rGradU,pElement);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainAssembler(
                        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                        unsigned numQuadPoints)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,CARDIAC>(pMesh,pTissue,numQuadPoints),
      mMassMatrixAssembler(pMesh, HeartConfig::Instance()->GetUseMassLumping(),
                           HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio()*HeartConfig::Instance()->GetCapacitance()),
      mStiffnessMatrixAssembler(pMesh, pTissue)
{
    assert(pTissue);
    mpConfig = HeartConfig::Instance();
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////


template class MonodomainAssembler<1,1>;
template class MonodomainAssembler<1,2>;
template class MonodomainAssembler<1,3>;
template class MonodomainAssembler<2,2>;
template class MonodomainAssembler<3,3>;

#endif /*MONODOMAINASSEMBLER_CPP_*/
