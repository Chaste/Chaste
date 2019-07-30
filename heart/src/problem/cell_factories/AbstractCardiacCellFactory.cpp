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

#include "AbstractCardiacCellFactory.hpp"
#include "FakeBathCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "HeartConfig.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellInterface*  AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::CreateCardiacCellForNode(
    Node<SPACE_DIM>* pNode)
{
    if (HeartRegionCode::IsRegionBath( pNode->GetRegion() ))
    {
        return new FakeBathCell(this->mpSolver, this->mpZeroStimulus);
    }
    else
    {
        AbstractCardiacCellInterface* p_cell = CreateCardiacCellForTissueNode(pNode);
#ifdef CHASTE_CVODE
        if (dynamic_cast<AbstractCvodeCell*>(p_cell))
        {
#if CHASTE_SUNDIALS_VERSION >= 20400
            // Tell the single cell model solver not to re-initialise on subsequent Solve calls.
            // Unfortunately the oldest CVODE we currently support (2.5.0 in Sundials 2.3.0) doesn't
            // work very well this this enabled, so we will go without and take a performance hit.
            static_cast<AbstractCvodeCell*>(p_cell)->SetMinimalReset(true);
#endif // SUNDIALS_VERSION
            // Use the PDE timestep as the [maximum] CVODE timestep.
            static_cast<AbstractCvodeCell*>(p_cell)->SetTimestep(HeartConfig::Instance()->GetPdeTimeStep());
        }
#endif // CHASTE_CVODE
        return p_cell;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::FinaliseCellCreation(
    std::vector< AbstractCardiacCellInterface* >* pCellsDistributed,
    unsigned lo,
    unsigned hi)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::FillInCellularTransmuralAreas()
{
    EXCEPTION("To get here you have probably asked for Epi/Mid/Endo CellularHeterogeneities in your HeartConfig "
              "options or configuration .xml file, to use this you will need to provide a method"
              " `FillInCellularTransmuralAreas()` in your cell factory to override this one.");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::GetNumberOfCells()
{
    assert(mpMesh != NULL);
    return mpMesh->GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacCellFactory(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver)
    : mpMesh(NULL),
      mpHeartGeometryInformation(NULL),
      mpZeroStimulus(new ZeroStimulus),
      mpSolver(pSolver)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::~AbstractCardiacCellFactory()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::GetMesh()
{
    if (mpMesh == NULL)
    {
        EXCEPTION("The mesh object has not been set in the cell factory");
    }
    return mpMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetHeartGeometryInformation(HeartGeometryInformation<SPACE_DIM>* pHeartGeometryInformation)
{
    mpHeartGeometryInformation = pHeartGeometryInformation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>* AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::GetHeartGeometryInformation()
{
    if (mpHeartGeometryInformation == NULL)
    {
        EXCEPTION("HeartGeometryInformation object has not been set in the cell factory");
    }
    return mpHeartGeometryInformation;
}

// Explicit instantiation
template class AbstractCardiacCellFactory<1,1>;
template class AbstractCardiacCellFactory<1,2>;
template class AbstractCardiacCellFactory<1,3>;
template class AbstractCardiacCellFactory<2,2>;
template class AbstractCardiacCellFactory<3,3>;
