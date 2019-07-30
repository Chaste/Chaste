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


#include "SolidMechanicsProblemDefinition.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"

template<unsigned DIM>
SolidMechanicsProblemDefinition<DIM>::SolidMechanicsProblemDefinition(AbstractTetrahedralMesh<DIM,DIM>& rMesh)
    : ContinuumMechanicsProblemDefinition<DIM>(rMesh),
      mSolveUsingSnes(false)
{
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetFixedNodes(std::vector<unsigned>& rFixedNodes, std::vector<c_vector<double,DIM> >& rFixedNodeLocations)
{
    assert(rFixedNodes.size()==rFixedNodeLocations.size());
    this->mDirichletNodes = rFixedNodes;

    this->mDirichletNodeValues.clear();
    for (unsigned i=0; i<this->mDirichletNodes.size(); i++)
    {
        unsigned index = this->mDirichletNodes[i];
        c_vector<double,DIM> displacement;
        for (unsigned j=0; j<DIM; j++)
        {
            double location = rFixedNodeLocations[i](j);

            // Compute the displacement, assuming the node is not free in this direction
            if (location != this->FREE)
            {
                displacement(j) = location - this->mrMesh.GetNode(index)->rGetLocation()[j];
            }
            else
            {
                displacement(j) = this->FREE;
            }
        }
        this->mDirichletNodeValues.push_back(displacement);
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetMaterialLaw(CompressibilityType compressibilityType,
                                                          AbstractMaterialLaw<DIM>* pMaterialLaw)
{
    mIsHomogeneousMaterial = true;
    mCompressibilityType = compressibilityType;

    mIncompressibleMaterialLaws.clear();
    mCompressibleMaterialLaws.clear();

    assert(pMaterialLaw);

    if (compressibilityType == INCOMPRESSIBLE)
    {
        AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(pMaterialLaw);
        CheckCastSuccess(compressibilityType, p_law);
        mIncompressibleMaterialLaws.push_back(p_law);
    }
    else
    {
        AbstractCompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractCompressibleMaterialLaw<DIM>*>(pMaterialLaw);
        CheckCastSuccess(compressibilityType, p_law);
        mCompressibleMaterialLaws.push_back(p_law);
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetMaterialLaw(CompressibilityType compressibilityType,
                                                          std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws)
{
    mIsHomogeneousMaterial = false;
    mCompressibilityType = compressibilityType;

    mIncompressibleMaterialLaws.clear();
    mCompressibleMaterialLaws.clear();

    assert(this->mrMesh.GetNumElements() == rMaterialLaws.size());

    if (compressibilityType == INCOMPRESSIBLE)
    {
        for (unsigned i=0; i<rMaterialLaws.size(); i++)
        {
            assert(rMaterialLaws[i]);
            AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(rMaterialLaws[i]);
            CheckCastSuccess(compressibilityType, p_law);
            mIncompressibleMaterialLaws.push_back(p_law);
        }
    }
    else
    {
        for (unsigned i=0; i<rMaterialLaws.size(); i++)
        {
            assert(rMaterialLaws[i]);
            AbstractCompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractCompressibleMaterialLaw<DIM>*>(rMaterialLaws[i]);
            CheckCastSuccess(compressibilityType, p_law);
            mCompressibleMaterialLaws.push_back(p_law);
        }
    }
}

template<unsigned DIM>
bool SolidMechanicsProblemDefinition<DIM>::IsHomogeneousMaterial()
{
    // If this fails, SetMaterialLaw() hasn't been called
    assert(mIncompressibleMaterialLaws.size()!=0  ||  mCompressibleMaterialLaws.size()!=0 );
    return mIsHomogeneousMaterial;
}

template<unsigned DIM>
CompressibilityType SolidMechanicsProblemDefinition<DIM>::GetCompressibilityType()
{
    // If this fails, SetMaterialLaw() hasn't been called
    assert(mIncompressibleMaterialLaws.size()!=0  ||  mCompressibleMaterialLaws.size()!=0 );
    return mCompressibilityType;
}

template<unsigned DIM>
AbstractIncompressibleMaterialLaw<DIM>* SolidMechanicsProblemDefinition<DIM>::GetIncompressibleMaterialLaw(unsigned elementIndex)
{
    assert(mCompressibilityType==INCOMPRESSIBLE);
    assert(mIncompressibleMaterialLaws.size()>0);
    assert(mCompressibleMaterialLaws.size()==0);

    if (mIsHomogeneousMaterial)
    {
        return mIncompressibleMaterialLaws[0];
    }
    else
    {
        assert(elementIndex < this->mrMesh.GetNumNodes());
        return mIncompressibleMaterialLaws[elementIndex];
    }
}

template<unsigned DIM>
AbstractCompressibleMaterialLaw<DIM>* SolidMechanicsProblemDefinition<DIM>::GetCompressibleMaterialLaw(unsigned elementIndex)
{
    assert(mCompressibilityType == COMPRESSIBLE);
    assert(mIncompressibleMaterialLaws.size() == 0);
    assert(mCompressibleMaterialLaws.size() > 0);

    if (mIsHomogeneousMaterial)
    {
        return mCompressibleMaterialLaws[0];
    }
    else
    {
        assert(elementIndex < this->mrMesh.GetNumNodes());
        return mCompressibleMaterialLaws[elementIndex];
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::CheckCastSuccess(CompressibilityType compressibilityType,
                                                            AbstractMaterialLaw<DIM>* pMaterialLaw)
{
    if ((compressibilityType==INCOMPRESSIBLE) && (pMaterialLaw==nullptr))
    {
        // then dynamic_cast to AbstractIncompressibleMaterialLaw failed
        EXCEPTION("Compressibility type was declared as INCOMPRESSIBLE but a compressible material law was given");
    }

    if ((compressibilityType==COMPRESSIBLE) && (pMaterialLaw==nullptr))
    {
        // then dynamic_cast to AbstractCompressibleMaterialLaw failed
        EXCEPTION("Incompressibility type was declared as COMPRESSIBLE but an incompressible material law was given");
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::Validate()
{
    ContinuumMechanicsProblemDefinition<DIM>::Validate();

    if ((mIncompressibleMaterialLaws.size()==0)  &&  (mCompressibleMaterialLaws.size()==0))
    {
        EXCEPTION("No material law has been set");
    }
}

// Explicit instantiation
template class SolidMechanicsProblemDefinition<2>;
template class SolidMechanicsProblemDefinition<3>;
