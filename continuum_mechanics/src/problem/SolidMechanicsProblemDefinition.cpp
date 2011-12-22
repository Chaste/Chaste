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


#include "SolidMechanicsProblemDefinition.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"


template<unsigned DIM>
SolidMechanicsProblemDefinition<DIM>::SolidMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh)
    : ContinuumMechanicsProblemDefinition<DIM>(rMesh)
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
        for(unsigned j=0; j<DIM; j++)
        {
            double location = rFixedNodeLocations[i](j);

            // compute the displacement, assuming the node
            // is not free in this direction
            if(location != this->FREE)
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

    if(compressibilityType==INCOMPRESSIBLE)
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

    assert(this->mrMesh.GetNumElements()==rMaterialLaws.size());

    if(compressibilityType==INCOMPRESSIBLE)
    {
        for(unsigned i=0; i<rMaterialLaws.size(); i++)
        {
            assert(rMaterialLaws[i]);
            AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(rMaterialLaws[i]);
            CheckCastSuccess(compressibilityType, p_law);
            mIncompressibleMaterialLaws.push_back(p_law);
        }
    }
    else
    {
        for(unsigned i=0; i<rMaterialLaws.size(); i++)
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
    // if this fails, SetMaterialLaw() hasn't been called
    assert(mIncompressibleMaterialLaws.size()!=0  ||  mCompressibleMaterialLaws.size()!=0 );
    return mIsHomogeneousMaterial;
}

template<unsigned DIM>
CompressibilityType SolidMechanicsProblemDefinition<DIM>::GetCompressibilityType()
{
    // if this fails, SetMaterialLaw() hasn't been called
    assert(mIncompressibleMaterialLaws.size()!=0  ||  mCompressibleMaterialLaws.size()!=0 );
    return mCompressibilityType;
}




template<unsigned DIM>
AbstractIncompressibleMaterialLaw<DIM>* SolidMechanicsProblemDefinition<DIM>::GetIncompressibleMaterialLaw(unsigned elementIndex)
{
    assert(mCompressibilityType==INCOMPRESSIBLE);
    assert(mIncompressibleMaterialLaws.size()>0);
    assert(mCompressibleMaterialLaws.size()==0);

    if(mIsHomogeneousMaterial)
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
    assert(mCompressibilityType==COMPRESSIBLE);
    assert(mIncompressibleMaterialLaws.size()==0);
    assert(mCompressibleMaterialLaws.size()>0);

    if(mIsHomogeneousMaterial)
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
    if(compressibilityType==INCOMPRESSIBLE && pMaterialLaw==NULL)
    {
        // then dynamic_cast to AbstractIncompressibleMaterialLaw failed
        EXCEPTION("Compressibility type was declared as INCOMPRESSIBLE but a compressible material law was given");
    }

    if(compressibilityType==COMPRESSIBLE && pMaterialLaw==NULL)
    {
        // then dynamic_cast to AbstractCompressibleMaterialLaw failed
        EXCEPTION("Incompressibility type was declared as COMPRESSIBLE but an incompressible material law was given");
    }
}


template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::Validate()
{
    ContinuumMechanicsProblemDefinition<DIM>::Validate();

    if((mIncompressibleMaterialLaws.size()==0)  &&  (mCompressibleMaterialLaws.size()==0))
    {
        EXCEPTION("No material law has been set");
    }
}
//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class SolidMechanicsProblemDefinition<2>;
template class SolidMechanicsProblemDefinition<3>;

