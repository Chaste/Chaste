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


#include "ContinuumMechanicsProblemDefinition.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"


template<unsigned DIM>
const double ContinuumMechanicsProblemDefinition<DIM>::FREE = DBL_MAX;

template<unsigned DIM>
ContinuumMechanicsProblemDefinition<DIM>::ContinuumMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh)
    : mrMesh(rMesh),
      mDensity(1.0),
      mBodyForceType(CONSTANT_BODY_FORCE),
      mConstantBodyForce(zero_vector<double>(DIM)),
      mTractionBoundaryConditionType(NO_TRACTIONS)
{
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetDensity(double density)
{
    assert(density>0.0);
    mDensity = density;
}

template<unsigned DIM>
double ContinuumMechanicsProblemDefinition<DIM>::GetDensity()
{
    return mDensity;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetBodyForce(c_vector<double,DIM> bodyForce)
{
    mBodyForceType = CONSTANT_BODY_FORCE;
    mConstantBodyForce = bodyForce;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t))
{
    mBodyForceType = FUNCTIONAL_BODY_FORCE;
    mpBodyForceFunction = pFunction;
}


template<unsigned DIM>
BodyForceType ContinuumMechanicsProblemDefinition<DIM>::GetBodyForceType()
{
    return mBodyForceType;
}

template<unsigned DIM>
c_vector<double,DIM> ContinuumMechanicsProblemDefinition<DIM>::GetConstantBodyForce()
{
    assert(mBodyForceType==CONSTANT_BODY_FORCE);
    return mConstantBodyForce;
}

template<unsigned DIM>
c_vector<double,DIM> ContinuumMechanicsProblemDefinition<DIM>::EvaluateBodyForceFunction(c_vector<double,DIM>& rX, double t)
{
    assert(mBodyForceType==FUNCTIONAL_BODY_FORCE);
    return (*mpBodyForceFunction)(rX,t);
}

template<unsigned DIM>
c_vector<double,DIM> ContinuumMechanicsProblemDefinition<DIM>::GetBodyForce(c_vector<double,DIM>& rX, double t)
{
    switch(mBodyForceType)
    {
        case CONSTANT_BODY_FORCE:
        {
            return mConstantBodyForce;
        }
        case FUNCTIONAL_BODY_FORCE:
        {
            return (*mpBodyForceFunction)(rX,t);
        }
        default:
            NEVER_REACHED;
    }
}


template<unsigned DIM>
TractionBoundaryConditionType ContinuumMechanicsProblemDefinition<DIM>::GetTractionBoundaryConditionType()
{
    return mTractionBoundaryConditionType;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                                             std::vector<c_vector<double,DIM> >& rElementwiseTractions)
{

    assert(rTractionBoundaryElements.size()==rElementwiseTractions.size());
    mTractionBoundaryConditionType = ELEMENTWISE_TRACTION;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mElementwiseTractions = rElementwiseTractions;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                                             c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t))
{
    mTractionBoundaryConditionType=FUNCTIONAL_TRACTION;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mpTractionBoundaryConditionFunction = pFunction;
}


template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                                                       double normalPressure)
{
    mTractionBoundaryConditionType = PRESSURE_ON_DEFORMED;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mNormalPressure = normalPressure;

}



template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetZeroDirichletNodes(std::vector<unsigned>& rZeroDirichletNodes)
{
    mDirichletNodes = rZeroDirichletNodes;

    for (unsigned i=0; i<mDirichletNodes.size(); i++)
    {
        assert(mDirichletNodes[i] < mrMesh.GetNumNodes());
    }

    mDirichletNodeValues.clear();
    for (unsigned i=0; i<mDirichletNodes.size(); i++)
    {
        mDirichletNodeValues.push_back(zero_vector<double>(DIM));
    }
}

template<unsigned DIM>
std::vector<unsigned>& ContinuumMechanicsProblemDefinition<DIM>::rGetDirichletNodes()
{
    return mDirichletNodes;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& ContinuumMechanicsProblemDefinition<DIM>::rGetDirichletNodeValues()
{
    return mDirichletNodeValues;
}

template<unsigned DIM>
std::vector<BoundaryElement<DIM-1,DIM>*>& ContinuumMechanicsProblemDefinition<DIM>::rGetTractionBoundaryElements()
{
    return mTractionBoundaryElements;
}


template<unsigned DIM>
std::vector<c_vector<double,DIM> >& ContinuumMechanicsProblemDefinition<DIM>::rGetElementwiseTractions()
{
    assert(mTractionBoundaryConditionType==ELEMENTWISE_TRACTION);
    return mElementwiseTractions;
}


template<unsigned DIM>
double ContinuumMechanicsProblemDefinition<DIM>::GetNormalPressure()
{
    assert(mTractionBoundaryConditionType==PRESSURE_ON_DEFORMED);
    return mNormalPressure;
}

template<unsigned DIM>
c_vector<double,DIM> ContinuumMechanicsProblemDefinition<DIM>::EvaluateTractionFunction(c_vector<double,DIM>& rX, double t)
{
    assert(mTractionBoundaryConditionType==FUNCTIONAL_TRACTION);
    return (*mpTractionBoundaryConditionFunction)(rX,t);
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::Validate()
{
    if(mDirichletNodes.size()==0)
    {
        EXCEPTION("No Dirichlet boundary conditions (eg fixed displacement or fixed flow) have been set");
    }
}



//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class ContinuumMechanicsProblemDefinition<2>;
template class ContinuumMechanicsProblemDefinition<3>;

