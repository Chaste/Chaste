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

#include <limits>

#include "ContinuumMechanicsProblemDefinition.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"


template<unsigned DIM>
const double ContinuumMechanicsProblemDefinition<DIM>::FREE = std::numeric_limits<double>::max();

template<unsigned DIM>
ContinuumMechanicsProblemDefinition<DIM>::ContinuumMechanicsProblemDefinition(AbstractTetrahedralMesh<DIM,DIM>& rMesh)
    : mrMesh(rMesh),
      mDensity(1.0),
      mBodyForceType(CONSTANT_BODY_FORCE),
      mConstantBodyForce(zero_vector<double>(DIM)),
      mTractionBoundaryConditionType(NO_TRACTIONS),
      mVerboseDuringSolve(false)
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
void ContinuumMechanicsProblemDefinition<DIM>::SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                                             c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t))
{
    mTractionBoundaryConditionType=FUNCTIONAL_TRACTION;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mpTractionBoundaryConditionFunction = pFunction;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                                                       double normalPressure)
{
    mTractionBoundaryConditionType = PRESSURE_ON_DEFORMED;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mNormalPressure = normalPressure;
    mOriginalNormalPressure = normalPressure;

}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                                                       double (*pFunction)(double t))
{
    mTractionBoundaryConditionType = FUNCTIONAL_PRESSURE_ON_DEFORMED;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mpNormalPressureFunction = pFunction;
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
double ContinuumMechanicsProblemDefinition<DIM>::EvaluateNormalPressureFunction(double t)
{
    assert(mTractionBoundaryConditionType==FUNCTIONAL_PRESSURE_ON_DEFORMED);
    return (*mpNormalPressureFunction)(t);
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::SetPressureScaling(double scaleFactor)
{
    assert(mTractionBoundaryConditionType==PRESSURE_ON_DEFORMED);
    mNormalPressure = mOriginalNormalPressure*scaleFactor;
}

template<unsigned DIM>
void ContinuumMechanicsProblemDefinition<DIM>::Validate()
{
    if (mDirichletNodes.size() == 0)
    {
        EXCEPTION("No Dirichlet boundary conditions (eg fixed displacement or fixed flow) have been set");
    }
}

// Explicit instantiation
template class ContinuumMechanicsProblemDefinition<2>;
template class ContinuumMechanicsProblemDefinition<3>;
