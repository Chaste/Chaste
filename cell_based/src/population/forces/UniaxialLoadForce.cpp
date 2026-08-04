/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "UniaxialLoadForce.hpp"

#include "ChasteCuboid.hpp"
#include "Exception.hpp"

template<unsigned DIM>
UniaxialLoadForce<DIM>::UniaxialLoadForce()
    : AbstractForce<DIM>(),
      mLoad(0.0),
      mLoadingAxis(DIM - 1u),
      mSlabThickness(0.0),
      mIsCompressive(true),
      mSlabsIdentified(false)
{
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::SetLoad(double load)
{
    if (load <= 0.0)
    {
        EXCEPTION("UniaxialLoadForce: load must be positive");
    }
    mLoad = load;
}

template<unsigned DIM>
double UniaxialLoadForce<DIM>::GetLoad() const
{
    return mLoad;
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::SetLoadingAxis(unsigned loadingAxis)
{
    if (loadingAxis >= DIM)
    {
        EXCEPTION("UniaxialLoadForce: loadingAxis must be less than DIM");
    }
    mLoadingAxis = loadingAxis;
}

template<unsigned DIM>
unsigned UniaxialLoadForce<DIM>::GetLoadingAxis() const
{
    return mLoadingAxis;
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::SetSlabThickness(double slabThickness)
{
    if (slabThickness <= 0.0)
    {
        EXCEPTION("UniaxialLoadForce: slabThickness must be positive");
    }
    mSlabThickness = slabThickness;
}

template<unsigned DIM>
double UniaxialLoadForce<DIM>::GetSlabThickness() const
{
    return mSlabThickness;
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::SetIsCompressive(bool isCompressive)
{
    mIsCompressive = isCompressive;
}

template<unsigned DIM>
bool UniaxialLoadForce<DIM>::GetIsCompressive() const
{
    return mIsCompressive;
}

template<unsigned DIM>
unsigned UniaxialLoadForce<DIM>::GetUpperSlabNodeCount() const
{
    return mUpperSlabNodeIndices.size();
}

template<unsigned DIM>
unsigned UniaxialLoadForce<DIM>::GetLowerSlabNodeCount() const
{
    return mLowerSlabNodeIndices.size();
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::IdentifySlabs(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Find the extent of the specimen along the loading axis
    const ChasteCuboid<DIM> bounding_box = rCellPopulation.rGetMesh().CalculateBoundingBox();
    const double lower_end = bounding_box.rGetLowerCorner()[mLoadingAxis];
    const double upper_end = bounding_box.rGetUpperCorner()[mLoadingAxis];

    // Slabs thicker than half the specimen would share nodes, so that the two loads would cancel
    // node by node rather than deforming anything
    if (2.0 * mSlabThickness > bounding_box.GetWidth(mLoadingAxis))
    {
        EXCEPTION("UniaxialLoadForce: slabThickness is more than half the extent of the "
                  "population along the loading axis, so the two slabs would overlap");
    }

    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        const double coordinate = node_iter->rGetLocation()[mLoadingAxis];

        if (coordinate > upper_end - mSlabThickness)
        {
            mUpperSlabNodeIndices.push_back(node_iter->GetIndex());
        }
        else if (coordinate < lower_end + mSlabThickness)
        {
            mLowerSlabNodeIndices.push_back(node_iter->GetIndex());
        }
    }

    mSlabsIdentified = true;
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (mLoad <= 0.0)
    {
        EXCEPTION("UniaxialLoadForce: SetLoad() must be called before the force is used");
    }

    if (mSlabThickness <= 0.0)
    {
        EXCEPTION("UniaxialLoadForce: SetSlabThickness() must be called before the force is used");
    }

    if (!mSlabsIdentified)
    {
        IdentifySlabs(rCellPopulation);
    }

    // Each slab carries a total load of mLoad, so dividing by that slab's own node count leaves
    // the two loads exactly equal and opposite however many nodes each contains
    const double direction = mIsCompressive ? -1.0 : 1.0;
    const double upper_force = direction * mLoad / static_cast<double>(mUpperSlabNodeIndices.size());
    const double lower_force = -direction * mLoad / static_cast<double>(mLowerSlabNodeIndices.size());

    c_vector<double, DIM> force_contribution = zero_vector<double>(DIM);

    force_contribution[mLoadingAxis] = upper_force;
    for (unsigned node_index : mUpperSlabNodeIndices)
    {
        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force_contribution);
    }

    force_contribution[mLoadingAxis] = lower_force;
    for (unsigned node_index : mLowerSlabNodeIndices)
    {
        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void UniaxialLoadForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Load>" << mLoad << "</Load>\n";
    *rParamsFile << "\t\t\t<LoadingAxis>" << mLoadingAxis << "</LoadingAxis>\n";
    *rParamsFile << "\t\t\t<SlabThickness>" << mSlabThickness << "</SlabThickness>\n";
    *rParamsFile << "\t\t\t<IsCompressive>" << mIsCompressive << "</IsCompressive>\n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class UniaxialLoadForce<1>;
template class UniaxialLoadForce<2>;
template class UniaxialLoadForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(UniaxialLoadForce)
