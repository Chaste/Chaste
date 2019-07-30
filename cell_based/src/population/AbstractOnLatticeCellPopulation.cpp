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

#include "AbstractOnLatticeCellPopulation.hpp"

template<unsigned DIM>
AbstractOnLatticeCellPopulation<DIM>::AbstractOnLatticeCellPopulation(AbstractMesh<DIM, DIM>& rMesh,
                                                                    std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices,
                                                                  bool deleteMesh)
    : AbstractCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh),
      mUpdateNodesInRandomOrder(true),
      mIterateRandomlyOverUpdateRuleCollection(false)
{
    std::list<CellPtr>::iterator it = this->mCells.begin();
    for (unsigned i=0; it != this->mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        AbstractCellPopulation<DIM>::AddCellUsingLocationIndex(index, *it);
    }
}

template<unsigned DIM>
AbstractOnLatticeCellPopulation<DIM>::AbstractOnLatticeCellPopulation(AbstractMesh<DIM, DIM>& rMesh)
    : AbstractCellPopulation<DIM>(rMesh),
      mDeleteMesh(true),
      mUpdateNodesInRandomOrder(true),
      mIterateRandomlyOverUpdateRuleCollection(false)
{
}

template<unsigned DIM>
AbstractOnLatticeCellPopulation<DIM>::~AbstractOnLatticeCellPopulation()
{
}

template<unsigned DIM>
bool AbstractOnLatticeCellPopulation<DIM>::GetUpdateNodesInRandomOrder()
{
    return mUpdateNodesInRandomOrder;
}

template<unsigned DIM>
void AbstractOnLatticeCellPopulation<DIM>::SetUpdateNodesInRandomOrder(bool updateNodesInRandomOrder)
{
    mUpdateNodesInRandomOrder = updateNodesInRandomOrder;
}

template<unsigned DIM>
void AbstractOnLatticeCellPopulation<DIM>::SetIterateRandomlyOverUpdateRuleCollection(bool iterateRandomly)
{
    mIterateRandomlyOverUpdateRuleCollection = iterateRandomly;
}

template<unsigned DIM>
bool AbstractOnLatticeCellPopulation<DIM>::GetIterateRandomlyOverUpdateRuleCollection()
{
    return mIterateRandomlyOverUpdateRuleCollection;
}

template<unsigned DIM>
void AbstractOnLatticeCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    EXCEPTION("SetNode() cannot be called on a subclass of AbstractOnLatticeCellPopulation.");
}

template<unsigned DIM>
std::set<unsigned> AbstractOnLatticeCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    EXCEPTION("Cannot call GetNeighbouringNodeIndices() on a subclass of AbstractOnLatticeCellPopulation, need to go through the PottsMesh instead");
    std::set<unsigned> neighbouring_node_indices;
    return neighbouring_node_indices;
}

template<unsigned DIM>
void AbstractOnLatticeCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UpdateNodesInRandomOrder>" << mUpdateNodesInRandomOrder << "</UpdateNodesInRandomOrder>\n";
    *rParamsFile << "\t\t<IterateRandomlyOverUpdateRuleCollection>" << mIterateRandomlyOverUpdateRuleCollection << "</IterateRandomlyOverUpdateRuleCollection>\n";

    // Call method on direct parent class
    AbstractCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double AbstractOnLatticeCellPopulation<DIM>::GetDefaultTimeStep()
{
    return 0.1;
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > AbstractOnLatticeCellPopulation<DIM>::GetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void AbstractOnLatticeCellPopulation<DIM>::RemoveAllUpdateRules()
{
    mUpdateRuleCollection.clear();
}

// Explicit instantiation
template class AbstractOnLatticeCellPopulation<1>;
template class AbstractOnLatticeCellPopulation<2>;
template class AbstractOnLatticeCellPopulation<3>;
