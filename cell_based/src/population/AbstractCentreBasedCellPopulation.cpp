/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::AbstractCentreBasedCellPopulation( AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                    std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh, rCells, locationIndices),
      mMeinekeDivisionSeparation(0.3) // educated guess
{
    // If no location indices are specified, associate with nodes from the mesh.
    std::list<CellPtr>::iterator it = this->mCells.begin();
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();

    for (unsigned i=0; it != this->mCells.end(); ++it, ++i, ++node_iter)
    {
        unsigned index = locationIndices.empty() ? node_iter->GetIndex() : locationIndices[i]; // assume that the ordering matches
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(index,*it);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::AbstractCentreBasedCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>(rMesh),
      mMeinekeDivisionSeparation(0.3) // educated guess

{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return GetNodeCorrespondingToCell(pCell)->rGetLocation();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetNodeCorrespondingToCell(CellPtr pCell)
{
    unsigned index = this->GetLocationIndexUsingCell(pCell);
    return this->GetNode(index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCell(CellPtr pNewCell, const c_vector<double,SPACE_DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Create a new node
    Node<SPACE_DIM>* p_new_node = new Node<SPACE_DIM>(this->GetNumNodes(), rCellDivisionVector, false);   // never on boundary
    unsigned new_node_index = this->AddNode(p_new_node); // use copy constructor so it doesn't matter that new_node goes out of scope

    // Update cells vector
    this->mCells.push_back(pNewCell);

    // Update mappings between cells and location indices
    this->SetCellUsingLocationIndex(new_node_index, pNewCell);
    this->mCellLocationMap[pNewCell.get()] = new_node_index;

    return pNewCell;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<CellPtr,CellPtr> AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CreateCellPair(CellPtr pCell1, CellPtr pCell2)
{
    assert(pCell1);
    assert(pCell2);

    std::pair<CellPtr,CellPtr> cell_pair;

    if (pCell1->GetCellId() < pCell2->GetCellId())
    {
        cell_pair.first = pCell1;
        cell_pair.second = pCell2;
    }
    else
    {
        cell_pair.first = pCell2;
        cell_pair.second = pCell1;
    }
    return cell_pair;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsMarkedSpring(const std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    return mMarkedSprings.find(rCellPair) != mMarkedSprings.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::MarkSpring(std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    mMarkedSprings.insert(rCellPair);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::UnmarkSpring(std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    mMarkedSprings.erase(rCellPair);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetNodeCorrespondingToCell(pCell)->IsDeleted();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);
    return this->GetNeighbouringNodeIndices(node_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::UpdateNodeLocations(double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // Get damping constant for node
        double damping_const = this->GetDampingConstant(node_index);

        // Get displacement
        c_vector<double,SPACE_DIM> displacement=dt*this->GetNode(node_index)->rGetAppliedForce()/damping_const;

        // Throws an exception if the cell movement goes beyond mAbsoluteMovementThreshold
        if (norm_2(displacement) > this->mAbsoluteMovementThreshold)
        {
            EXCEPTION("Cells are moving by: " << norm_2(displacement) <<
                    ", which is more than the AbsoluteMovementThreshold: "
                    << this->mAbsoluteMovementThreshold <<
                    ". Use a smaller timestep to avoid this exception.");
        }

        // Get new node location
        c_vector<double, SPACE_DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        // Create ChastePoint for new node location
        ChastePoint<SPACE_DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetDampingConstant(unsigned nodeIndex)
{
    CellPtr p_cell = this->GetCellUsingLocationIndex(nodeIndex);
    if (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>() && !p_cell->HasCellProperty<CellLabel>())
    {
        return this->GetDampingConstantNormal();
    }
    else
    {
        return this->GetDampingConstantMutant();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsGhostNode(unsigned index)
{
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsParticle(unsigned index)
{
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetMeinekeDivisionSeparation()
{
    return mMeinekeDivisionSeparation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetMeinekeDivisionSeparation(double divisionSeparation)
{
    assert(divisionSeparation <= 1.0);
    assert(divisionSeparation >= 0.0);
    mMeinekeDivisionSeparation = divisionSeparation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AcceptCellWritersAcrossPopulation()
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            CellPtr cell_from_node = this->GetCellUsingLocationIndex(node_iter->GetIndex());
            this->AcceptCellWriter(*cell_writer_iter, cell_from_node);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    // Location of parent and daughter cells
    c_vector<double, SPACE_DIM> parent_coords = this->GetLocationOfCellCentre(pParentCell);
    c_vector<double, SPACE_DIM> daughter_coords;

    // Get separation parameter
    double separation = GetMeinekeDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, SPACE_DIM> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */
    switch (SPACE_DIM)
    {
        case 1:
        {
            double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

            random_vector(0) = 0.5*separation*random_direction;
            break;
        }
        case 2:
        {
            double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

            random_vector(0) = 0.5*separation*cos(random_angle);
            random_vector(1) = 0.5*separation*sin(random_angle);
            break;
        }
        case 3:
        {
            /*
             * Note that to pick a random point on the surface of a sphere, it is incorrect
             * to select spherical coordinates from uniform distributions on [0, 2*pi) and
             * [0, pi) respectively, since points picked in this way will be 'bunched' near
             * the poles. See #2230.
             */
            double u = RandomNumberGenerator::Instance()->ranf();
            double v = RandomNumberGenerator::Instance()->ranf();

            double random_azimuth_angle = 2*M_PI*u;
            double random_zenith_angle = std::acos(2*v - 1);

            random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
            random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
            random_vector(2) = 0.5*separation*cos(random_zenith_angle);
            break;
        }
        default:
            // This can't happen
        NEVER_REACHED;
    }

    parent_coords = parent_coords - random_vector;
    daughter_coords = parent_coords + random_vector;

    // Set the parent to use this location
    ChastePoint<SPACE_DIM> parent_coords_point(parent_coords);
    unsigned node_index = this->GetLocationIndexUsingCell(pParentCell);
    this->SetNode(node_index, parent_coords_point);

    return daughter_coords;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MeinekeDivisionSeparation>" << mMeinekeDivisionSeparation << "</MeinekeDivisionSeparation>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCentreBasedCellPopulation<1,1>;
template class AbstractCentreBasedCellPopulation<1,2>;
template class AbstractCentreBasedCellPopulation<2,2>;
template class AbstractCentreBasedCellPopulation<1,3>;
template class AbstractCentreBasedCellPopulation<2,3>;
template class AbstractCentreBasedCellPopulation<3,3>;
