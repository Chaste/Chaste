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

#include "AveragedSourceParabolicPde.hpp"
#include "ApoptoticCellProperty.hpp"

template<unsigned DIM>
AveragedSourceParabolicPde<DIM>::AveragedSourceParabolicPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double constantCellSourceCoefficient, 
                                                            double linearCellSourceCoefficient, 
                                                            double constantSourceCoefficient, 
                                                            double linearSourceCoefficient, 
                                                            double diffusionCoefficient,
                                                            double duDtCoefficient,
                                                            bool scaleByCellVolume)
    : mrCellPopulation(rCellPopulation),
      mConstantCellSourceCoefficient(constantCellSourceCoefficient),
      mLinearCellSourceCoefficient(linearCellSourceCoefficient),
      mConstantSourceCoefficient(constantSourceCoefficient),
      mLinearSourceCoefficient(linearSourceCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mDuDtCoefficient(duDtCoefficient),
      mScaleByCellVolume(scaleByCellVolume)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& AveragedSourceParabolicPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetConstantCellCoefficient() const
{
    return mConstantCellSourceCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetLinearCellCoefficient() const
{
    return mLinearCellSourceCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetConstantCoefficient() const
{
    return mConstantSourceCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetLinearCoefficient() const
{
    return mLinearSourceCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetDiffusionCoefficient() const
{
    return mDiffusionCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetDuDtCoefficient() const
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
bool AveragedSourceParabolicPde<DIM>::GetScaleByCellVolume() const
{
    return mScaleByCellVolume;
}

template<unsigned DIM>
void AveragedSourceParabolicPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (double & mCellDensityOnCoarseElement : mCellDensityOnCoarseElements)
    {
        mCellDensityOnCoarseElement = 0.0;
    }

    // bool mSmallMesh = false;

    // if(mSmallMesh)
    // {
    //     std::vector<unsigned> nearest_cell;
    //     std::vector< std::set<unsigned> > voronoi_cell_node_indices;

    //     // Allocate memory
    //     nearest_cell.resize(rCoarseMesh.GetNumElements());
    //     for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    //     {
    //         nearest_cell[elem_index] = UNSIGNED_UNSET;
    //     }

    //     // Double check set cleared
    //     voronoi_cell_node_indices.resize(mrCellPopulation.GetNumAllCells());
    //     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
    //         cell_iter != mrCellPopulation.End();
    //         ++cell_iter)
    //     {
    //         unsigned cell_location_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);   
    //         voronoi_cell_node_indices[cell_location_index].clear(); 
    //     }

    //     // Loop over elements of the finite element mesh and work out which voronoi region the centre of the element is in.
    //     for (typename TetrahedralMesh<DIM,DIM>::ElementIterator elem_iter = rCoarseMesh.GetElementIteratorBegin();
    //      elem_iter != rCoarseMesh.GetElementIteratorEnd();
    //      ++elem_iter)
    //     {
    //         unsigned element_index = elem_iter->GetIndex();

    //         c_vector<double,DIM> element_centre_location = elem_iter->CalculateCentroid();

    //         double closest_separation = DBL_MAX;
    //         unsigned local_nearest_cell = UNSIGNED_UNSET;

    //         for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
    //             cell_iter != mrCellPopulation.End();
    //             ++cell_iter)
    //         {
    //             c_vector<double, DIM> cell_location = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //             double separation = norm_2(element_centre_location - cell_location);

    //             if (separation < closest_separation)
    //             {
    //                 closest_separation = separation;
    //                 local_nearest_cell = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //             }                
    //         }
    //         assert(closest_separation<DBL_MAX);

    //         nearest_cell[element_index] = local_nearest_cell;

    //         voronoi_cell_node_indices[local_nearest_cell].insert(element_index);
    //     }   

    //     // Loop over the cells and assign the weights to the FE Mesh elements
    //     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
    //     cell_iter != mrCellPopulation.End();
    //     ++cell_iter)
    //     {
    //         unsigned cell_location_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);   

    //         unsigned number_of_contained_nodes = voronoi_cell_node_indices[cell_location_index].size();
            
    //         if (number_of_contained_nodes==0)
    //         {
    //             //EXCEPTION("One or more of the cells doesnt contain any pde nodes can't use mSmallMesh=true");

    //             /*
    //              * Cell doesnt contain any element centroids, probably as cells smaller than elements. 
    //              * In this case assign all the weight to the element that contains the cell centre
    //              */
    //             unsigned elem_index = UNSIGNED_UNSET;
    //             const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //             if (pCellPdeElementMap != nullptr)
    //             {
    //                 elem_index = (*pCellPdeElementMap)[*cell_iter];
    //             }
    //             else
    //             {
    //                 elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
    //             }

    //             voronoi_cell_node_indices[cell_location_index].insert(elem_index);
    //             number_of_contained_nodes = 1;
    //         }

    //         double cell_weight = 1.0;

    //         if (mScaleByCellVolume)
    //         {   
    //             cell_weight = mrCellPopulation.GetVolumeOfCell(*cell_iter);

    //             if (cell_weight <1e-6)
    //             {
    //                 EXCEPTION("The volume of one of the cells is " << cell_weight << 
    //                         " and you are scaling by cell volume. Either turn scaling off or use"  
    //                         " a cell model with non zero areas (i.e. a Bounded Voronoi Tesselation model).");
    //             }
    //         }

    //         // Partition the cell between the containing elements
    //         for (std::set<unsigned>::iterator iter =  voronoi_cell_node_indices[cell_location_index].begin();
    //             iter !=  voronoi_cell_node_indices[cell_location_index].end();
    //             ++iter)
    //         {
    //             mCellDensityOnCoarseElements[*iter] += cell_weight/(double) number_of_contained_nodes;
    //         }
            
    //     }
    // }
    // else
    {
        // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
            cell_iter != mrCellPopulation.End();
            ++cell_iter)
        {
            unsigned elem_index = UNSIGNED_UNSET;
            const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

            if (pCellPdeElementMap != nullptr)
            {
                elem_index = (*pCellPdeElementMap)[*cell_iter];
            }
            else
            {
                elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
            }

            bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

            if (!cell_is_apoptotic)
            {
                double cell_weight = 1.0;

                if (mScaleByCellVolume)
                {   
                    // If scaling by cell volume then use volume here instead of cell count 
                    cell_weight = mrCellPopulation.GetVolumeOfCell(*cell_iter);

                    if (cell_weight <1e-6)
                    {
                        EXCEPTION("The volume of one of the cells is " << cell_weight << 
                                " and you are scaling by cell volume. Either turn scaling off or use"  
                                " a cell model with non zero areas (i.e. a Bounded Voronoi Tesselation model).");
                    }
                }

                mCellDensityOnCoarseElements[elem_index] += cell_weight;
            }
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    assert(!mCellDensityOnCoarseElements.empty());

    // The source term is (a*density + c)*u + b*density + d
    double constant_source_term = mConstantCellSourceCoefficient * mCellDensityOnCoarseElements[pElement->GetIndex()] + mConstantSourceCoefficient;
    double linear_source_term_coeficient = mLinearCellSourceCoefficient * mCellDensityOnCoarseElements[pElement->GetIndex()] + mLinearSourceCoefficient;
    
    double source_term = constant_source_term + linear_source_term_coeficient * u;

    return source_term;
    //return 25*mCellDensityOnCoarseElements[pElement->GetIndex()] - u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
c_matrix<double,DIM,DIM> AveragedSourceParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    NEVER_REACHED;
    return this->mCellDensityOnCoarseElements[elementIndex];
}

// Explicit instantiation
template class AveragedSourceParabolicPde<1>;
template class AveragedSourceParabolicPde<2>;
template class AveragedSourceParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AveragedSourceParabolicPde)
