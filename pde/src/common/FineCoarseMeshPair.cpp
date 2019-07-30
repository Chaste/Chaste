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

#include "FineCoarseMeshPair.hpp"

template<unsigned DIM>
FineCoarseMeshPair<DIM>::FineCoarseMeshPair(AbstractTetrahedralMesh<DIM,DIM>& rFineMesh, AbstractTetrahedralMesh<DIM,DIM>& rCoarseMesh)
    : mrFineMesh(rFineMesh),
      mrCoarseMesh(rCoarseMesh),
      mpFineMeshBoxCollection(nullptr),
      mpCoarseMeshBoxCollection(nullptr)
{
    ResetStatisticsVariables();
}

template<unsigned DIM>
const AbstractTetrahedralMesh<DIM,DIM>& FineCoarseMeshPair<DIM>::GetFineMesh() const
{
    return  mrFineMesh;
}

template<unsigned DIM>
const AbstractTetrahedralMesh<DIM,DIM>& FineCoarseMeshPair<DIM>::GetCoarseMesh() const
{
    return  mrCoarseMesh;
}

template<unsigned DIM>
FineCoarseMeshPair<DIM>::~FineCoarseMeshPair()
{
    DeleteFineBoxCollection();
    DeleteCoarseBoxCollection();
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::DeleteFineBoxCollection()
{
    if (mpFineMeshBoxCollection != nullptr)
    {
        delete mpFineMeshBoxCollection;
        mpFineMeshBoxCollection = nullptr;
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::DeleteCoarseBoxCollection()
{
    if (mpCoarseMeshBoxCollection != nullptr)
    {
        delete mpCoarseMeshBoxCollection;
        mpCoarseMeshBoxCollection = nullptr;
    }
}

////////////////////////////////////////////////////////////////////////////////////
//   Setting up boxes methods
////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::SetUpBoxesOnFineMesh(double boxWidth)
{
    SetUpBoxes(mrFineMesh, boxWidth, mpFineMeshBoxCollection);
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::SetUpBoxesOnCoarseMesh(double boxWidth)
{
    SetUpBoxes(mrCoarseMesh, boxWidth, mpCoarseMeshBoxCollection);
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::SetUpBoxes(AbstractTetrahedralMesh<DIM, DIM>& rMesh,
                                         double boxWidth,
                                         DistributedBoxCollection<DIM>*& rpBoxCollection)
{
    if (rpBoxCollection)
    {
        delete rpBoxCollection;
        rpBoxCollection = nullptr;
    }

    // Compute min and max values for the fine mesh nodes
    ChasteCuboid<DIM> bounding_box = rMesh.CalculateBoundingBox();

    // Set up the boxes using a domain that is slightly larger than the fine mesh
    c_vector<double,2*DIM> extended_min_and_max;
    for (unsigned i=0; i<DIM; i++)
    {
        double width = bounding_box.GetWidth(i);

        // Subtract from the minima
        extended_min_and_max(2*i) = bounding_box.rGetLowerCorner()[i] - 0.05*width;

        // Add to the maxima
        extended_min_and_max(2*i+1) = bounding_box.rGetUpperCorner()[i] + 0.05*width;
    }

    if (boxWidth < 0)
    {
        /*
         * Use default value = max(max_edge_length, w20),  where w20 is the width
         * corresponding to 20 boxes in the x-direction.
         *
         * BoxCollection creates an extra box so divide by 19 not 20. Add a little
         * bit on to ensure minor numerical fluctuations don't change the answer.
         */
        boxWidth = (extended_min_and_max(1) - extended_min_and_max(0))/19.000000001;

        // Determine the maximum edge length
        c_vector<double, 2> min_max_edge_length = rMesh.CalculateMinMaxEdgeLengths();

        if (boxWidth < min_max_edge_length[1])
        {
            boxWidth = 1.1*min_max_edge_length[1];
        }
    }

    rpBoxCollection = new DistributedBoxCollection<DIM>(boxWidth, extended_min_and_max);
    rpBoxCollection->SetupAllLocalBoxes();

    // For each element, if ANY of its nodes are physically in a box, put that element in that box
    for (unsigned i=0; i<rMesh.GetNumElements(); i++)
    {
        Element<DIM,DIM>* p_element = rMesh.GetElement(i);

        std::set<unsigned> box_indices_each_node_this_elem;
        for (unsigned j=0; j<DIM+1; j++) // num vertices per element
        {
            Node<DIM>* p_node = p_element->GetNode(j);
            // Only take note of box inclusions which are in our domain
            if (rpBoxCollection->IsOwned(p_node))
            {
                unsigned box_index = rpBoxCollection->CalculateContainingBox(p_node);
                box_indices_each_node_this_elem.insert(box_index);
            }
        }
        for (std::set<unsigned>::iterator iter = box_indices_each_node_this_elem.begin();
            iter != box_indices_each_node_this_elem.end();
            ++iter)
        {
            assert(rpBoxCollection->IsBoxOwned( *iter ));
            rpBoxCollection->rGetBox( *iter ).AddElement(p_element);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////
// ComputeFineElementsAndWeightsForCoarseQuadPoints()
// and
// ComputeFineElementsAndWeightsForCoarseNodes()
// and
// common method
////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule,
                                                                               bool safeMode)
{
    if (mpFineMeshBoxCollection == nullptr)
    {
        EXCEPTION("Call SetUpBoxesOnFineMesh() before ComputeFineElementsAndWeightsForCoarseQuadPoints()");
    }

    // Get the quad point (physical) positions
    QuadraturePointsGroup<DIM> quad_point_posns(mrCoarseMesh, rQuadRule);

    // Resize the elements and weights vector.
    mFineMeshElementsAndWeights.resize(quad_point_posns.Size());
    // Make sure that, in parallel, silent processes have their structs initialised to zero values
    for (unsigned i=0; i<mFineMeshElementsAndWeights.size(); i++)
    {
        mFineMeshElementsAndWeights[i].ElementNum = 0u;
        mFineMeshElementsAndWeights[i].Weights = zero_vector<double>(DIM+1);
    }


    // LCOV_EXCL_START
    if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
    {
        std::cout << "\nComputing fine elements and weights for coarse quad points\n";
    }
    // LCOV_EXCL_STOP


    ResetStatisticsVariables();
    for (unsigned i=0; i<quad_point_posns.Size(); i++)
    {
        // LCOV_EXCL_START
        if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
        {
            std::cout << "\t" << i << " of " << quad_point_posns.Size() << std::flush;
        }
        // LCOV_EXCL_STOP

        // Get the box this point is in
        unsigned box_for_this_point = mpFineMeshBoxCollection->CalculateContainingBox( quad_point_posns.rGet(i) );
        if (mpFineMeshBoxCollection->IsBoxOwned(box_for_this_point))
        {
            // A chaste point version of the c-vector is needed for the GetContainingElement call.
            ChastePoint<DIM> point(quad_point_posns.rGet(i));

            ComputeFineElementAndWeightForGivenPoint(point, safeMode, box_for_this_point, i);
        }
        else
        {
            assert(mFineMeshElementsAndWeights[i].ElementNum == 0u);
            assert(norm_2(mFineMeshElementsAndWeights[i].Weights) == 0.0 );
        }
    }
    ShareFineElementData();
    if (mStatisticsCounters[1] > 0)
    {
        WARNING(mStatisticsCounters[1] << " of " << quad_point_posns.Size() << " coarse-mesh quadrature points were outside the fine mesh");
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementsAndWeightsForCoarseNodes(bool safeMode)
{
    if (mpFineMeshBoxCollection==nullptr)
    {
        EXCEPTION("Call SetUpBoxesOnFineMesh() before ComputeFineElementsAndWeightsForCoarseNodes()");
    }

    // Resize the elements and weights vector.
    mFineMeshElementsAndWeights.resize(mrCoarseMesh.GetNumNodes());
    // Make sure that, in parallel, silent processes have their structs initialised to zero values
    for (unsigned i=0; i<mFineMeshElementsAndWeights.size(); i++)
    {
        mFineMeshElementsAndWeights[i].ElementNum = 0u;
        mFineMeshElementsAndWeights[i].Weights = zero_vector<double>(DIM+1);
    }

    // LCOV_EXCL_START
    if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
    {
        std::cout << "\nComputing fine elements and weights for coarse nodes\n";
    }
    // LCOV_EXCL_STOP


    ResetStatisticsVariables();
    for (unsigned i=0; i<mrCoarseMesh.GetNumNodes(); i++)
    {
        // LCOV_EXCL_START
        if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
        {
            std::cout << "\t" << i << " of " << mrCoarseMesh.GetNumNodes() << std::flush;
        }
        // LCOV_EXCL_STOP

        Node<DIM>* p_node = mrCoarseMesh.GetNode(i);

        // Get the box this point is in
        unsigned box_for_this_point = mpFineMeshBoxCollection->CalculateContainingBox( p_node->rGetModifiableLocation() );
        if (mpFineMeshBoxCollection->IsBoxOwned(box_for_this_point))
        {
            // A chaste point version of the c-vector is needed for the GetContainingElement call
            ChastePoint<DIM> point(p_node->rGetLocation());

            ComputeFineElementAndWeightForGivenPoint(point, safeMode, box_for_this_point, i);
        }
    }
    ShareFineElementData();
}

/**
 * \todo: could possibly merge with ComputeCoarseElementForGivenPoint(). Difference between
 * the methods are: this uses fine mesh and fine mesh box, computes weights as well (and sets
 * the element and weight in the vec), rather than returning the element, and this method
 * saves information in mStatisticsCounters.
 */
template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementAndWeightForGivenPoint(ChastePoint<DIM>& rPoint,
                                                                       bool safeMode,
                                                                       unsigned boxForThisPoint,
                                                                       unsigned index)
{
    std::set<unsigned> test_element_indices;

    /*
     * The elements to try (initially) are those contained in the box the point is in.
     *
     * Note: it is possible the point to be in an element not 'in' this box, as it is
     * possible for all element nodes to be in different boxes.
     */
    CollectElementsInContainingBox(mpFineMeshBoxCollection, boxForThisPoint, test_element_indices);

    unsigned elem_index;
    c_vector<double,DIM+1> weight;

    try
    {
        //std::cout << "\n" << "# test elements initially " << test_element_indices.size() << "\n";
        // Try these elements only, initially
        elem_index = mrFineMesh.GetContainingElementIndex(rPoint,
                                                          false,
                                                          test_element_indices,
                                                          true /* quit if not in test_elements */);
        weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

        mStatisticsCounters[0]++;
    }
    catch(Exception&) //not_in_box
    {
        // Now try all elements, trying the elements contained in the boxes local to this element first
        test_element_indices.clear();

        CollectElementsInLocalBoxes(mpFineMeshBoxCollection, boxForThisPoint, test_element_indices);

        try
        {
            elem_index = mrFineMesh.GetContainingElementIndex(rPoint, false, test_element_indices, true);
            weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);
            mStatisticsCounters[0]++;
        }
        catch(Exception&) //not_in_local_boxes
        {
            if (safeMode)
            {
                // Try the remaining elements
                try
                {
                    elem_index = mrFineMesh.GetContainingElementIndex(rPoint, false);
                    weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);
                    mStatisticsCounters[0]++;

                }
                catch (Exception&) // not_in_mesh
                {
                    // The point is not in ANY element, so store the nearest element and corresponding weights
                    elem_index = mrFineMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                    weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

                    mNotInMesh.push_back(index);
                    mNotInMeshNearestElementWeights.push_back(weight);
                    mStatisticsCounters[1]++;
                }
            }
            else
            {
                assert(test_element_indices.size() > 0); // boxes probably too small if this fails

                /*
                 * Immediately assume it isn't in the rest of the mesh - this should be the
                 * case assuming the box width was chosen suitably. Store the nearest element
                 * and corresponding weights.
                 */
                elem_index = mrFineMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

                mNotInMesh.push_back(index);
                mNotInMeshNearestElementWeights.push_back(weight);
                mStatisticsCounters[1]++;
            }
        }
    }

    mFineMeshElementsAndWeights[index].ElementNum = elem_index;
    mFineMeshElementsAndWeights[index].Weights = weight;
}

////////////////////////////////////////////////////////////////////////////////////
// ComputeCoarseElementsForFineNodes
// and
// ComputeCoarseElementsForFineElementCentroids
// and
// common method
////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeCoarseElementsForFineNodes(bool safeMode)
{
    if (mpCoarseMeshBoxCollection==nullptr)
    {
        EXCEPTION("Call SetUpBoxesOnCoarseMesh() before ComputeCoarseElementsForFineNodes()");
    }

    // LCOV_EXCL_START
    if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
    {
        std::cout << "\nComputing coarse elements for fine nodes\n";
    }
    // LCOV_EXCL_STOP
    mCoarseElementsForFineNodes.clear();
    mCoarseElementsForFineNodes.resize(mrFineMesh.GetNumNodes(), 0.0);

    ResetStatisticsVariables();
    for (unsigned i=0; i<mCoarseElementsForFineNodes.size(); i++)
    {
        // LCOV_EXCL_START
        if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
        {
            std::cout << "\t" << i << " of " << mCoarseElementsForFineNodes.size() << std::flush;
        }
        // LCOV_EXCL_STOP

        ChastePoint<DIM> point = mrFineMesh.GetNode(i)->GetPoint();

        // Get the box this point is in
        unsigned box_for_this_point = mpCoarseMeshBoxCollection->CalculateContainingBox(mrFineMesh.GetNode(i)->rGetModifiableLocation());
        if (mpCoarseMeshBoxCollection->IsBoxOwned(box_for_this_point))
        {
            mCoarseElementsForFineNodes[i] = ComputeCoarseElementForGivenPoint(point, safeMode, box_for_this_point);
        }
    }
    ShareCoarseElementData();
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeCoarseElementsForFineElementCentroids(bool safeMode)
{
    if (mpCoarseMeshBoxCollection==nullptr)
    {
        EXCEPTION("Call SetUpBoxesOnCoarseMesh() before ComputeCoarseElementsForFineElementCentroids()");
    }

    // LCOV_EXCL_START
    if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
    {
        std::cout << "\nComputing coarse elements for fine element centroids\n";
    }
    // LCOV_EXCL_STOP

    mCoarseElementsForFineElementCentroids.clear();
    mCoarseElementsForFineElementCentroids.resize(mrFineMesh.GetNumElements(), 0.0);

    ResetStatisticsVariables();
    for (unsigned i=0; i<mrFineMesh.GetNumElements(); i++)
    {
        // LCOV_EXCL_START
        if (CommandLineArguments::Instance()->OptionExists("-mesh_pair_verbose"))
        {
            std::cout << "\t" << i << " of " << mrFineMesh.GetNumElements() << std::flush;
        }
        // LCOV_EXCL_STOP

        c_vector<double,DIM> point_cvec = mrFineMesh.GetElement(i)->CalculateCentroid();
        ChastePoint<DIM> point(point_cvec);

        // Get the box this point is in
        unsigned box_for_this_point = mpCoarseMeshBoxCollection->CalculateContainingBox( point_cvec );

        if (mpCoarseMeshBoxCollection->IsBoxOwned(box_for_this_point))
        {
            mCoarseElementsForFineElementCentroids[i] = ComputeCoarseElementForGivenPoint(point, safeMode, box_for_this_point);
        }
    }
    ShareCoarseElementData();
}

template<unsigned DIM>
unsigned FineCoarseMeshPair<DIM>::ComputeCoarseElementForGivenPoint(ChastePoint<DIM>& rPoint,
                                                                    bool safeMode,
                                                                    unsigned boxForThisPoint)
{
    /**
     * \todo: could possibly merge with ComputeFineElementAndWeightForGivenPoint(). Differences
     * between the methods are: the other method uses fine mesh and fine mesh box, computes
     * weights as well (and sets the element and weight in the vec), rather than returning the
     * element, and that method saves information in mStatisticsCounters.
     */
    std::set<unsigned> test_element_indices;
    CollectElementsInContainingBox(mpCoarseMeshBoxCollection, boxForThisPoint, test_element_indices);

    unsigned elem_index;

    try
    {
        elem_index = mrCoarseMesh.GetContainingElementIndex(rPoint,
                                                            false,
                                                            test_element_indices,
                                                            true /* quit if not in test_elements */);

        mStatisticsCounters[0]++;
    }
    catch(Exception&) // not_in_box
    {
        // Now try all elements, trying the elements contained in the boxes local to this element first
        test_element_indices.clear();
        CollectElementsInLocalBoxes(mpCoarseMeshBoxCollection, boxForThisPoint, test_element_indices);

        try
        {
            elem_index = mrCoarseMesh.GetContainingElementIndex(rPoint, false, test_element_indices, true);
            mStatisticsCounters[0]++;
        }
        catch(Exception&) // not_in_local_boxes
        {
            if (safeMode)
            {
                // Try the remaining elements
                try
                {
                    elem_index = mrCoarseMesh.GetContainingElementIndex(rPoint, false);

                    mStatisticsCounters[0]++;
                }
                catch (Exception&) // not_in_mesh
                {
                    // The point is not in ANY element, so store the nearest element and corresponding weights
                    elem_index = mrCoarseMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                    mStatisticsCounters[1]++;
                }
            }
            else
            {
                assert(test_element_indices.size() > 0); // boxes probably too small if this fails

                /*
                 * Immediately assume it isn't in the rest of the mesh - this should be the
                 * case assuming the box width was chosen suitably. Store the nearest element
                 * and corresponding weights.
                 */
                elem_index = mrCoarseMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                mStatisticsCounters[1]++;
            }
        }
    }

    return elem_index;
}

////////////////////////////////////////////////////////////////////////////////////
// Helper methods for code
////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::CollectElementsInContainingBox(DistributedBoxCollection<DIM>*& rpBoxCollection,
                                                             unsigned boxIndex,
                                                             std::set<unsigned>& rElementIndices)
{
    for (typename std::set<Element<DIM,DIM>*>::iterator elem_iter = rpBoxCollection->rGetBox(boxIndex).rGetElementsContained().begin();
         elem_iter != rpBoxCollection->rGetBox(boxIndex).rGetElementsContained().end();
         ++elem_iter)
    {
        rElementIndices.insert((*elem_iter)->GetIndex());
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::CollectElementsInLocalBoxes(DistributedBoxCollection<DIM>*& rpBoxCollection,
                                                          unsigned boxIndex,
                                                          std::set<unsigned>& rElementIndices)
{
    std::set<unsigned> local_boxes = rpBoxCollection->rGetLocalBoxes(boxIndex);
    for (std::set<unsigned>::iterator local_box_iter = local_boxes.begin();
         local_box_iter != local_boxes.end();
         ++local_box_iter)
    {
        for (typename std::set<Element<DIM,DIM>*>::iterator elem_iter = rpBoxCollection->rGetBox(*local_box_iter).rGetElementsContained().begin();
             elem_iter != rpBoxCollection->rGetBox(*local_box_iter).rGetElementsContained().end();
             ++elem_iter)
        {
            rElementIndices.insert((*elem_iter)->GetIndex());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////
// Statistics related methods
////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ResetStatisticsVariables()
{
    mNotInMesh.clear();
    mNotInMeshNearestElementWeights.clear();
    mStatisticsCounters.resize(2, 0u);
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ShareFineElementData()
{
    if (PetscTools::IsSequential())
    {
        return;
    }

    // This sums the results so it isn't idempotent: you get a different result if you call this method twice
    // Should not matter: the methods which call this helper method have reset everything which is about to be shared
    std::vector<unsigned> local_counters = mStatisticsCounters;
    MPI_Allreduce(&local_counters[0], &mStatisticsCounters[0], 2u, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);


    // Get all the element number and weights into a contiguous format
    unsigned elements_size = mFineMeshElementsAndWeights.size();
    std::vector<unsigned> all_element_indices(elements_size);
    unsigned weights_size = elements_size*(DIM+1);
    std::vector<double> all_weights(weights_size);
    for (unsigned index=0; index<mFineMeshElementsAndWeights.size(); index++)
    {
        all_element_indices[index] = mFineMeshElementsAndWeights[index].ElementNum;
        for (unsigned j=0; j<DIM+1; j++)
        {
            all_weights[index*(DIM+1)+j] = mFineMeshElementsAndWeights[index].Weights[j];
        }
    }

    // Copy and share
    std::vector<unsigned> local_all_element_indices = all_element_indices;
    std::vector<double> local_all_weights = all_weights;
    MPI_Allreduce(&local_all_element_indices[0], &all_element_indices[0], elements_size, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce( &local_all_weights[0], &all_weights[0], weights_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    // Put back into the regular data structure
    for (unsigned index=0; index<mFineMeshElementsAndWeights.size(); index++)
    {
        mFineMeshElementsAndWeights[index].ElementNum = all_element_indices[index];
        for (unsigned j=0; j<DIM+1; j++)
        {
            mFineMeshElementsAndWeights[index].Weights[j] = all_weights[index*(DIM+1)+j] ;
        }
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ShareCoarseElementData()
{
    if (PetscTools::IsSequential())
    {
        return;
    }

    // This sums the results so it isn't idempotent: you get a different result if you call this method twice
    // Should not matter: the methods which call this helper method have reset #mStatisticsCounters
    std::vector<unsigned> local_counters = mStatisticsCounters;
    MPI_Allreduce(&local_counters[0], &mStatisticsCounters[0], 2u, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

    // The rest uses "max" so it is idempotent.  You can safely re-share results between processes without them changing.
    if (mCoarseElementsForFineNodes.empty() == false)
    {
        std::vector<unsigned> temp_coarse_elements = mCoarseElementsForFineNodes;
        MPI_Allreduce( &temp_coarse_elements[0], &mCoarseElementsForFineNodes[0], mCoarseElementsForFineNodes.size(), MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
    }
    if (mCoarseElementsForFineElementCentroids.empty() == false)
    {
        std::vector<unsigned> temp_coarse_elements = mCoarseElementsForFineElementCentroids;
        MPI_Allreduce( &temp_coarse_elements[0], &mCoarseElementsForFineElementCentroids[0], mCoarseElementsForFineElementCentroids.size(), MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::PrintStatistics()
{
    std::cout << "\nFineCoarseMeshPair statistics for the last-called method:\n";

//    std::cout << "\tNum points for which containing element was found, using box containing that point = " << mStatisticsCounters[0] << "\n";
//    std::cout << "\tNum points for which containing element was in local box = " << mStatisticsCounters[1] << "\n";
//    std::cout << "\tNum points for which containing element was in an element in a non-local box = " << mStatisticsCounters[2] << "\n";
//    std::cout << "\tNum points for which no containing element was found = " << mStatisticsCounters[3] << "\n";

    std::cout << "\tNum points for which containing element was found: " << mStatisticsCounters[0] << "\n";
    std::cout << "\tNum points for which no containing element was found = " << mStatisticsCounters[1] << "\n";

    if (mNotInMesh.size() > 0)
    {
        std::cout << "\tIndices and weights for points (nodes/quad points) for which no containing element was found:\n";
        for (unsigned i=0; i<mNotInMesh.size(); i++)
        {
            std::cout << "\t\t" << mNotInMesh[i] << ", " << mNotInMeshNearestElementWeights[i] << "\n";
        }
    }
}

///////// Explicit instantiation///////

template class FineCoarseMeshPair<1>;
template class FineCoarseMeshPair<2>;
template class FineCoarseMeshPair<3>;
