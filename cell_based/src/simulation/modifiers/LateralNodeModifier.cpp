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

#include "LateralNodeModifier.hpp"
#include <set>
#include <vector>
#include "Exception.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include "MutableVertexMesh.hpp"
#include "Node.hpp"
#include "VertexElement.hpp"

#include "Debug.hpp"

void LateralNodeModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<3, 3>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

void LateralNodeModifier::SetupSolve(AbstractCellPopulation<3, 3>& rCellPopulation, std::string outputDirectory)
{
}

void LateralNodeModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}

void LateralNodeModifier::UpdateCellData(AbstractCellPopulation<3, 3>& rCellPopulation)
{
    AbstractMesh<3, 3>* p_tmp_mesh = &(rCellPopulation.rGetMesh());

    if (dynamic_cast<MutableVertexMesh<3, 3>*>(p_tmp_mesh) == nullptr)
    {
        EXCEPTION("only for vertex mesh"); //LCOV_EXCL_LINE
    }
    MutableVertexMesh<3, 3>* const p_mesh = static_cast<MutableVertexMesh<3, 3>*>(p_tmp_mesh);

    for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
    {
        Node<3>* p_lateral_node = p_mesh->GetNode(i);
        if (!IsLateralNode(p_lateral_node))
        {
            continue;
        }

        const std::set<VertexElement<2, 3>*> lateral_faces = GetFacesWithIndices(p_lateral_node->rGetContainingFaceIndices(), p_mesh, Monolayer::LateralValue);
        std::vector<Node<3>*> basal_nodes;
        std::vector<Node<3>*> apical_nodes;

        for (VertexElement<2, 3>* p_face_tmp : lateral_faces)
        {
            if (p_face_tmp->GetNumNodes() != 3u)
            {
                continue;
            }

            const std::vector<Node<3>*> apical_nodes_tmp = GetNodesWithType(p_face_tmp, Monolayer::ApicalValue);
            apical_nodes.insert(apical_nodes.end(), apical_nodes_tmp.begin(), apical_nodes_tmp.end());
            const std::vector<Node<3>*> basal_nodes_tmp = GetNodesWithType(p_face_tmp, Monolayer::BasalValue);
            basal_nodes.insert(basal_nodes.end(), basal_nodes_tmp.begin(), basal_nodes_tmp.end());
        }

        if (basal_nodes.size() != 2u || apical_nodes.size() != 2u)
        {
            NEVER_REACHED;
        }

        const double apical_length = p_mesh->GetDistanceBetweenNodes(apical_nodes[0]->GetIndex(), apical_nodes[1]->GetIndex());
        const double basal_length = p_mesh->GetDistanceBetweenNodes(basal_nodes[0]->GetIndex(), basal_nodes[1]->GetIndex());

        /*
         * #480 The interface (lateral) node is no longer repositioned to a geometric midpoint here.
         * Instead it is left to move under the ordinary force-based dynamics, which minimise the
         * GeneralMonolayerVertexMeshForce energy. This represents a scutoid whose partial (asynchronous)
         * T1 swap can "unzipper" up or down the lateral direction: when the resulting apical or basal
         * edge falls below the cell rearrangement threshold, the swap is resolved into a full T1 below.
         */

        // Resolve the asynchronous T1 (unzipper to a full T1) once the other edge falls below mCellRearrangementThreshold
        if (std::min(apical_length, basal_length) < p_mesh->GetCellRearrangementThreshold())
        {
            const bool reverse_t1_on_basal = basal_length < p_mesh->GetCellRearrangementThreshold();
            const Monolayer::v_type T1_type = reverse_t1_on_basal ? Monolayer::BasalValue : Monolayer::ApicalValue;
            Node<3>* p_node_a = reverse_t1_on_basal ? basal_nodes[0] : apical_nodes[0];
            Node<3>* p_node_b = reverse_t1_on_basal ? basal_nodes[1] : apical_nodes[1];
            const double distance_ab = reverse_t1_on_basal ? basal_length : apical_length;

            Node<3>* p_node_x = reverse_t1_on_basal ? apical_nodes[0] : basal_nodes[0];
            Node<3>* p_node_y = reverse_t1_on_basal ? apical_nodes[1] : basal_nodes[1];
            const double distance_xy = reverse_t1_on_basal ? apical_length : basal_length;

            VertexElement<2, 3>* p_small_triangular_face = GetSharedLateralFace(p_node_a, p_node_b, p_mesh);
            if (p_small_triangular_face == nullptr || p_small_triangular_face->GetNumNodes() != 3u || p_small_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
            {
                NEVER_REACHED;
            }

            VertexElement<2, 3>* p_big_triangular_face = GetSharedLateralFace(p_node_x, p_node_y, p_mesh);
            if (p_big_triangular_face == nullptr || p_big_triangular_face->GetNumNodes() != 3u || p_big_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
            {
                NEVER_REACHED;
            }

            MARK;
            TRACE("=============================== Reverse Asynchronous ============================================");
            PRINT_2_VARIABLES(p_small_triangular_face->GetIndex(), p_big_triangular_face->GetIndex());

            const c_vector<double, 3> vector_ab = p_mesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            const c_vector<double, 3> vector_xy = p_mesh->GetVectorFromAtoB(p_node_x->rGetLocation(), p_node_y->rGetLocation());

            if (!((distance_ab < p_mesh->GetCellRearrangementThreshold()) && (distance_xy > p_mesh->GetCellRearrangementThreshold())))
            {
                NEVER_REACHED;
            }

            // Using more complicated way to calculate midpoint rather than
            // 0.5 * (p_node_a->rGetLocation() + p_node_b->rGetLocation())
            // for more general cases (like cylindrical mesh).
            const c_vector<double, 3> mid_ab = p_node_a->rGetLocation()
                + 0.5 * p_mesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            // Compute and store the location of the T1 swap, which is at the midpoint of nodes A and B
            // mLocationsOfT1Swaps.push_back(mid_ab);

            std::vector<VertexElement<3, 3>*> elems(5, nullptr);
            {
                const std::set<unsigned>& elem_24 = p_small_triangular_face->rFaceGetContainingElementIndices();
                const std::set<unsigned>& elem_a_124 = p_node_a->rGetContainingElementIndices();
                const std::set<unsigned>& elem_b_234 = p_node_b->rGetContainingElementIndices();

                if (elem_24.size() == 0 || elem_24.size() > 2)
                {
                    // in case of 0, it should be T3, not T1
                    // in case of >2, error somewhere.
                    NEVER_REACHED;
                }

                elems[2] = p_mesh->GetElement(no1(elem_24));
                if (elem_24.size() == 2)
                {
                    elems[4] = p_mesh->GetElement(no2(elem_24));
                }

                // Assign element 1 and 3 if exist.
                if (elem_a_124.size() - elem_24.size() == 1)
                {
                    const std::set<unsigned> s_tmp = elem_a_124 - elem_24;
                    assert(s_tmp.size() == 1);
                    elems[1] = p_mesh->GetElement(no1(s_tmp));
                }
                if (elem_b_234.size() - elem_24.size() == 1)
                {
                    const std::set<unsigned> s_tmp = elem_b_234 - elem_24;
                    assert(s_tmp.size() == 1);
                    elems[3] = p_mesh->GetElement(no1(s_tmp));
                }
            }

            // Now since we know the lateral swap face, we can move the nodes using the face normal.
            {
                c_vector<double, 3> vector_cd = p_mesh->GetCellRearrangementRatio() * p_mesh->GetCellRearrangementThreshold() * vector_xy / norm_2(vector_xy);

                const c_vector<double, 3> centroid_2 = GetFacesWithType(elems[2], T1_type)[0]->GetCentroid();
                const c_vector<double, 3> centroid_4 = (elems[4] != nullptr) ? GetFacesWithType(elems[4], T1_type)[0]->GetCentroid()
                                                                          : p_small_triangular_face->GetCentroid();
                const c_vector<double, 3> vector_from2 = p_mesh->GetVectorFromAtoB(centroid_2, centroid_4);

                // change sign accordingly
                if (inner_prod(vector_cd, vector_from2) < 0)
                {
                    vector_cd *= -1;
                }

                // Move nodes A and B to C and D respectively
                p_node_a->rGetModifiableLocation() = mid_ab - 0.5 * vector_cd;
                p_node_b->rGetModifiableLocation() = mid_ab + 0.5 * vector_cd;
            }

            // Start modifications

            // Remove the lateral node
            {
                const std::set<unsigned> faces_tmp = p_lateral_node->rGetContainingFaceIndices();
                for (unsigned face_index : faces_tmp)
                {
                    p_mesh->GetFace(face_index)->FaceDeleteNode(p_lateral_node);
                }

                const std::set<unsigned> elems_tmp = p_lateral_node->rGetContainingElementIndices();
                for (unsigned elem_index : elems_tmp)
                {
                    p_mesh->GetElement(elem_index)->DeleteNode(p_lateral_node);
                }
                p_mesh->DeleteNodePriorToReMesh(p_lateral_node->GetIndex());
            }

            // Settle the swap face.
            p_big_triangular_face->FaceAddNode(p_node_a);
            p_big_triangular_face->FaceAddNode(p_node_b);

            FaceRearrangeNodesInMesh(p_mesh, p_big_triangular_face);

            p_small_triangular_face->FaceDeleteNode(p_node_a);
            p_small_triangular_face->FaceDeleteNode(p_node_b);

            // Modify lateral face "share" by element 2 and element 3
            {
                const std::set<unsigned> tmp_face_ids = p_node_b->rGetContainingFaceIndices();

                std::set<VertexElement<2, 3>*> s_tmp = GetFacesWithIndices(tmp_face_ids, elems[2], Monolayer::LateralValue);
                // p_node_b is already removed from p_small_triangular_face
                if (s_tmp.erase(p_small_triangular_face) != 0)
                {
                    NEVER_REACHED;
                }
                assert(s_tmp.size() == 1);
                VertexElement<2, 3>* const p_lateral_face_23 = no1(s_tmp);

                p_lateral_face_23->FaceUpdateNode(p_node_b, p_node_a);
                FaceRearrangeNodesInMesh(p_mesh, p_lateral_face_23);
            }

            // Modify lateral face "share" by element 1 and element 4
            {
                VertexElement<2, 3>* p_lateral_face_14 = nullptr;
                const std::set<unsigned>& tmp_face_ids = p_node_a->rGetContainingFaceIndices();
                if (elems[4] != nullptr)
                {
                    std::set<VertexElement<2, 3>*> s_tmp = GetFacesWithIndices(tmp_face_ids, elems[4], Monolayer::LateralValue);
                    // p_node_a is already removed from p_small_triangular_face
                    if (s_tmp.erase(p_small_triangular_face) != 0)
                    {
                        NEVER_REACHED;
                    }
                    assert(s_tmp.size() == 1);
                    p_lateral_face_14 = no1(s_tmp);
                }
                else
                {
                    const std::set<VertexElement<2, 3>*> s_tmp = GetFacesWithIndices(tmp_face_ids, p_mesh, Monolayer::LateralValue)
                        - GetFacesWithIndices(tmp_face_ids, elems[2], Monolayer::LateralValue);
                    assert(s_tmp.size() == 1);
                    p_lateral_face_14 = no1(s_tmp);
                }
                p_lateral_face_14->FaceUpdateNode(p_node_a, p_node_b);
                FaceRearrangeNodesInMesh(p_mesh, p_lateral_face_14);
            }

            for (unsigned i = 1; i <= 4; ++i)
            {
                if (elems[i] == nullptr)
                {
                    continue;
                }

                VertexElement<3, 3>* p_elem = elems[i];
                VertexElement<2, 3>* p_face = GetFacesWithType(p_elem, T1_type)[0];

                switch (i)
                {
                    case 1:
                    {
                        p_face->FaceAddNode(p_node_b, p_face->GetNumNodes() - 1);
                        FaceRearrangeNodesInMesh(p_mesh, p_face);
                        p_elem->AddNode(p_node_b, p_elem->GetNumNodes() - 1);
                        break;
                    }
                    case 2:
                    {
                        p_face->FaceDeleteNode(p_node_b);
                        p_elem->DeleteNode(p_node_b);
                        p_elem->DeleteFace(p_small_triangular_face);
                        break;
                    }
                    case 3:
                    {
                        p_face->FaceAddNode(p_node_a, p_face->GetNumNodes() - 1);
                        FaceRearrangeNodesInMesh(p_mesh, p_face);
                        p_elem->AddNode(p_node_a, p_elem->GetNumNodes() - 1);
                        break;
                    }
                    case 4:
                    {
                        p_face->FaceDeleteNode(p_node_a);
                        p_elem->DeleteNode(p_node_a);
                        p_elem->DeleteFace(p_small_triangular_face);
                        break;
                    }
                    default:
                        NEVER_REACHED;
                }

                p_elem->MonolayerElementRearrangeFacesNodes();
            }

            // Sort out boundary nodes
            if (p_node_a->IsBoundaryNode() || p_node_b->IsBoundaryNode())
            {
                if (p_node_a->GetNumContainingElements() == 3)
                {
                    p_node_a->SetAsBoundaryNode(false);
                }
                else
                {
                    p_node_a->SetAsBoundaryNode(true);
                }
                if (p_node_b->GetNumContainingElements() == 3)
                {
                    p_node_b->SetAsBoundaryNode(false);
                }
                else
                {
                    p_node_b->SetAsBoundaryNode(true);
                }
            }

            p_mesh->DeleteFacePriorToReMesh(p_small_triangular_face->GetIndex());
        } // end of async
    }
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(LateralNodeModifier)
