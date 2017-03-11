/*
 * LateralNodeModifier.cpp
 *
 *  Created on: 11 Mar 2017
 *      Author: Weijie
 */

#include "LateralNodeModifier.hpp"
#include "Node.hpp"
#include "VertexElement.hpp"
#include "MutableVertexMesh.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "Exception.hpp"
#include <set>
#include <vector>


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
    AbstractMesh<3, 3>& r_tmp_mesh = rCellPopulation.rGetMesh();

    if (dynamic_cast<MutableVertexMesh<3, 3>*>(&r_tmp_mesh) == NULL)
    {
        EXCEPTION("only to vertex mesh");
    }
    MutableVertexMesh<3, 3>* const p_mesh = static_cast<MutableVertexMesh<3, 3>*>(&r_tmp_mesh);

    for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
    {
        Node<3>* p_lateral_node = p_mesh->GetNode(i);
        if (IsLateralNode(p_lateral_node))
        {
            const std::set<VertexElement<2, 3>*> lateral_faces = GetFacesWithIndices(p_lateral_node->rGetContainingFaceIndices(), p_mesh, Monolayer::LateralValue);
            std::vector<Node<3>*> basal_nodes;
            std::vector<Node<3>*> apical_nodes;

            for (std::set<VertexElement<2, 3>*>::const_iterator it = lateral_faces.begin();
                 it != lateral_faces.end(); ++it)
            {
                VertexElement<2, 3>* p_face_tmp = *it;
                if (p_face_tmp->GetNumNodes() != 3u)
                    continue;

                for (unsigned j = 0; j < p_face_tmp->GetNumNodes(); ++j)
                {
                    Node<3>* p_tmp_node = p_face_tmp->GetNode(j);
                    switch (GetNodeType(p_tmp_node))
                    {
                        case Monolayer::ApicalValue:
                            apical_nodes.push_back(p_tmp_node);
                            break;
                        case Monolayer::BasalValue:
                            basal_nodes.push_back(p_tmp_node);
                            break;
                        default:
                            assert(p_tmp_node == p_lateral_node);
                    }
                }
            }

            if (basal_nodes.size() != 2u || apical_nodes.size() != 2u)
            {
                NEVER_REACHED;
            }

            const c_vector<double, 3> mid_apical = (apical_nodes[0]->rGetLocation() + apical_nodes[1]->rGetLocation()) / 2;
            const c_vector<double, 3> mid_basal = (basal_nodes[0]->rGetLocation() + basal_nodes[1]->rGetLocation()) / 2;

            const double apical_length = norm_2(apical_nodes[0]->rGetLocation() - apical_nodes[1]->rGetLocation());
            const double basal_length = norm_2(basal_nodes[0]->rGetLocation() - basal_nodes[1]->rGetLocation());

            // HOW_MANY_TIMES_HERE("bla bla bla");
            // PRINT_CONTAINER(p_lateral_node->rGetLocation());
            p_lateral_node->rGetModifiableLocation() = (basal_length * mid_apical + apical_length * mid_basal) / (apical_length + basal_length);
            // PRINT_CONTAINER(p_lateral_node->rGetLocation());

            // Reverse AsynchronousT1 when the other edge is shorter than mCellRearrangementThreshold
            if (std::min(apical_length, basal_length) < p_mesh->GetCellRearrangementThreshold())
            {
                bool reverse_t1_on_basal = basal_length < p_mesh->GetCellRearrangementThreshold();
                Node<3>* p_node_a = reverse_t1_on_basal ? basal_nodes[0] : apical_nodes[0];
                Node<3>* p_node_b = reverse_t1_on_basal ? basal_nodes[1] : apical_nodes[1];

                Node<3>* p_node_x = reverse_t1_on_basal ? apical_nodes[0] : basal_nodes[0];
                Node<3>* p_node_y = reverse_t1_on_basal ? apical_nodes[1] : basal_nodes[1];

                VertexElement<2, 3>* p_small_triangular_face = GetSharedLateralFace(p_node_a, p_node_b, p_mesh);
                if (p_small_triangular_face == NULL || p_small_triangular_face->GetNumNodes() != 3u || p_small_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
                {
                    NEVER_REACHED;
                }

                VertexElement<2, 3>* p_big_triangular_face = GetSharedLateralFace(p_node_x, p_node_y, p_mesh);
                if (p_big_triangular_face == NULL || p_big_triangular_face->GetNumNodes() != 3u || p_big_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
                {
                    NEVER_REACHED;
                }

                MARK;
                TRACE("=============================== Reverse Asynchronous ============================================");
                PRINT_2_VARIABLES(p_small_triangular_face->GetIndex(), p_big_triangular_face->GetIndex());

                const c_vector<double, 3> vector_ab = p_mesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
                const c_vector<double, 3> vector_xy = p_mesh->GetVectorFromAtoB(p_node_x->rGetLocation(), p_node_y->rGetLocation());

                const unsigned node_a_index(p_node_a->GetIndex());
                const unsigned node_b_index(p_node_b->GetIndex());

                if (!((norm_2(vector_ab) < p_mesh->GetCellRearrangementThreshold()) && (norm_2(vector_xy) > p_mesh->GetCellRearrangementThreshold())))
                {
                    NEVER_REACHED;
                }

                // Compute and store the location of the T1 swap, which is at the midpoint of nodes A and B
                // mLocationsOfT1Swaps.push_back(p_node_a->rGetLocation() + 0.5 * vector_ab);

                std::vector<VertexElement<3, 3>*> elems(5, NULL);

                {
                    const std::set<unsigned>& elem_24 = p_small_triangular_face->rFaceGetContainingElementIndices();
                    const std::set<unsigned>& elem_a_124 = p_node_a->rGetContainingElementIndices();
                    const std::set<unsigned>& elem_b_234 = p_node_b->rGetContainingElementIndices();
                    MARK;
                    PRINT_CONTAINER(elem_24);
                    PRINT_CONTAINER(elem_a_124);
                    PRINT_CONTAINER(elem_b_234);

                    // Assign element 1 and 3 if exist.
                    if (elem_a_124.size() - elem_24.size() == 1)
                    {
                        std::set<unsigned> s_tmp;
                        std::set_difference(elem_a_124.begin(), elem_a_124.end(), elem_24.begin(),
                                            elem_24.end(), std::inserter(s_tmp, s_tmp.begin()));
                        assert(s_tmp.size() == 1);
                        PRINT_CONTAINER(s_tmp);
                        elems[1] = p_mesh->GetElement(no1(s_tmp));
                    }
                    if (elem_b_234.size() - elem_24.size() == 1)
                    {
                        std::set<unsigned> s_tmp;
                        std::set_difference(elem_b_234.begin(), elem_b_234.end(), elem_24.begin(),
                                            elem_24.end(), std::inserter(s_tmp, s_tmp.begin()));
                        assert(s_tmp.size() == 1);
                        PRINT_CONTAINER(s_tmp);
                        elems[3] = p_mesh->GetElement(no1(s_tmp));
                    }

                    MARK;
                    PRINT_VARIABLE(no1(elem_24));

                    if (elem_24.size() == 0)
                    {
                        // Should be T3, not T1
                        NEVER_REACHED;
                    }
                    VertexElement<3, 3>* p_elem = p_mesh->GetElement(no1(elem_24));
                    VertexElement<2, 3>* p_this_face = reverse_t1_on_basal ? GetBasalFace(p_elem) : GetApicalFace(p_elem);
                    const unsigned face_num_nodes = p_this_face->GetNumNodes();
                    const unsigned node_a_local_index = p_this_face->GetNodeLocalIndex(node_a_index);
                    const unsigned node_b_local_index = p_this_face->GetNodeLocalIndex(node_b_index);

                    PRINT_2_VARIABLES(node_a_local_index, node_b_local_index)
                    /*
                         * Locate local index of node_a and node_b and use the ordering to
                         * identify the element, if node_b_index > node_a_index then element 4
                         * and if node_a_index > node_b_index then element 2
                         */
                    if (node_a_local_index == plus1(node_b_local_index, face_num_nodes))
                    {
                        elems[2] = p_elem;
                        if (elem_24.size() == 2)
                        {
                            elems[4] = p_mesh->GetElement(*(++elem_24.begin()));
                        }
                    }
                    else if (node_b_local_index == plus1(node_a_local_index, face_num_nodes))
                    {
                        elems[4] = p_elem;
                        if (elem_24.size() == 2)
                        {
                            elems[2] = p_mesh->GetElement(*(++elem_24.begin()));
                        }
                    }
                    else
                    {
                        NEVER_REACHED;
                    }
                    // Assign 0th element as the 4th since pointer is cheap.
                    elems[0] = elems[4];
                }

                // Start modifications
                // Settle the swap face.
                p_big_triangular_face->FaceAddNode(p_node_a);
                p_big_triangular_face->FaceAddNode(p_node_b);
                p_big_triangular_face->FaceDeleteNode(p_lateral_node);

                FaceRearrangeNodesInMesh(p_mesh, p_big_triangular_face);

                p_small_triangular_face->FaceDeleteNode(p_node_a);
                p_small_triangular_face->FaceDeleteNode(p_node_b);
                p_small_triangular_face->FaceDeleteNode(p_lateral_node);

                {
                    std::vector<VertexElement<2, 3>*> lateral_faces(4, NULL);
                    for (unsigned i = 0; i < 4; ++i)
                    {
                        lateral_faces[i] = GetSharedLateralFace(elems[i], elems[i + 1]);
                    }

                    lateral_faces[2]->FaceUpdateNode(p_node_b, p_node_a);
                    lateral_faces[0]->FaceUpdateNode(p_node_a, p_node_b);

                    for (unsigned i = 0; i < 4; ++i)
                    {
                        if (lateral_faces[i] != NULL)
                        {
                            lateral_faces[i]->FaceDeleteNode(p_lateral_node);
                            FaceRearrangeNodesInMesh(p_mesh, lateral_faces[i]);
                        }
                    }
                }

                for (unsigned i = 1; i <= 4; ++i)
                {
                    if (elems[i] == NULL)
                    {
                        continue;
                    }

                    VertexElement<3, 3>* p_elem = elems[i];
                    VertexElement<2, 3>* p_face = reverse_t1_on_basal ? GetBasalFace(p_elem) : GetApicalFace(p_elem);

                    p_elem->DeleteNode(p_elem->GetNodeLocalIndex(p_lateral_node->GetIndex()));

                    switch (i)
                    {
                        case 1:
                        {
                            MARK;
                            PRINT_VARIABLE(elems[1]->GetIndex());
                            const unsigned node_a_local_index = p_face->GetNodeLocalIndex(node_a_index);
                            assert(node_a_local_index != UINT_MAX);
                            p_face->FaceAddNode(p_node_b, node_a_local_index);

                            p_elem->AddNode(p_node_b, p_elem->GetNumNodes() - 1);

                            break;
                        }
                        case 2:
                        {
                            MARK;
                            PRINT_VARIABLE(elems[2]->GetIndex());

                            p_face->FaceDeleteNode(p_node_b);
                            p_elem->DeleteNode(p_elem->GetNodeLocalIndex(p_node_b->GetIndex()));
                            p_elem->DeleteFace(p_elem->GetFaceLocalIndex(p_small_triangular_face->GetIndex()));

                            break;
                        }
                        case 3:
                        {
                            MARK;
                            PRINT_VARIABLE(elems[3]->GetIndex());
                            const unsigned node_b_local_index = p_face->GetNodeLocalIndex(node_b_index);
                            assert(node_b_local_index != UINT_MAX);
                            p_face->FaceAddNode(p_node_a, node_b_local_index);

                            p_elem->AddNode(p_node_a, p_elem->GetNumNodes() - 1);

                            break;
                        }
                        case 4:
                        {
                            MARK;
                            PRINT_VARIABLE(elems[4]->GetIndex());

                            p_face->FaceDeleteNode(p_node_a);
                            p_elem->DeleteNode(p_elem->GetNodeLocalIndex(p_node_a->GetIndex()));
                            p_elem->DeleteFace(p_elem->GetFaceLocalIndex(p_small_triangular_face->GetIndex()));

                            break;
                        }
                        default:
                            NEVER_REACHED;
                    }

                    MARK;
                    p_elem->MonolayerElementRearrangeFacesNodes();
                }

                // Now since we know the lateral swap face, we can move the nodes using the face normal.
                {
                    c_vector<double, 3> vector_CD = p_mesh->GetCellRearrangementRatio() * p_mesh->GetCellRearrangementThreshold() * vector_xy / norm_2(vector_xy);

                    // change sign accordingly
                    ///\todo #2850 condition based on element 2 and 4 locations
                    if (inner_prod(vector_CD, p_mesh->GetVectorFromAtoB(elems[2]->GetCentroid(), elems[4]->GetCentroid())) < 0)
                    {
                        vector_CD *= -1;
                    }

                    // Move nodes A and B to C and D respectively
                    p_node_a->rGetModifiableLocation() += 0.5 * vector_ab - 0.5 * vector_CD;
                    p_node_b->rGetModifiableLocation() += -0.5 * vector_ab + 0.5 * vector_CD;
                }

                p_mesh->DeleteFacePriorToReMesh(p_small_triangular_face->GetIndex());
                p_mesh->DeleteNodePriorToReMesh(p_lateral_node->GetIndex());
                // std::swap(p_mesh->mFaces[p_big_triangular_face->GetIndex()], p_mesh->mFaces[p_big_triangular_face->GetIndex()]);

            } // end of async
        }
    }
}


#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(LateralNodeModifier)
