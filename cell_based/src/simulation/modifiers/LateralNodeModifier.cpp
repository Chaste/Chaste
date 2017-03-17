/*
 * LateralNodeModifier.cpp
 *
 *  Created on: 11 Mar 2017
 *      Author: Weijie
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

    if (dynamic_cast<MutableVertexMesh<3, 3>*>(p_tmp_mesh) == NULL)
    {
        EXCEPTION("only for vertex mesh");
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

        for (std::set<VertexElement<2, 3>*>::const_iterator it = lateral_faces.begin();
             it != lateral_faces.end(); ++it)
        {
            VertexElement<2, 3>* p_face_tmp = *it;
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

        const c_vector<double, 3> mid_apical = (apical_nodes[0]->rGetLocation() + apical_nodes[1]->rGetLocation()) / 2;
        const c_vector<double, 3> mid_basal = (basal_nodes[0]->rGetLocation() + basal_nodes[1]->rGetLocation()) / 2;

        const double apical_length = p_mesh->GetDistanceBetweenNodes(apical_nodes[0]->GetIndex(), apical_nodes[1]->GetIndex());
        const double basal_length = p_mesh->GetDistanceBetweenNodes(basal_nodes[0]->GetIndex(), basal_nodes[1]->GetIndex());

        // HOW_MANY_TIMES_HERE("bla bla bla");
        // PRINT_CONTAINER(p_lateral_node->rGetLocation());
        p_lateral_node->rGetModifiableLocation() = (basal_length * mid_apical + apical_length * mid_basal) / (apical_length + basal_length);
        // PRINT_CONTAINER(p_lateral_node->rGetLocation());

        // Reverse AsynchronousT1 when the other edge is shorter than mCellRearrangementThreshold
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

            std::vector<VertexElement<3, 3>*> elems(5, NULL);
            {
                const std::set<unsigned>& elem_24 = p_small_triangular_face->rFaceGetContainingElementIndices();
                const std::set<unsigned>& elem_a_124 = p_node_a->rGetContainingElementIndices();
                const std::set<unsigned>& elem_b_234 = p_node_b->rGetContainingElementIndices();
                MARK;
                PRINT_CONTAINER(elem_24);
                PRINT_CONTAINER(elem_a_124);
                PRINT_CONTAINER(elem_b_234);

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
                const c_vector<double, 3> centroid_4 = (elems[4] != NULL) ? GetFacesWithType(elems[4], T1_type)[0]->GetCentroid()
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
                for (std::set<unsigned>::const_iterator face_it = faces_tmp.begin(); face_it != faces_tmp.end(); ++face_it)
                {
                    p_mesh->GetFace(*face_it)->FaceDeleteNode(p_lateral_node);
                }

                const std::set<unsigned> elems_tmp = p_lateral_node->rGetContainingElementIndices();
                for (std::set<unsigned>::const_iterator elem_it = elems_tmp.begin(); elem_it != elems_tmp.end(); ++elem_it)
                {
                    p_mesh->GetElement(*elem_it)->DeleteNode(p_lateral_node);
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
                VertexElement<2, 3>* p_lateral_face_14 = NULL;
                const std::set<unsigned>& tmp_face_ids = p_node_a->rGetContainingFaceIndices();
                if (elems[4] != NULL)
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
                if (elems[i] == NULL)
                {
                    continue;
                }

                VertexElement<3, 3>* p_elem = elems[i];
                VertexElement<2, 3>* p_face = GetFacesWithType(p_elem, T1_type)[0];

                switch (i)
                {
                    case 1:
                    {
                        MARK;
                        PRINT_VARIABLE(elems[1]->GetIndex());

                        p_face->FaceAddNode(p_node_b, p_face->GetNumNodes() - 1);
                        FaceRearrangeNodesInMesh(p_mesh, p_face);
                        p_elem->AddNode(p_node_b, p_elem->GetNumNodes() - 1);
                        break;
                    }
                    case 2:
                    {
                        MARK;
                        PRINT_VARIABLE(elems[2]->GetIndex());

                        p_face->FaceDeleteNode(p_node_b);
                        p_elem->DeleteNode(p_node_b);
                        p_elem->DeleteFace(p_small_triangular_face);
                        break;
                    }
                    case 3:
                    {
                        MARK;
                        PRINT_VARIABLE(elems[3]->GetIndex());

                        p_face->FaceAddNode(p_node_a, p_face->GetNumNodes() - 1);
                        FaceRearrangeNodesInMesh(p_mesh, p_face);
                        p_elem->AddNode(p_node_a, p_elem->GetNumNodes() - 1);
                        break;
                    }
                    case 4:
                    {
                        MARK;
                        PRINT_VARIABLE(elems[4]->GetIndex());

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
