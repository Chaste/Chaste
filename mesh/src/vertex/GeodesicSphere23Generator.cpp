/*
 * GeodesicSphere23Generator.cpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Weijie
 */

#include "GeodesicSphere23Generator.hpp"

#include <utility>              // for pair
#include <set>
#include <algorithm>            // sort
#include "UblasCustomFunctions.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include <cmath>                                    // for atan2

#include "Debug.hpp"


class SortWithIndex
{
public:
    bool operator()(const std::pair<double, unsigned>& a, const std::pair<double, unsigned>& b) const
    {
        return a.first < b.first;
    }
};

GeodesicSphere23Generator::GeodesicSphere23Generator(const unsigned numDivision)
: mDepth(0)
{
    const double X = 0.525731112119133606;
    const double Z = 0.850650808352039932;

    const double vdata[12][3] = {
            {-X, 0.0, Z}, { X, 0.0, Z }, { -X, 0.0, -Z }, { X, 0.0, -Z },
            { 0.0, Z, X }, { 0.0, Z, -X }, { 0.0, -Z, X }, { 0.0, -Z, -X },
            { Z, X, 0.0 }, { -Z, X, 0.0 }, { Z, -X, 0.0 }, { -Z, -X, 0.0 }};
    for (unsigned i=0; i<12; ++i)
    {
        mNodes.push_back(new Node<3>(i, false, vdata[i][0], vdata[i][1], vdata[i][2]));
    }

    int tindices[20][3] = {
            {0, 4, 1}, { 0, 9, 4 }, { 9, 5, 4 }, { 4, 5, 8 }, { 4, 8, 1 },
            { 8, 10, 1 }, { 8, 3, 10 }, { 5, 3, 8 }, { 5, 2, 3 }, { 2, 7, 3 },
            { 7, 10, 3 }, { 7, 6, 10 }, { 7, 11, 6 }, { 11, 0, 6 }, { 0, 1, 6 },
            { 6, 1, 10 }, { 9, 0, 11 }, { 9, 11, 2 }, { 9, 2, 5 }, { 7, 2, 11 }};

    // Creating each face (edges at the same time.)
    for (unsigned i = 0; i<20; ++i)
    {
        std::vector<Node<3>*> this_face_nodes;
        std::vector<VertexElement<1,3>*> this_face_edges;
        // well, orientation doesn't really matter.
        std::vector<bool> this_face_orientations(3);

        for (unsigned ja=0; ja<3; ++ja)
        {
            this_face_nodes.push_back(mNodes[tindices[i][ja]]);

            const unsigned jb = (ja+1)%3;
            Node<3>* node_a = mNodes[tindices[i][ja]];
            Node<3>* node_b = mNodes[tindices[i][jb]];
            //PRINT_3_VARIABLES(i, tindices[i][ja], tindices[i][jb])
            //PRINT_CONTAINER(node_a_edges)
            //PRINT_CONTAINER(node_b_edges)
            std::set<unsigned> shared_edge = GetSharedFaceIndices(node_a, node_b);

            unsigned current_edge_index = UINT_MAX;
            if (shared_edge.size() == 0)
            {
                std::vector<Node<3>*> this_edge_nodes;
                this_edge_nodes.push_back(node_a);
                this_edge_nodes.push_back(node_b);
                current_edge_index = mEdges.size();
                mEdges.push_back(new VertexElement<1,3>(current_edge_index, this_edge_nodes));
            }
            else
            {
                //TRACE("FOUND shared edge")
                //PRINT_3_VARIABLES(i, tindices[i][ja], tindices[i][jb])
                //PRINT_CONTAINER(node_a_edges)
                //PRINT_CONTAINER(node_b_edges)
                assert(shared_edge.size() == 1);
                current_edge_index = *(shared_edge.begin());
            }

            this_face_edges.push_back(mEdges[current_edge_index]);
        }
        assert(mFaces.size() == i);
        mFaces.push_back(new VertexElement<2,3>(mFaces.size(), this_face_edges, this_face_orientations, this_face_nodes));
    }

    PRINT_3_VARIABLES(mNodes.size(), mEdges.size(), mFaces.size());

    for (unsigned i=0; i<numDivision; ++i)
    {
        SubDivide();
        PRINT_3_VARIABLES(mNodes.size(), mEdges.size(), mFaces.size());
    }
}

void GeodesicSphere23Generator::SubDivide()
{
    const unsigned num_old_nodes = mNodes.size();
    const unsigned num_old_edges = mEdges.size();
    const unsigned num_old_faces = mFaces.size();
    if (num_old_nodes+num_old_faces-num_old_edges != 2)
    {
        TRACE("EXCEPTION EULER IS THE BOSS")
                                    EXCEPTION("Breaking Euler's Polyhedron Formula!");
    }

    // For sanity check later
    const unsigned num_new_total_nodes = num_old_nodes + num_old_edges;
    const unsigned num_new_edges = 2*num_old_edges + 3*num_old_faces;
    const unsigned num_new_faces = 4*num_old_faces;

    // Each edge will be "bisected", and each face will become 4 new faces.
    std::vector<Node<3>*> map_edge_to_new_nodes(mEdges.size());
    for (unsigned i=0; i<num_old_edges; ++i)
    {
        Node<3>* p_tmp_node = new Node<3>(mNodes.size(), normalise(mEdges[i]->GetCentroid()), false);
        map_edge_to_new_nodes[i] = p_tmp_node;
        mNodes.push_back(p_tmp_node);
    }
    if (mNodes.size() != num_new_total_nodes) 
    {
        NEVER_REACHED;
    }

    const unsigned cyc[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
    // Creating new faces and edges
    for (unsigned old_face_index=0; old_face_index<num_old_faces; ++old_face_index)
    {
        /*
         *           1
         *           /\
         *          /  \
         *         /  b \
         *      0'/______\2'
         *       /\      /\
         *      /  \  d /  \
         *     /  c \  /  a \
         *    /______\/______\
         *   2        1'      0
         */
        VertexElement<2, 3>* p_face = mFaces[old_face_index];
        if (old_face_index != p_face->GetIndex()) 
        {
            NEVER_REACHED;
        }

        Node<3>* old_nodes[3] = {p_face->GetNode(0), p_face->GetNode(1), p_face->GetNode(2)};
        std::vector<Node<3>*> mid_nodes(3);
        for (unsigned i=0; i<3; ++i)
        {
            const std::set<unsigned> shared_edge =
                    GetSharedFaceIndices(old_nodes[cyc[i][0]], old_nodes[cyc[i][1]]);
            if (shared_edge.size() != 1) 
            {
                NEVER_REACHED;
            }
            if(*(shared_edge.begin()) != p_face->GetFace(i)->GetIndex()) 
            {
                NEVER_REACHED;
            }
            mid_nodes[cyc[i][2]] = map_edge_to_new_nodes[*(shared_edge.begin())];
        }

        // create face a,b,c
        for (unsigned i=0; i<4; ++i)
        {
            std::vector<Node<3>*> this_face_nodes(3);
            if (i<3)
            {
                this_face_nodes[0] = old_nodes[cyc[i][0]];
                this_face_nodes[1] = mid_nodes[cyc[i][1]];
                this_face_nodes[2] = mid_nodes[cyc[i][2]];
            }
            else
            {
                this_face_nodes[0] = mid_nodes[0];
                this_face_nodes[1] = mid_nodes[1];
                this_face_nodes[2] = mid_nodes[2];
            }


            std::vector<VertexElement<1,3>*> this_face_edges(3);
            std::vector<bool> this_orientation(3);
            for (unsigned j=0; j<3; ++j)
            {
                const std::set<unsigned> shared_edge = GetSharedFaceIndices(this_face_nodes[cyc[j][0]], this_face_nodes[cyc[j][1]]);
                VertexElement<1,3>* p_tmp_edge;
                if (shared_edge.size() == 0)
                {
                    std::vector<Node<3>*> tmp_edge_nodes;
                    tmp_edge_nodes.push_back(this_face_nodes[cyc[j][0]]);
                    tmp_edge_nodes.push_back(this_face_nodes[cyc[j][1]]);

                    p_tmp_edge = new VertexElement<1,3>(mEdges.size(), tmp_edge_nodes);
                    p_tmp_edge->RegisterFaceWithNodes();
                    mEdges.push_back(p_tmp_edge);
                }
                else
                {
                    if (shared_edge.size() != 1) 
                    {
                        NEVER_REACHED;
                    }
                    p_tmp_edge = mEdges[*(shared_edge.begin())];
                }
                this_face_edges[j] = p_tmp_edge;
            }

            if (i<3)
            {
                mFaces.push_back(new VertexElement<2,3>(mFaces.size(), this_face_edges,
                        this_orientation, this_face_nodes));
            }
            else
            {
                p_face->MarkAsDeleted();
                delete p_face;
                mFaces[old_face_index] = new VertexElement<2,3>(old_face_index, this_face_edges,
                        this_orientation, this_face_nodes);
            }
        }
    }


    assert(mEdges.size() == num_old_edges+num_new_edges);
    for (unsigned i=0; i<num_old_edges; ++i)
    {
        mEdges[i]->MarkFaceAsDeleted();
        delete mEdges[i];

        mEdges[i] = mEdges[num_new_edges+i];
        mEdges[i]->FaceResetIndex(i);
    }

    for (unsigned i=num_new_edges, j=0; i<mEdges.size(); ++i, ++j)
    {
        if (mEdges[i]->GetIndex() != j) 
        {
            NEVER_REACHED;
        }
    }
    mEdges.resize(num_new_edges);

    if (mFaces.size() != num_new_faces) 
    {
        NEVER_REACHED;
    }
    MARK;
}

MutableVertexMesh<2, 3>* GeodesicSphere23Generator::GetDual()
{
    std::vector<Node<3>*> dual_nodes(mFaces.size());
    for (unsigned dual_node_index=0; dual_node_index<mFaces.size(); ++dual_node_index)
    {
        dual_nodes[dual_node_index] = new Node<3>(dual_node_index, mFaces[dual_node_index]->GetCentroid(), false);
    }

    std::vector<VertexElement<2, 3>*> dual_faces(mNodes.size());
    for (unsigned dual_face_index=0; dual_face_index<mNodes.size(); ++dual_face_index)
    {
        const c_vector<double, 3> normal_v = normalise(mNodes[dual_face_index]->rGetLocation());

        const std::set<unsigned> s_tmp = mNodes[dual_face_index]->rGetContainingElementIndices();
        // Switch to vector for better access
        const std::vector<unsigned> this_dual_node_global_indices(s_tmp.begin(), s_tmp.end());
        const unsigned this_num_dual_nodes = this_dual_node_global_indices.size();


        std::vector<Node<3>*> unsorted_this_dual_face_nodes(this_num_dual_nodes);
        std::vector<std::pair<double, unsigned> > angles_with_index;
        unsorted_this_dual_face_nodes[0] = dual_nodes[this_dual_node_global_indices[0]];
        //            angles_with_index.push_back(std::make_pair(0.0, 0));

        // Create orthonormal basis with the face normal, first dual node location of the face.
        // e1 is obtained by Gram-Schmidt process.
        const c_vector<double, 3> e1 = normalise(unsorted_this_dual_face_nodes[0]->rGetLocation()
                - inner_prod(unsorted_this_dual_face_nodes[0]->rGetLocation(), normal_v)*normal_v);
        // e2 simply by cross prod of normal_v and e1
        const c_vector<double, 3> e2 = VectorProduct(normal_v, e1);
        assert(abs(norm_2(e2)-1) < 1e-5);

        for (unsigned dual_node_local_index=0; dual_node_local_index<this_num_dual_nodes; ++dual_node_local_index)
        {
            Node<3>* p_node = dual_nodes[this_dual_node_global_indices[dual_node_local_index]];
            unsorted_this_dual_face_nodes[dual_node_local_index] = p_node;

            const c_vector<double, 3> r_node = p_node->rGetLocation();
            const double tmp_angle = atan2(inner_prod(r_node, e1), inner_prod(r_node, e2));
            angles_with_index.push_back(std::make_pair(tmp_angle, dual_node_local_index));
        }

        std::sort(angles_with_index.begin(), angles_with_index.end(), SortWithIndex());
        std::vector<Node<3>*> this_dual_face_nodes(this_num_dual_nodes);

        std::vector<double> Debug_anglesss;
        for (unsigned iii=0; iii<this_num_dual_nodes; ++iii)
        {
            this_dual_face_nodes[iii] = unsorted_this_dual_face_nodes[angles_with_index[iii].second];

            Debug_anglesss.push_back(angles_with_index[iii].first);
        }
        //PRINT_VECTOR(Debug_anglesss);

        dual_faces[dual_face_index] = new VertexElement<2, 3>(dual_face_index, this_dual_face_nodes);
    }


    return new MutableVertexMesh<2,3>(dual_nodes, dual_faces);
}






