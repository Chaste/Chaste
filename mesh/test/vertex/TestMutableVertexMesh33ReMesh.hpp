#ifndef TESTMUTABLEVERTEXMESH33REMESH_HPP_
#define TESTMUTABLEVERTEXMESH33REMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "Debug.hpp"
#include "VertexMeshWriter.hpp"
#include "FileComparison.hpp"
#include "Warnings.hpp"
#include "MutableVertexMesh.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

// I need some helper function here. otherwise just the code to generate mesh would be too much!!
/**
 * nodeIndicesThisElem  array of global node indices which belongs to this to-be-created element in CCW
 * rLowerNodes  all the lower nodes
 * rUpperNodes  maybe all the upper nodes, if nothing here, it will be populated
 * rExistingFaces  all the created lateral faces so that no repeating faces, non-const as newly created face will be push back
 */
// so in the end a helper class is easier lol
class MeshBuilderHelper
{
private:
    unsigned mNumLowerNodes;
    std::vector<Node<3>*> mLowerNodes;
    std::vector<Node<3>*> mUpperNodes;
    // since number of nodes is known a priori, vector is good enough
    std::vector<std::vector<unsigned> > mNodeToLateralFaceIndices;
    std::vector<VertexElement<2, 3>*> mFaces;
    std::vector<VertexElement<3, 3>*> mElements;
    MutableVertexMesh<3, 3>* mpMesh;

public:
    MeshBuilderHelper(std::vector<Node<3>*>& rLowerNodes, unsigned zHeight = 1)
                :mNumLowerNodes(rLowerNodes.size()),
                 mLowerNodes(rLowerNodes),
                 mpMesh(NULL)
    {
        // mUpperNodes uses copy constructor, need some updates
        for (unsigned i=0; i<mNumLowerNodes; ++i)
        {
            mLowerNodes[i]->AddNodeAttribute(1.1);

            c_vector<double, 3> tmp = mLowerNodes[i]->rGetLocation();
            Node<3>* p_node_tmp = new Node<3>(i+mNumLowerNodes, true, tmp[0], tmp[1], tmp[2] + zHeight);
            p_node_tmp->AddNodeAttribute(2.1);
            mUpperNodes.push_back(p_node_tmp);
        }
        mNodeToLateralFaceIndices.reserve(mNumLowerNodes);

        // just checking, not sure if they copy the pointer address or actually copy the real nodes
        assert(mLowerNodes[0]->GetIndex()==0);
    }


    MutableVertexMesh<3, 3>* GenerateMesh()
    {
        // Combine the upper and lower nodes by adding upper_nodes into lower_nodes.
        // The index of upper nodes need to be modified.
        // unsigned lowerNodeLength = lower_nodes.size();
        unsigned upperNodeLength = mUpperNodes.size();

        for (unsigned upperRunningIndex=0; upperRunningIndex<upperNodeLength; ++upperRunningIndex)
        {
            mLowerNodes.push_back(mUpperNodes[upperRunningIndex]);
        }
        mpMesh = new MutableVertexMesh<3, 3>(mLowerNodes, mElements);
        return mpMesh;
    }


    void buildElementWith(const unsigned numNodesThis, const unsigned nodeIndicesThis[] )
    {
        // Initializing vectors which are required for the generation of the VertexElement<3, 3>
        std::vector<VertexElement<2, 3>*> faces_this_elem;
        std::vector<bool> faces_orientation;
//        unsigned nums = sizeof(nodeIndicesThis)/sizeof(*nodeIndicesThis);
//PRINT_3_VARIABLES(sizeof(nodeIndicesThis), sizeof(*nodeIndicesThis), nums);

        std::vector<Node<3>*> lower_nodes_this_elem(numNodesThis);
        std::vector<Node<3>*> upper_nodes_this_elem(numNodesThis);
        // Populate lower&upper_nodes_this_elem
        for (unsigned j=0; j<numNodesThis; ++j)
        {
            lower_nodes_this_elem[j] = mLowerNodes[ nodeIndicesThis[j] ];
            upper_nodes_this_elem[j] = mUpperNodes[ nodeIndicesThis[j] ];
        }


        // Creating the lower face
        VertexElement<2, 3>* p_lower_face = new VertexElement<2, 3>(mFaces.size(), lower_nodes_this_elem);
        // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
        // 1.1 instead of 1.0 as it will be casted into unsigned for simpler comparison.
        p_lower_face->AddElementAttribute(1.1);
        mFaces.push_back(p_lower_face);
        faces_this_elem.push_back(p_lower_face);
        faces_orientation.push_back(true);

        // Creating the upper face
        VertexElement<2,3>* p_upper_face = new VertexElement<2,3>(mFaces.size(), upper_nodes_this_elem);
        // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
        // 2.1 instead of 2.0 as it will be casted into unsigned for simpler comparison.
        p_upper_face->AddElementAttribute(2.1);
        mFaces.push_back(p_upper_face);
        faces_this_elem.push_back(p_upper_face);
        faces_orientation.push_back(false);

        // Creating all the lateral faces in CCW
        for (unsigned local_node_index=0; local_node_index<numNodesThis; ++local_node_index )
        {
            unsigned node1Index = nodeIndicesThis[local_node_index];
            unsigned node2Index = nodeIndicesThis[(local_node_index+1) % numNodesThis];

            // The values of the maps are called here because they will be used both existing and creating branch.
            // They are called by reference as they will be modified if they enter creating branch.
            std::vector<unsigned>& r_face1_indices = mNodeToLateralFaceIndices[node1Index];
            std::vector<unsigned>& r_face2_indices = mNodeToLateralFaceIndices[node2Index];

            // If both nodes already exist, the lateral face MIGHT have been created.
            unsigned existing_face_index = UINT_MAX;

            // Now need to search for the same lateral face index in both vector
            // not a too complicated and resource intensive (as r_faces_index vectors have length of at most 3 or 4
            // therefore not using existing function in <algorithm>
            for (unsigned i1 = 0; i1 < r_face1_indices.size() && existing_face_index==UINT_MAX; ++i1)
            {
                for (unsigned i2 = 0; i2 < r_face2_indices.size() && existing_face_index==UINT_MAX; ++i2)
                {
                    if (r_face1_indices[i1] == r_face2_indices[i2])
                    {
                        existing_face_index = r_face1_indices[i1];
                        break;
                    }
                }
            }

            if (existing_face_index != UINT_MAX) // meaning it's found
            {
                faces_this_elem.push_back(mFaces[existing_face_index]);
                // Face orientation is false as it was created by another element. CCW for another will be CW when
                // viewing from the other side as rotation is pseudovectorial
                faces_orientation.push_back(false);
            }


            if (existing_face_index == UINT_MAX)
            {
                // Create new lateral rectangular face
                std::vector<Node<3>*> nodes_of_lateral_face;
                nodes_of_lateral_face.push_back(mLowerNodes[node1Index]);
                nodes_of_lateral_face.push_back(mLowerNodes[node2Index]);
                nodes_of_lateral_face.push_back(mUpperNodes[node2Index]);
                nodes_of_lateral_face.push_back(mUpperNodes[node1Index]);

                unsigned newFaceIndex = mFaces.size();
                VertexElement<2, 3>* p_lateral_face = new VertexElement<2, 3>(newFaceIndex, nodes_of_lateral_face);
                // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
                // 3.1 instead of 3.0 as it will be casted into unsigned for simpler comparison.
                p_lateral_face->AddElementAttribute(3.1);
                mFaces.push_back(p_lateral_face);
                faces_this_elem.push_back(p_lateral_face);
                faces_orientation.push_back(true);

                // Update node_to_lateral_face_indices
                r_face1_indices.push_back(newFaceIndex);
                r_face2_indices.push_back(newFaceIndex);
            }
        }

        VertexElement<3, 3>* p_elem = new VertexElement<3, 3>(mElements.size(), faces_this_elem, faces_orientation);
        mElements.push_back( p_elem );
    }

    ~MeshBuilderHelper()
    {
        std::cout << "Destructorrr!!! \n";
        if (mpMesh)
        {
            delete mpMesh;
            std::cout << "mesh deleted!\n";
        }

    }
};

class TestMutableVertexMeshReMesh : public CxxTest::TestSuite
{
public:

    void TestPerformT1SwapAndIdentifySwapType() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};


        MeshBuilderHelper builder(nodes);
        builder.buildElementWith(3, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        builder.buildElementWith(3, node_indices_elem_2);
        builder.buildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        VertexMeshWriter<3, 3> writer("ReMesh33", "First_T1", false);
        writer.WriteVtkUsingMeshWithCellId(vertex_mesh, "", false);




        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);
MARK
        // Perform a T1 swap on nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(10), vertex_mesh.GetNode(11));
MARK; TRACE("Finish IdentifySwapType")

VertexMeshWriter<3, 3> writer2("ReMesh33", "First_T1_swap", false);
        writer2.WriteVtkUsingMeshWithCellId(vertex_mesh, "", false);

        // Test that each moved node has the correct location following the rearrangement
//        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
//        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
//        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
//        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
//
//        // Test that each element contains the correct nodes following the rearrangement
//        unsigned node_indices_element_0[4] = {2, 3, 5, 4};
//        unsigned node_indices_element_1[3] = {2, 4, 1};
//        unsigned node_indices_element_2[4] = {1, 4, 5, 0};
//        unsigned node_indices_element_3[3] = {0, 5, 3};
//        for (unsigned i=0; i<4; i++)
//        {
//            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
//            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
//            if (i < 3)
//            {
//                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
//                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
//            }
//        }



















    }
};

#endif /*TESTMUTABLEVERTEXMESH33REMESH_HPP_*/
