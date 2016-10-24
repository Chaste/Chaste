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

#define OUTPUT_NAME "ReMesh33"

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
    std::string mName;
    unsigned mNumLowerNodes;
    std::vector<Node<3>*> mLowerNodes;
    std::vector<Node<3>*> mUpperNodes;
    // since number of nodes is known a priori, vector is good enough
    std::vector<std::vector<unsigned> > mNodeToLateralFaceIndices;
    std::vector<VertexElement<2, 3>*> mFaces;
    std::vector<VertexElement<3, 3>*> mElements;
    MutableVertexMesh<3, 3>* mpMesh;
    VertexMeshWriter<3, 3>* mpWriter;

public:
    MeshBuilderHelper(const std::vector<Node<3>*>& rLowerNodes, const std::string& Name, const unsigned zHeight = 1)
                        : mName(Name),
                          mNumLowerNodes(rLowerNodes.size()),
                          mLowerNodes(rLowerNodes),
                          mUpperNodes(mNumLowerNodes),
                          mNodeToLateralFaceIndices(mNumLowerNodes),
                          mFaces(),
                          mElements(),
                          mpMesh(NULL),
                          mpWriter(NULL)
    {
        // mUpperNodes uses copy constructor, need some updates
        for (unsigned i=0; i<mNumLowerNodes; ++i)
        {
            mLowerNodes[i]->AddNodeAttribute(1.1);

            const c_vector<double, 3> tmp = mLowerNodes[i]->rGetLocation();
            Node<3>* p_node_tmp = new Node<3>(i+mNumLowerNodes, mLowerNodes[i]->IsBoundaryNode(),
                    tmp[0], tmp[1], tmp[2] + zHeight);
            p_node_tmp->AddNodeAttribute(2.1);
            mUpperNodes[i] = p_node_tmp;
        }
    }

    MeshBuilderHelper(const std::string& Name)
    : mName(Name),
      mNumLowerNodes(0),
      mLowerNodes(),
      mUpperNodes(),
      mNodeToLateralFaceIndices(),
      mFaces(),
      mElements(),
      mpMesh(NULL),
      mpWriter(NULL)
    {}

    MutableVertexMesh<3, 3>* MakeMeshUsing2dMesh(const MutableVertexMesh<2, 2>& mesh2, const double zHeight=1)
    {
        mNumLowerNodes = mesh2.GetNumNodes();
        mLowerNodes.resize(mNumLowerNodes);
        mUpperNodes.resize(mNumLowerNodes);
        mNodeToLateralFaceIndices.resize(mNumLowerNodes);

        for (unsigned i=0 ; i<mNumLowerNodes ; ++i)
        {
            const Node<2>* p_2node = mesh2.GetNode(i);
            assert( i == p_2node->GetIndex() );
            const c_vector<double, 2> loc = p_2node->rGetLocation();
            const bool is_boundary = p_2node->IsBoundaryNode();
            Node<3>* p_lower = new Node<3>(i, is_boundary, loc[0], loc[1], 0);
            Node<3>* p_upper = new Node<3>(i+mNumLowerNodes, is_boundary, loc[0], loc[1], zHeight);
            p_lower->AddNodeAttribute(1.1);
            p_upper->AddNodeAttribute(2.1);
            mLowerNodes[i] = p_lower;
            mUpperNodes[i] = p_upper;
        }
        mElements.reserve(mesh2.GetNumElements());

        const unsigned num_elem = mesh2.GetNumElements();
        for (unsigned elem_index=0 ; elem_index<num_elem ; ++elem_index)
        {
            const VertexElement<2, 2>* p_2elem = mesh2.GetElement(elem_index);
            std::vector<unsigned> node_index_this_elem;
            for (unsigned i=0 ; i<p_2elem->GetNumNodes() ; ++i)
            {
                node_index_this_elem.push_back( p_2elem->GetNode(i)->GetIndex() );
            }
            this->buildElementWith(node_index_this_elem);
        }
        return this->GenerateMesh();
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

    void PrintMesh() const
    {
        for (unsigned i=0; i<mpMesh->GetNumElements(); ++i)
        {
            VertexElement<3,3>& elem = *(mpMesh->GetElement(i));
            std::cout << "ELEMENT (" << i<< ") : " << elem.GetIndex() << std::endl;
            std::cout << "number of Faces : " << elem.GetNumFaces() << " {  ";
            for (unsigned j=0; j<elem.GetNumFaces(); ++j)
            {
                std::cout << elem.GetFace(j)->GetIndex() << "  ";
            }
            std::cout << "}" << std::endl;
            std::cout << "number of Nodes : " << elem.GetNumNodes() << " {  ";
            for (unsigned j=0; j<elem.GetNumNodes(); ++j)
            {
                std::cout << elem.GetNode(j)->GetIndex() << "  ";
            }
            std::cout << "}" << std::endl;

            VertexElement<2,3>& basal = *(elem.GetFace(0));
            std::cout << "Nodes for basal face " << basal.GetIndex() << " {  ";
            for (unsigned j=0; j<basal.GetNumNodes(); ++j)
            {
                std::cout << basal.GetNode(j)->GetIndex() << "  ";
            }
            std::cout << "}" << std::endl;
        }
        std::cout << std::endl;
    }

    void WriteVtk(const std::string& AdditionalTag = "")
    {
        if (mpWriter == NULL)
        {
            mpWriter = new VertexMeshWriter<3, 3>(OUTPUT_NAME, mName, false);
        }
        else
        {
            // current workaround
            delete mpWriter;
            mpWriter = new VertexMeshWriter<3, 3>(OUTPUT_NAME, mName, false);
        }
        mpWriter->WriteVtkUsingMeshWithCellId(*mpMesh, AdditionalTag, false);
    }

    void buildElementWith(const unsigned numNodesThis, const unsigned nodeIndicesThis[] )
    {
        std::vector<unsigned> node_indices_this_elem(numNodesThis);
        for (unsigned id=0 ; id<numNodesThis ;  node_indices_this_elem[id] = nodeIndicesThis[id], ++id);

        buildElementWith(node_indices_this_elem);
    }

    void buildElementWith(const std::vector<unsigned>& nodeIndicesThisElem)
    {
        const unsigned num_nodes_this_elem = nodeIndicesThisElem.size();
        // Initializing vectors which are required for the generation of the VertexElement<3, 3>
        std::vector<VertexElement<2, 3>*> faces_this_elem;
        std::vector<bool> faces_orientation;

        std::vector<Node<3>*> lower_nodes_this_elem(num_nodes_this_elem);
        std::vector<Node<3>*> upper_nodes_this_elem(num_nodes_this_elem);
        std::vector<Node<3>*> all_nodes_this_elem(2*num_nodes_this_elem);
        // Populate lower & upper_nodes_this_elem
        for (unsigned j=0; j<num_nodes_this_elem; ++j)
        {
            lower_nodes_this_elem[j] = mLowerNodes[ nodeIndicesThisElem[j] ];
            upper_nodes_this_elem[j] = mUpperNodes[ nodeIndicesThisElem[j] ];
            all_nodes_this_elem[j] = mLowerNodes[ nodeIndicesThisElem[j] ];
            all_nodes_this_elem[j+num_nodes_this_elem] = mUpperNodes[ nodeIndicesThisElem[j] ];
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
        for (unsigned local_node_index=0; local_node_index<num_nodes_this_elem; ++local_node_index )
        {
            unsigned node1Index = nodeIndicesThisElem[local_node_index];
            unsigned node2Index = nodeIndicesThisElem[(local_node_index+1) % num_nodes_this_elem];

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
        VertexElement<3, 3>* p_elem = new VertexElement<3, 3>(mElements.size(), faces_this_elem, faces_orientation, all_nodes_this_elem);
        mElements.push_back( p_elem );
    }

    ~MeshBuilderHelper()
    {
        if (mpMesh)
            delete mpMesh;
        if (mpWriter)
            delete mpWriter;
    }
};

class TestMutableVertexMesh33ReMesh : public CxxTest::TestSuite
{
public:

    void TestCheckForSwapsAndIdentifySwapType() throw(Exception)
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
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {4, 1, 2, 5};
        unsigned node_indices_elem_2[3] = {0, 1, 4};
        unsigned node_indices_elem_3[4] = {4, 5, 3, 0};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "T1SwapWith4Elements", height);
        builder.buildElementWith(3, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        builder.buildElementWith(3, node_indices_elem_2);
        builder.buildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Set the threshold distance between vertices for a T1 swap as follows
        // so that it will trigger CheckForSwapsFromShortEdges
        vertex_mesh.SetCellRearrangementThreshold(0.3);
        vertex_mesh.CheckForSwapsFromShortEdges();
        builder.WriteVtk("AfterOnce");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.725, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.275, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.725, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.275, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1 , 1e-8);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = {2, 3, 5, 4};
        unsigned node_indices_element_1[3] = {4, 1, 2};
        unsigned node_indices_element_2[4] = {0, 1, 4, 5};
        unsigned node_indices_element_3[3] = {5, 3, 0};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);
        // Perform a T1 swap on nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        builder.WriteVtk("AfterTwice");

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0)+2*0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.2+0.2*sqrt(41.0)+2*0.3, 1e-6);

        // Test T1 swap location tracking
        std::vector< c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 4u);
        for (unsigned i=0 ; i<4 ; ++i)
        {
            TS_ASSERT_DELTA(t1_locations[i][0], 0.5, 1e-6);
            TS_ASSERT_DELTA(t1_locations[i][1], 0.5, 1e-6);
            TS_ASSERT_DELTA(t1_locations[i][2], i%2, 1e-6);
        }

        // Keep testing...
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 17u);

    }

    void TestPerformT1SwapOnBoundary() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all nodes are
         * boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *  _____
         * |\   /
         * | \ /
         * |  |
         * | / \
         * |/___\
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6));

        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[3] = {1, 4, 0};
        unsigned node_indices_elem_2[4] = {0, 4, 5, 3};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "T1SwapWithOnBoundary", height);
        builder.buildElementWith(3, node_indices_elem_0);
        builder.buildElementWith(3, node_indices_elem_1);
        builder.buildElementWith(4, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);
        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtk("After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0 , 1e-8);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 14u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = {2, 3, 5, 4};
        unsigned node_indices_element_1[4] = {1, 4, 5, 0};
        unsigned node_indices_element_2[4] = {0, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0)+2*0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0)+2*0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = !(i==5 || i==11);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT1SwapOnBoundary2() throw(Exception)
            {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all but one node
         * are boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *
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
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_1[3] = {1, 4, 0};
        unsigned node_indices_elem_2[4] = {0, 4, 5, 3};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "T1SwapWithOnBoundary2", height);
        builder.buildElementWith(4, node_indices_elem_0);
        builder.buildElementWith(3, node_indices_elem_1);
        builder.buildElementWith(4, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 14u);

        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtk("After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0 , 1e-8);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumFaces(), 5u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = {1, 2, 4};
        unsigned node_indices_element_1[4] = {1, 4, 5, 0};
        unsigned node_indices_element_2[3] = {0, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0)+2*0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapWhenVoidForms() throw(Exception)
    {
        /*
         * Create a mesh containing six nodes containing in two elements. We will test that
         * a T1 swap is correctly performed in the case where a void forms as a result of
         * the rearrangement, as shown below.
         *
         * |\   /|     |\      /|
         * | \ / |     | \    / |
         * |  |  |  => | /    \ |
         * | / \ |     |/      \|
         * |/   \|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[4] = {0, 4, 5, 3};
        unsigned node_indices_elem_1[4] = {4, 1, 2, 5};
        const double height = 1;
        MeshBuilderHelper builder(nodes, "T1SwapWhenVoidForms", height);
        builder.buildElementWith(4, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);
        // Perform a T1 swap on nodes 5 and 4.
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtk("After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0 , 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0 , 1e-8);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 10u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = {0, 5, 3};
        unsigned node_indices_element_1[3] = {4, 1, 2};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0)+2*0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapExceptions() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes containing in two triangle and two rhomboid elements,
         * where two nodes (those with indices 4 and 5) have the same location. We will test that
         * trying to perform a T1 swap on these nodes throws the correct exception.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.5, 0.0));

        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "T1SwapWith4Elements", height);
        builder.buildElementWith(3, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        builder.buildElementWith(3, node_indices_elem_2);
        builder.buildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Test that trying to perform a T1 swap on nodes 4 and 5 throws the correct exception
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5)), "Nodes are too close together, this shouldn't happen");
    }

    void TestDoNotPerforT1SwapWithRemovingEdgeFromTriangularElement() throw(Exception)
    {
        /**
         * In this test we check that a T1 swap does not occur if one of the elements is triangular
         * and would loose an edge by swapping nodes. The mesh looks like this
         *
         *       ______________
         *      |\             |
         *      | \ _________  |
         *      |  |          \| ...where the funny shaped element in the middle is supposed to be
         *      |  |_________ /|    a very long triangle that has the third vertex on the right hand boundary.
         *      | /            |
         *      |/_____________|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  2.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  2.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.3, 0.95, 0.0));
        nodes.push_back(new Node<3>(5, true, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(6, false, 0.3, 1.05, 0.0));

        unsigned node_indices_elem_0[4] = {0, 1, 5, 4};
        unsigned node_indices_elem_1[4] = {5, 2, 3, 6};
        unsigned node_indices_elem_2[4] = {0, 4, 6, 3};
        unsigned node_indices_elem_3[3] = { 4, 5, 6};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "NoT1SwapWithTriangularPrism", height);
        builder.buildElementWith(4, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        builder.buildElementWith(4, node_indices_elem_2);
        builder.buildElementWith(3, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);

        // Check for T1 swaps and carry them out if allowed - the short edge should not swap!
        vertex_mesh.CheckForSwapsFromShortEdges();
        builder.WriteVtk("After");

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 6u);

        // Test that each element still contains the correct nodes following the rearrangement
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_elem_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_elem_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_elem_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_elem_3[i]);
            }
        }
    }

    void TestExceptionForVoidRemovalWithRemovingEdgeFromTriangularElement() throw(Exception)
    {
        /**
         * In this test we check that void removal does not occur if one of the adjacent elements is triangular
         * and would loose an edge by swapping nodes. The code should throw and exception in this case.
         * The mesh looks like this
         *
         *       ______________./This corner is not a node.
         *      |\      1      |
         *      | \ _________  |
         *      |  |   void   \| ...where elements 1, and 2 are triangles that share the right hand vertex
         *      |  |_________ /|    with the triangular void in the middle.
         *      | /     2      |
         *      |/_____________|.This corner is not a node either.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.3, 0.95, 0.0));
        nodes.push_back(new Node<3>(3, true, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.3, 1.05, 0.0));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[4] = {0, 2, 4, 1};
        unsigned node_indices_elem_1[3] = {1, 4, 3};
        unsigned node_indices_elem_2[3] = {0, 3, 2};

        const double height = 1;
        MeshBuilderHelper builder(nodes, "NoT1SwapWithTriangularVoid", height);
        builder.buildElementWith(4, node_indices_elem_0);
        builder.buildElementWith(3, node_indices_elem_1);
        builder.buildElementWith(3, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Before");

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);
        builder.WriteVtk("After");

        // Check for possible swaps and carry them out if allowed - the short edge should not swap and
        // the void should not be removed!
        TS_ASSERT_THROWS_THIS(vertex_mesh.CheckForSwapsFromShortEdges(),
                "Triangular element next to triangular void, not implemented yet.");
    }

    void TestReMeshForT1Swaps() throw(Exception)
    {
        /*
         * Read in a vertex mesh that contains several pairs of nodes that are close enough for
         * T1 swaps to be performed, as shown below. The mesh consists of six elements and all
         * T1 swaps are performed on all horizontal edges. We will test that the ReMesh() method
         * correctly performs T1 swaps for internal and boundary elements, and correctly updates
         * which nodes are labelled as boundary nodes.
         *
         *      /\    /\
         *     /  \__/  \
         *    /   /  \   \
         *    \__/\__/\__/
         *    /  \/  \/  \
         *    \   \__/   /
         *     \  /  \  /
         *      \/    \/
         */
        VertexMeshReader<2,2> mesh_reader("cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1");
        MutableVertexMesh<2,2> vertex_2mesh;
        vertex_2mesh.ConstructFromMeshReader(mesh_reader);
        vertex_2mesh.SetCellRearrangementThreshold(0.1);

        MeshBuilderHelper builder("TestReMesh");
        MutableVertexMesh<3, 3>& vertex_mesh = *(builder.MakeMeshUsing2dMesh(vertex_2mesh) );
        builder.WriteVtk("");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 44u);

        // Calls ReMesh() to identify and perform any T1 swaps
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 44u);

        std::string dirname = "TestVertexMeshReMesh";
        std::string mesh_filename = "vertex_remesh_T1";

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
//        mesh_writer.WriteFilesUsingMesh(vertex_mesh);
//
//        // Check the positions are updated correctly
//        OutputFileHandler handler("TestVertexMeshReMesh", false);
//        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.node";
//        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.cell";
//
//        FileComparison comparer1(results_file1, "cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.node");
//        TS_ASSERT(comparer1.CompareFiles());
//        FileComparison comparer2(results_file2, "cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.cell");
//        TS_ASSERT(comparer2.CompareFiles());
    }
};

#endif /*TESTMUTABLEVERTEXMESH33REMESH_HPP_*/
