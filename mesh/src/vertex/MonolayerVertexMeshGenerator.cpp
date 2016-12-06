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

#include "MonolayerVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include <algorithm>

MonolayerVertexMeshGenerator::MonolayerVertexMeshGenerator(const std::string& name)
    : mName(name),
      mpMesh(NULL),
      mpWriter(NULL)
{
}

MonolayerVertexMeshGenerator::MonolayerVertexMeshGenerator(const std::vector<Node<3>*>& rLowerNodes,
                                                           const std::string& name,
                                                           const unsigned zHeight)
    : mName(name),
      mBasalNodes(rLowerNodes),
      mApicalNodes(mBasalNodes.size()),
      mNodeLateralFaceMap(mBasalNodes.size()),
      mpMesh(NULL),
      mpWriter(NULL)
{
    const unsigned num_lower_nodes = mBasalNodes.size();
    for (unsigned node_index=0; node_index<num_lower_nodes; ++node_index)
    {
        mBasalNodes[node_index]->AddNodeAttribute(1.1);

        // Generate new apical nodes for each basal node
        c_vector<double, 3> tmp_location;
        tmp_location = mBasalNodes[node_index]->rGetLocation();
        Node<3>* p_node_tmp = new Node<3>(node_index+num_lower_nodes, mBasalNodes[node_index]->IsBoundaryNode(),
                                          tmp_location[0], tmp_location[1], tmp_location[2] + zHeight);
        p_node_tmp->AddNodeAttribute(2.1);
        mApicalNodes[node_index] = p_node_tmp;
    }
}

MonolayerVertexMeshGenerator::~MonolayerVertexMeshGenerator()
{
    if (mpMesh)
    {
        delete mpMesh;
    }
    if (mpWriter)
    {
        delete mpWriter;
    }
}

MutableVertexMesh<3, 3>* MonolayerVertexMeshGenerator::MakeMeshUsing2dMesh(const MutableVertexMesh<2, 2>& mesh2d,
                                                                           const double zHeight)
{
    const unsigned num_lower_nodes = mesh2d.GetNumNodes();
    mBasalNodes.resize(num_lower_nodes);
    mApicalNodes.resize(num_lower_nodes);
    mNodeLateralFaceMap.resize(num_lower_nodes);

    for (unsigned i=0 ; i<num_lower_nodes ; ++i)
    {
        const Node<2>* p_2node = mesh2d.GetNode(i);
        assert(i == p_2node->GetIndex());
        c_vector<double, 2> loc;
        loc = p_2node->rGetLocation();
        const bool is_boundary = p_2node->IsBoundaryNode();
        Node<3>* p_lower = new Node<3>(i, is_boundary, loc[0], loc[1], 0);
        Node<3>* p_upper = new Node<3>(i+num_lower_nodes, is_boundary, loc[0], loc[1], zHeight);
        p_lower->AddNodeAttribute(1.1);
        p_upper->AddNodeAttribute(2.1);
        mBasalNodes[i] = p_lower;
        mApicalNodes[i] = p_upper;
    }
    mElements.reserve(mesh2d.GetNumElements());

    const unsigned num_elem = mesh2d.GetNumElements();
    for (unsigned elem_index=0 ; elem_index<num_elem ; ++elem_index)
    {
        const VertexElement<2, 2>* p_2elem = mesh2d.GetElement(elem_index);
        std::vector<unsigned> node_index_this_elem;
        for (unsigned i=0 ; i<p_2elem->GetNumNodes() ; ++i)
        {
            node_index_this_elem.push_back( p_2elem->GetNode(i)->GetIndex() );
        }
        this->BuildElementWith(node_index_this_elem);
    }
    return this->GenerateMesh();
}

void MonolayerVertexMeshGenerator::BuildElementWith(const unsigned numBasalNodes,
                                                    const unsigned basalNodeIndices[])
{
    std::vector<unsigned> node_indices_this_elem(numBasalNodes);
    for (unsigned id=0 ; id<numBasalNodes ; ++id)
    {
        node_indices_this_elem[id] = basalNodeIndices[id];
    }

    BuildElementWith(node_indices_this_elem);
}

void MonolayerVertexMeshGenerator::BuildElementWith(const std::vector<unsigned>& basalNodeIndices)
{
    const unsigned num_nodes_this_elem = basalNodeIndices.size();
    // Initializing vectors which are required for the generation of the VertexElement<3, 3>
    std::vector<VertexElement<2, 3>*> faces_this_elem;
    std::vector<bool> faces_orientation;

    std::vector<Node<3>*> lower_nodes_this_elem(num_nodes_this_elem);
    std::vector<Node<3>*> upper_nodes_this_elem(num_nodes_this_elem);
    std::vector<Node<3>*> all_nodes_this_elem(2*num_nodes_this_elem);
    // Populate lower & upper_nodes_this_elem
    for (unsigned j=0; j<num_nodes_this_elem; ++j)
    {
        lower_nodes_this_elem[j] = mBasalNodes[ basalNodeIndices[j] ];
        upper_nodes_this_elem[j] = mApicalNodes[ basalNodeIndices[j] ];
        all_nodes_this_elem[j] = mBasalNodes[ basalNodeIndices[j] ];
        all_nodes_this_elem[j+num_nodes_this_elem] = mApicalNodes[ basalNodeIndices[j] ];
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
        unsigned node1Index = basalNodeIndices[local_node_index];
        unsigned node2Index = basalNodeIndices[(local_node_index+1) % num_nodes_this_elem];

        // The values of the maps are called here because they will be used both existing and creating branch.
        // They are called by reference as they will be modified if they enter creating branch.
        std::set<unsigned>& r_face1_indices = mNodeLateralFaceMap[node1Index];
        std::set<unsigned>& r_face2_indices = mNodeLateralFaceMap[node2Index];

        // If both nodes already exist, the lateral face MIGHT have been created.
        std::vector<unsigned> common_index;
        // Now need to search for the same lateral face index in both sets
        std::set_intersection(r_face1_indices.begin(), r_face1_indices.end(),
                              r_face2_indices.begin(), r_face2_indices.end(),
                              std::back_inserter(common_index));

        assert(common_index.size()<=1);
        unsigned existing_face_index = common_index.size()==0 ? UINT_MAX : common_index[0];

        if (existing_face_index != UINT_MAX) // meaning it's found
        {
            faces_this_elem.push_back(mFaces[existing_face_index]);
            // Face orientation is false as it was created by another element. CCW for another will be CW when
            // viewing from the other side as rotation is pseudovectorial
            faces_orientation.push_back(true);
        }

        if (existing_face_index == UINT_MAX)
        {
            // Create new lateral rectangular face
            std::vector<Node<3>*> nodes_of_lateral_face;
            nodes_of_lateral_face.push_back(mBasalNodes[node1Index]);
            nodes_of_lateral_face.push_back(mBasalNodes[node2Index]);
            nodes_of_lateral_face.push_back(mApicalNodes[node2Index]);
            nodes_of_lateral_face.push_back(mApicalNodes[node1Index]);

            unsigned newFaceIndex = mFaces.size();
            VertexElement<2, 3>* p_lateral_face = new VertexElement<2, 3>(newFaceIndex, nodes_of_lateral_face);
            // Attribute is added so that it can be identified in simulation (as basal, apical and lateral faces have different contributions)
            // 3.1 instead of 3.0 as it will be casted into unsigned for simpler comparison.
            p_lateral_face->AddElementAttribute(3.1);
            mFaces.push_back(p_lateral_face);
            faces_this_elem.push_back(p_lateral_face);
            faces_orientation.push_back(false);
            // Update node_to_lateral_face_indices
            r_face1_indices.insert(newFaceIndex);
            r_face2_indices.insert(newFaceIndex);
        }
    }
    VertexElement<3, 3>* p_elem = new VertexElement<3, 3>(mElements.size(), faces_this_elem, faces_orientation, all_nodes_this_elem);
    mElements.push_back(p_elem);
}

MutableVertexMesh<3, 3>* MonolayerVertexMeshGenerator::GenerateMesh()
{
    // Combine the upper and lower nodes by adding upper_nodes into lower_nodes.
    // The index of upper nodes need to be modified.
    // unsigned lowerNodeLength = lower_nodes.size();
    assert(mApicalNodes.size() == mBasalNodes.size());

    mBasalNodes.insert(this->mBasalNodes.end(), mApicalNodes.begin(), mApicalNodes.end());

    mpMesh = new MutableVertexMesh<3, 3>(mBasalNodes, mElements);
    return mpMesh;
}

MutableVertexMesh<3, 3>* MonolayerVertexMeshGenerator::ConvertMeshToCylinder(const double widthX,
        const double widthY, const double radius, const double thickness, const double length)
{
    const double mTol = 1e-6;
    // normal mesh should already exist before call of this function.
    if (mpMesh == NULL)
    {
        EXCEPTION("No mesh available. Please generate flat mesh before using this function.");
    }

    /**
     * Outline:
     * I    first find the (boundary) nodes which are located along the 'glue' line.
     *      Replace nodes and mark one as deleted.
     * II   Update the (boundary) faces along the 'glue' line and delete the faces as well.
     * III  Convert the coordinates to 'pseudo'-cylindrical and update their new coordinates.
     */

    // Operation I
    std::vector<Node<3>*> boundary_nodes;
    for (unsigned node_index=0; node_index<mpMesh->GetNumAllNodes(); ++node_index)
    {
        Node<3>* node_tmp = mpMesh->GetNode(node_index);
        if (!node_tmp->IsBoundaryNode())
        {
            continue;
        }
        else
        {
            boundary_nodes.push_back(node_tmp);
        }
    }

    // Store the nodes and faces on the go, but mark them as deleted only after all modifications
    // are completed as otherwise we might not be able to call the nodes from the mesh due to node mapping.
    std::vector<Node<3>*> deleted_nodes;
    std::vector<Node<3>*> glued_nodes;
    std::vector<VertexElement<2, 3>*> deleted_faces;

    // Iterate over the boundary nodes and find the pair of congruent nodes.
    for (unsigned index_a=0; index_a<boundary_nodes.size(); ++index_a)
    {
        Node<3>* node_a = boundary_nodes[index_a];
        // Skip over identified nodes for efficiency
        if (node_a == NULL)
        {
            continue;
        }
        const c_vector<double, 3> location_a = node_a->rGetLocation();

        for (unsigned index_b=index_a+1; index_b<boundary_nodes.size(); ++index_b)
        {
            Node<3>* node_b = boundary_nodes[index_b];
            if (node_b == NULL)
            {
                continue;
            }
            const c_vector<double, 3> location_b = node_b->rGetLocation();
            // Avoid warning on -Wmaybe-uninitialized
            c_vector<double, 3> a_to_b;
            a_to_b = location_b - location_a;

            // as they are double precision floating number, their absolute values are compared
            if (abs(abs(a_to_b[0]) - widthX)<mTol && abs(a_to_b[1])<mTol && abs(a_to_b[2])<mTol)
            {
                // for simplicity just delete the node b
                boundary_nodes[index_a] = NULL;
                boundary_nodes[index_b] = NULL;
                glued_nodes.push_back(node_a);
                deleted_nodes.push_back(node_b);

                const std::set<unsigned> elements_of_b = node_b->rGetContainingElementIndices();
                for (std::set<unsigned>::iterator it=elements_of_b.begin(); it!=elements_of_b.end(); ++it)
                {
                    VertexElement<3, 3>* elem_tmp = mpMesh->GetElement(*it);
                    elem_tmp->ReplaceNode(node_b, node_a);
                }

                // Break the loop as the congruent node is found.
                break;
            }
        }
    }

    // Operation II
    std::vector<VertexElement<2, 3>*> boundary_faces;
    std::vector<std::set<unsigned> > nodes_for_boundary_face;
    for (unsigned face_index=0; face_index<mpMesh->GetNumFaces(); ++face_index)
    {
        VertexElement<2, 3>* face_tmp = mpMesh->GetFace(face_index);
        if (!IsLateralFace(face_tmp) || !IsFaceOnBoundary(face_tmp))
        {
            continue;
        }
        else
        {
            assert(face_tmp->GetNumNodes()==4u);
            std::set<unsigned> tmp_set;
            // save only the basal nodes
            for (unsigned i=0; i<2u; ++i)
            {
                tmp_set.insert(face_tmp->GetNodeGlobalIndex(i));
            }
            boundary_faces.push_back(face_tmp);
            nodes_for_boundary_face.push_back(tmp_set);
        }
    }

    for (unsigned index_a=0; index_a<boundary_faces.size(); ++index_a)
    {
        VertexElement<2, 3>* p_face_a = boundary_faces[index_a];
        if (p_face_a == NULL)
        {
            continue;
        }

        for (unsigned index_b=index_a+1; index_b<boundary_faces.size(); ++index_b)
        {
            if (boundary_faces[index_b] == NULL)
            {
                continue;
            }

            VertexElement<2, 3>* p_face_b = boundary_faces[index_b];

            if (nodes_for_boundary_face[index_a] == nodes_for_boundary_face[index_b])
            {
                boundary_faces[index_a] = NULL;
                boundary_faces[index_b] = NULL;
                deleted_faces.push_back(p_face_b);

                std::set<unsigned> elem_a_index = p_face_a->rGetContainingElementIndices();
                assert(elem_a_index.size() == 1);
                VertexElement<3, 3>* p_elem_a = mpMesh->GetElement(*elem_a_index.begin());
                ///\todo: #2850 get face orientation
                std::vector<unsigned> v_tmp = GetLateralFace(p_elem_a, p_face_a->GetNodeGlobalIndex(0), p_face_a->GetNodeGlobalIndex(1));
                assert(v_tmp[0] != UINT_MAX);
                const bool face_orientation_a = v_tmp[1];

                std::set<unsigned> elem_b_index = p_face_b->rGetContainingElementIndices();
                assert(elem_b_index.size() == 1);
                VertexElement<3, 3>* p_elem_b = mpMesh->GetElement(*elem_b_index.begin());
                v_tmp = GetLateralFace(p_elem_b, p_face_b->GetNodeGlobalIndex(0), p_face_b->GetNodeGlobalIndex(1));
                assert(v_tmp[0] == p_face_b->GetIndex());
                const unsigned index_tmp = v_tmp[2];

                p_elem_b->ReplaceFace(index_tmp, p_face_a, !face_orientation_a);
                p_elem_b->MonolayerElementRearrangeFacesNodes();
                //MARK
                // Break the loop as the congruent node is found.
                break;
            }
        }
    }

    for (unsigned ind=0; ind<glued_nodes.size(); ++ind)
    {
        // Since it is still during its creation, there should be no higher order junctions
        if (glued_nodes[ind]->GetNumContainingElements() > 2)
        {
            glued_nodes[ind]->SetAsBoundaryNode(false);
        }
    }

    for (unsigned ind=0; ind<deleted_nodes.size(); ++ind)
    {
        mpMesh->DeleteNodePriorToReMesh(deleted_nodes[ind]->GetIndex());
    }

    for (unsigned ind=0; ind<deleted_faces.size(); ++ind)
    {
        mpMesh->DeleteFacePriorToReMesh(deleted_faces[ind]->GetIndex());
    }
    mpMesh->ReMesh();

    // Operation III
    assert(mpMesh->GetNumNodes()%2 == 0);   // LCOV_EXCL_LINE
    const unsigned half_num_nodes = mpMesh->GetNumNodes()/2;

    for (unsigned node_index=0; node_index<half_num_nodes; ++node_index)
    {
        Node<3>* basal_node = mpMesh->GetNode(node_index);
        c_vector<double, 3>& basal_position_ref = basal_node->rGetModifiableLocation();
        basal_position_ref[1] = length*basal_position_ref[1]/widthY;
        const double basal_theta = 2*M_PI*basal_position_ref[0]/widthX;
        basal_position_ref[0] = radius*cos(basal_theta);
        basal_position_ref[2] = radius*sin(basal_theta);

        Node<3>* apical_node = mpMesh->GetNode(node_index+half_num_nodes);
        c_vector<double, 3>& apical_position_ref = apical_node->rGetModifiableLocation();
        apical_position_ref[1] = length*apical_position_ref[1]/widthY;
        const double apical_theta = apical_position_ref[0]/widthX*2*M_PI;
        apical_position_ref[0] = (radius+thickness)*cos(apical_theta);
        apical_position_ref[2] = (radius+thickness)*sin(apical_theta);
    }
    return mpMesh;
}


void MonolayerVertexMeshGenerator::ClearStoredMeshObjects()
{
    mBasalNodes.clear();
    mApicalNodes.clear();
    mNodeLateralFaceMap.clear();
    mFaces.clear();
    mElements.clear();
    if (mpMesh)
        delete mpMesh;
}

void MonolayerVertexMeshGenerator::WriteVtk(const std::string& outputFile, const std::string& additionalTag,
                                            const bool usingFaceId)
{
    if (mpWriter != NULL)
    {
        delete mpWriter;
    }
    mpWriter = new VertexMeshWriter<3, 3>(outputFile, mName, false);
    mpWriter->WriteVtkUsingMeshWithCellId(*mpMesh, additionalTag, usingFaceId);
}

void MonolayerVertexMeshGenerator::WriteVtkWithSubfolder(const std::string& outputFile,
                                                         const std::string& additionalTag,
                                                         const bool usingFaceId)
{
    if (mpWriter != NULL)
    {
        delete mpWriter;
    }
    mpWriter = new VertexMeshWriter<3, 3>(outputFile+"/"+mName, mName, false);
    mpWriter->WriteVtkUsingMeshWithCellId(*mpMesh, additionalTag, usingFaceId);
}

void MonolayerVertexMeshGenerator::PrintMesh(const bool printDeletedObjects) const
{
    const std::string TAB = "    " ;

    std::cout <<"=================================================================================" << std::endl;
    // Printing out each elements
    const unsigned num_elems = printDeletedObjects ? mpMesh->GetNumAllElements() : mpMesh->GetNumElements();
    for (unsigned i=0; i<num_elems; ++i)
    {
        VertexElement<3,3>& elem = *(mpMesh->GetElement(i));
        std::cout << "ELEMENT (" << i<< ") : " << elem.GetIndex() << (elem.IsDeleted()?" (DELETED)": "") << std::endl;

        std::cout << TAB << "number of Faces : " << elem.GetNumFaces() << " {";
        for (unsigned j=0; j<elem.GetNumFaces(); ++j)
        {
            std::cout << std::setw(3) << elem.GetFace(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "Face oriented.. : " << elem.GetNumFaces() << " {";
        for (unsigned j=0; j<elem.GetNumFaces(); ++j)
        {
            std::cout << std::setw(3) << elem.FaceIsOrientatedAntiClockwise(j) << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "number of Nodes : " << elem.GetNumNodes() << " {  ";
        for (unsigned j=0; j<elem.GetNumNodes(); ++j)
        {
            std::cout << elem.GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl;

        VertexElement<2,3>& basal = *(elem.GetFace(0));
        std::cout << TAB << "Nodes for basal face " << basal.GetIndex() << " {  ";
        for (unsigned j=0; j<basal.GetNumNodes(); ++j)
        {
            std::cout << basal.GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl << "---------------------------------------------------------" << std::endl;
    }

    std::cout <<"***************************************************************" << std::endl;
    // Now printing all faces
    const unsigned num_faces = printDeletedObjects ? mpMesh->GetNumAllFaces() : mpMesh->GetNumFaces();
    for (unsigned i=0; i<num_faces; ++i)
    {
        VertexElement<2, 3>& face = *(mpMesh->GetFace(i));
        std::cout << "FACE (" << i<< ") : " << face.GetIndex() << (face.IsDeleted()?" (DELETED)": "") << std::endl;
        std::cout << TAB << "Face Attribute : " << face.rGetElementAttributes()[0] << (IsFaceOnBoundary(&face)?" (BOUNDARY)": "") << std::endl;

        std::set<unsigned> set_tmp = face.rGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (std::set<unsigned>::iterator it=set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "number of Nodes : " << face.GetNumNodes() << " {  ";
        for (unsigned j=0; j<face.GetNumNodes(); ++j)
        {
            std::cout << face.GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl << "---------------------------------------------------------" << std::endl;
    }

    std::cout <<"***************************************************************" << std::endl;
    //Now printing all the nodes
    const unsigned num_nodes = printDeletedObjects ? mpMesh->GetNumAllNodes() : mpMesh->GetNumNodes();
    for (unsigned i=0; i<num_nodes; ++i)
    {
        Node<3>& node = *(mpMesh->GetNode(i));
        std::set<unsigned> set_tmp = node.rGetContainingElementIndices();
        std::cout << "NODE (" << i<< ") : " << node.GetIndex() << (node.IsDeleted()?" (DELETED)": "") << std::endl;
        std::cout << TAB << "Node Attribute : " << node.rGetNodeAttributes()[0] << (node.IsBoundaryNode()?" (BOUNDARY)": "") << std::endl;
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (std::set<unsigned>::iterator it=set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl << "---------------------------------------------------------" << std::endl;
    }
}

