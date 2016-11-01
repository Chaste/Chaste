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

#ifndef TESTVERTEXMESH33UNIAXIALLOAD_HPP_
#define TESTVERTEXMESH33UNIAXIALLOAD_HPP_

#include "Debug.hpp"
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractForce.hpp"

#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "MisraForce.hpp"
#include "OffLatticeSimulation.hpp"

#include "FakePetscSetup.hpp"


#include "UblasCustomFunctions.hpp"
#define ADDFORCEPARAMETER p_force3->SetVolumeParameter(0, 10);
#define CURRENT_TEST std::string("Volume")
#define END_TIME 0.05

#include "HoneycombVertexMeshGenerator.hpp"
#define OUTPUT_NAME "TestUniaxialLoad"
class MeshBuilderHelper
{
private:
    std::string mName;
    std::string mAdditionalPath;
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
    MeshBuilderHelper(const std::vector<Node<3>*>& rLowerNodes, const std::string& AdditionalPath = "",
                      const std::string& Name = "mesh",const unsigned zHeight = 1)
                        : mName(Name),
                          mAdditionalPath("/" + AdditionalPath),
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

    MeshBuilderHelper(const std::string& AdditionalPath = "", const std::string& Name = "mesh")
    : mName(Name),
      mAdditionalPath("/" + AdditionalPath),
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
        assert ( mUpperNodes.size() == mLowerNodes.size() );

        mLowerNodes.insert(this->mLowerNodes.end(), mUpperNodes.begin(), mUpperNodes.end());

        mpMesh = new MutableVertexMesh<3, 3>(mLowerNodes, mElements);
        return mpMesh;
    }

    void PrintMesh(const bool allElements=false) const
    {
        const unsigned num_elems = allElements ? mpMesh->GetNumAllElements() : mpMesh->GetNumElements();
        for (unsigned i=0; i<num_elems; ++i)
        {
            VertexElement<3,3>& elem = *(mpMesh->GetElement(i));
            std::cout << "ELEMENT (" << i<< ") : " << elem.GetIndex() << std::endl;
            std::cout << "number of Faces : " << elem.GetNumFaces() << " {  ";
            for (unsigned j=0; j<elem.GetNumFaces(); ++j)
            {
                std::cout << elem.GetFace(j)->GetIndex() << "  ";
            }
            std::cout << "}" << std::endl;
            std::cout << "Face oriented.. : " << elem.GetNumFaces() << " {  ";
            for (unsigned j=0; j<elem.GetNumFaces(); ++j)
            {
                std::cout << elem.FaceIsOrientatedAntiClockwise(j) << "  ";
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
            mpWriter = new VertexMeshWriter<3, 3>(OUTPUT_NAME + mAdditionalPath, mName, false);
        }
        else
        {
            // current workaround
            delete mpWriter;
            mpWriter = new VertexMeshWriter<3, 3>(OUTPUT_NAME + mAdditionalPath, mName, false);
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
//template<unsigned DIM>
//class HorizontalStretchForce : public AbstractForce<DIM>
//{
//private:
//
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
//        archive & mForceMagnitude;
//        archive & mRelativeWidth;
//    }
//
//    /**
//     * Internal variable to save the magnitude
//     */
//    double mForceMagnitude;
//
//    /**
//     * The relative width from the boundary upon with the force is acting
//     */
//    double mRelativeWidth;
//
//public:
//    /**
//     * default constructor
//     */
//    HorizontalStretchForce(const double ForceMagnitude=1, const double RelativeWidth=0.1)
//    : mForceMagnitude(ForceMagnitude),
//      mRelativeWidth(RelativeWidth)
//    {
//MARK
//    }
//
//    virtual ~HorizontalStretchForce()
//    {
//    }
//
//    /**
//     * Adding force contributions to the Populations
//     * @param rCellPopulation
//     */
//    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
//    {
//        // Throw an exception message if not using a VertexBasedCellPopulation
//        ///\todo check whether this line influences profiling tests - if so, we should remove it.
//        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
//        {
//            EXCEPTION("MisraForce is to be used with a VertexBasedCellPopulation only");
//        }
//
//        // Define some helper variables
//        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
//        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();
//
//        const unsigned num_nodes = p_cell_population->GetNumNodes();
//        c_vector<double, DIM> force_in_positive_y = zero_vector<double>(DIM);
//        force_in_positive_y[0] = mForceMagnitude;
//        c_vector<double, DIM> force_in_negative_y = -force_in_positive_y;
//
//        ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
//        const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
//        const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
//        assert(max_x > min_x);
//        const double width_x = max_x - min_x;
//        const double width = width_x*mRelativeWidth;
//        const double left_bound = min_x + width;
//        const double right_bound = max_x - width;
//        // could have built two more Cuboids to check if the point is inside,
//        // but will have a lot of wrapping which reduce efficiency.
//        // Iterate over nodes in the cell population
//        for (unsigned node_global_index=0; node_global_index<num_nodes; ++node_global_index)
//        {
//            Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
// //            if ( ! p_this_node->IsBoundaryNode())
// //                continue;
//            const double loc_x = p_this_node->rGetLocation()[0];
//            if (loc_x < left_bound)
//            {
//                p_this_node->AddAppliedForceContribution(force_in_negative_y);
//            }
//            else if (loc_x > right_bound)
//            {
//                p_this_node->AddAppliedForceContribution(force_in_positive_y);
//            }
//        }
//    }
//
//    /**
//     * For the compiler now
//     * @param rParamsFile
//     */
//    virtual void OutputForceParameters(out_stream& rParamsFile)
//    {
//MARK
//        *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
//        *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";
//
//        // Call method on direct parent class
//        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
//    }
//};
//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)
//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)

template<unsigned DIM>
class From2dGeneral3dVertexModelForce : public AbstractForce<DIM>
{
private:

    double mTargetApicalArea;
    double mApicalAreaParameter;
    double mApicalEdgeParameter;
    double mTargetBasalArea;
    double mBasalAreaParameter;
    double mBasalEdgeParameter;
    double mLateralEdgeParameter;
    double mTargetVolume;
    double mVolumeParameter;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mTargetApicalArea;
    }

    c_vector<double, DIM> GetVolumeGradientofElementAtNode(const VertexElement<DIM, DIM>* pElement, const unsigned globalIndex) const
    {
        ///\todo check this derivation again
        // If I have done it correctly, it should look like
        //
        // sum_{j, such that face j containing node i} { 1/6*VectorProduct( C_j, r_{i+1} - r_{i-1}) + A_j/(3*n_j).
        // where C_j (vector) is the centroid of face j, A_j (vector) is the area of face j,
        // r_{i+1}/r_{i-1} is the next/previous node of node i in that face.

        assert(DIM==3 && DIM==3);
        c_vector<double, DIM> volume_gradient = zero_vector<double>(3);

        for (unsigned face_index=0; face_index<pElement->GetNumFaces(); ++face_index)
        {
            const VertexElement<DIM-1, DIM>* p_face = pElement->GetFace(face_index);
            const bool face_orientation = pElement->FaceIsOrientatedAntiClockwise(face_index);
            const unsigned this_local_index = p_face->GetNodeLocalIndex(globalIndex);
            const unsigned num_nodes = p_face->GetNumNodes();
            c_vector<double, DIM> this_face_gradient_contribution = zero_vector<double>(3);

            // if this face contains the node, it will be some number, not UINT_MAX. Second statement just as safety precaution
            if ( this_local_index!=UINT_MAX && this_local_index<num_nodes)
            {
                const unsigned next_local_index = (this_local_index+1)%num_nodes;
                const unsigned previous_local_index = (this_local_index-1+num_nodes)%num_nodes;
                const c_vector<double, DIM> previous_to_next = p_face->GetNodeLocation(next_local_index) - p_face->GetNodeLocation(previous_local_index);

                c_vector<double, DIM> face_centroid = p_face->GetCentroid();
                this_face_gradient_contribution += VectorProduct(face_centroid, previous_to_next) / 6.0;

                c_vector<double, 3> face_area = zero_vector<double>(3);
                for (unsigned current_running_index=0; current_running_index<num_nodes; ++current_running_index)
                {
                    // We add an extra num_nodes_in_element in the line below as otherwise this term can be negative, which breaks the % operator
                    const unsigned previous_running_index = (current_running_index+num_nodes-1)%num_nodes;

                    const c_vector<double, DIM> current_node_relative_location = p_face->GetNodeLocation(current_running_index) - face_centroid;
                    const c_vector<double, DIM> previous_node_relative_location = p_face->GetNodeLocation(previous_running_index) - face_centroid;

                    face_area += VectorProduct(current_node_relative_location, previous_node_relative_location) / 2.0;
                }

                this_face_gradient_contribution += face_area/(3.0*num_nodes);
            }

            // If it is oriented the other way, then r_{i+1} and r_{i-1} would be exchanged,
            // and the face_area (vector) will have different sign.
            if (face_orientation)
            {
                this_face_gradient_contribution *= -1.0;
            }

            volume_gradient += this_face_gradient_contribution;
        }

        return volume_gradient;
    }

public:
    From2dGeneral3dVertexModelForce()
        : AbstractForce<DIM>(),
          mTargetApicalArea(0),
          mApicalAreaParameter(0),
          mApicalEdgeParameter(0),
          mTargetBasalArea(0),
          mBasalAreaParameter(0),
          mBasalEdgeParameter(0),
          mLateralEdgeParameter(0),
          mTargetVolume(0),
          mVolumeParameter(0)
    {
    }

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
TRACE("Starting force calculation")
        c_vector<double, DIM> force = zero_vector<double>(DIM);
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("3dVertexModelForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();
        const unsigned num_nodes = p_cell_population->GetNumNodes();
        const unsigned num_elements = p_cell_population->GetNumElements();

        // Begin by computing the volumes of each element in the mesh, to avoid having to do this multiple times
        std::vector<double> element_volumes(num_elements);
        std::vector<double> apical_areas(num_elements);
        std::vector<double> basal_areas(num_elements);
        for (unsigned elem_index = 0 ; elem_index<num_elements ; ++elem_index)
        {
            const VertexElement<DIM, DIM>* p_elem = p_cell_population->GetElement(elem_index);
            assert(elem_index == p_elem->GetIndex());
            element_volumes[elem_index] = rMesh.GetVolumeOfElement(elem_index);
            apical_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(1));
            basal_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(0));
        }

        // Iterate over nodes in the cell population
        for (unsigned node_index=0 ; node_index<num_nodes ; node_index++)
        {
            Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
            assert(node_index == p_this_node->GetIndex());
            // Get the type of node. 1=basal; 2=apical
            const unsigned node_type = unsigned(p_this_node->rGetNodeAttributes()[0]);
PRINT_2_VARIABLES(node_index, node_type)
            c_vector<double, DIM> basal_face_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> ab_edge_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> apical_face_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> lateral_edge_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> volume_contribution = zero_vector<double>(DIM);

            const unsigned opposite_node_index = node_index + num_nodes/2*(node_type==1u?1:-1);
            const Node<DIM>* p_opposite_node = p_cell_population->GetNode(opposite_node_index);
            const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_opposite_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_opposite_node->rGetLocation());
            lateral_edge_contribution -= edge_gradient*mLateralEdgeParameter;

            // a variable to store such that the apical/basal edge forces are not counted twice for non-boundary edges.
            std::set<unsigned> neighbour_node_indices;

            // Find the indices of the elements owned by this node
            const std::set<unsigned> containing_elem_indices = p_this_node->rGetContainingElementIndices();

            // Iterate over these elements
            for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                    iter != containing_elem_indices.end();
                    ++iter)
            {
PRINT_VARIABLE(*iter)
                // Get this element, its index and its number of nodes
                VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                std::vector<VertexElement<DIM-1, DIM>*> lateral_faces;

                // Populate the pointers/vector to different types of face
                for (unsigned face_index=0; face_index<p_element->GetNumFaces(); ++face_index)
                {
                    VertexElement<DIM-1,DIM>* p_tmp_face = p_element->GetFace(face_index);
                    switch (unsigned(p_tmp_face->rGetElementAttributes()[0]))
                    {
                    case 1:
                        assert(face_index == 0);
                        break;
                    case 2:
                        assert(face_index == 1);
                        break;
                    case 3:
                        lateral_faces.push_back(p_tmp_face);
                        break;
                    default:
                        NEVER_REACHED;
                    }
                }

                const unsigned elem_index = p_element->GetIndex();

                // Calculating volume contribution
                c_vector<double, DIM> element_volume_gradient = this->GetVolumeGradientofElementAtNode(p_element, node_index);
                // Add the force contribution from this cell's volume compressibility (note the minus sign)
                volume_contribution -= mVolumeParameter*element_volume_gradient*(element_volumes[elem_index] - mTargetVolume);
if (mVolumeParameter>0.001)
PRINT_VECTOR((mVolumeParameter*element_volume_gradient*(element_volumes[elem_index] - mTargetVolume)))
                // Pointer to the face having the same type as the node
                const VertexElement<DIM-1, DIM>* p_ab_face = p_element->GetFace(node_type - 1);
                const unsigned local_node_index_in_ab_face = p_ab_face->GetNodeLocalIndex(node_index);
                const c_vector<double, DIM> ab_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_ab_face, local_node_index_in_ab_face);
                // Calculating apical face contribution
                if (node_type == 2)
                {
if (mApicalAreaParameter>0.001)
PRINT_VECTOR((mApicalAreaParameter*ab_face_gradient*(apical_areas[elem_index] - mTargetApicalArea)))
                    apical_face_contribution -= mApicalAreaParameter*ab_face_gradient*(apical_areas[elem_index] - mTargetApicalArea);
                }
                // Computing basal face contribution
                if (node_type == 1)
                {
if (mBasalAreaParameter>0.001)
PRINT_VECTOR((mBasalAreaParameter*ab_face_gradient*(basal_areas[elem_index] - mTargetBasalArea)))
                    basal_face_contribution -= mBasalAreaParameter*ab_face_gradient*(basal_areas[elem_index] - mTargetBasalArea);
                }
                const unsigned num_nodes_in_ab_face = p_ab_face->GetNumNodes();
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face+1)%num_nodes_in_ab_face));
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face-1+num_nodes_in_ab_face)%num_nodes_in_ab_face));
            }

            for(std::set<unsigned>::iterator it = neighbour_node_indices.begin(); it!=neighbour_node_indices.end() ; ++it)
            {

                Node<DIM>* p_neighbour_node = p_cell_population->GetNode(*it);
                const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_neighbour_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_neighbour_node->rGetLocation());
                ab_edge_contribution -= edge_gradient*(node_type==1u ? mBasalEdgeParameter : mApicalEdgeParameter);
            }

            c_vector<double, DIM> force_on_node = basal_face_contribution + ab_edge_contribution + apical_face_contribution
                                                  + lateral_edge_contribution + volume_contribution;
            p_this_node->AddAppliedForceContribution(force_on_node);
        }
    }

    void SetApicalParameter(const double TargetArea, const double AreaParameter, const double LineParameter)
    {
        mTargetApicalArea = TargetArea;
        mApicalAreaParameter = AreaParameter;
        mApicalEdgeParameter = LineParameter;
    }
    void SetBasalParameter(const double TargetArea, const double AreaParameter, const double LineParameter)
    {
        mTargetBasalArea = TargetArea;
        mBasalAreaParameter = AreaParameter;
        mBasalEdgeParameter = LineParameter;
    }
    void SetLateralParameter(const double Parameter)
    {
        mLateralEdgeParameter = Parameter;
    }
    void SetVolumeParameter(const double TargetVolume, const double VolumeParameter)
    {
        mTargetVolume = TargetVolume;
        mVolumeParameter = VolumeParameter;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
        *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
        *rParamsFile << "\t\t\t<TargetBasalArea>" << mTargetBasalArea << "</TargetBasalArea>\n";
        *rParamsFile << "\t\t\t<TargetVolume>" << mTargetVolume << "</TargetVolume>\n";
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};
//template class From2dGeneral3dVertexModelForce<3>;
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(From2dGeneral3dVertexModelForce,3)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(From2dGeneral3dVertexModelForce,3)


class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
public:

//    void TestVertexMesh33UiaxialLoad() throw (Exception)
//    {
//        // Make a mesh of 10x5
//        const double z_height = 1;
//        const double target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer("TestUniaxialLoad/Test1", "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//
//        std::vector<CellPtr> cells;
//        CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<3> simulator(cell_population);
//        simulator.SetOutputDirectory("TestUniaxialLoad/Test1");
//        simulator.SetSamplingTimestepMultiple(10);
//        const double end_time = 0.5;
//        simulator.SetEndTime(end_time);
//
// //        MAKE_PTR(MisraForce<3>, p_force);
// //        p_force->SetLateralSurfaceEnergyParameter(0);
// //        p_force->SetApicalLineTensionParameter(0);
// //        p_force->SetBasalSurfaceEnergyParameter(0);
// //        p_force->SetTargetVolume(z_height*target_area*2);
// //        p_force->SetVolumeCompressibilityParameter(2);
// //        simulator.AddForce(p_force);
// //        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
// //        simulator.AddForce(p_force2);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//        p_force3->SetApicalParameter(0, 0, 1);
//        simulator.AddForce(p_force3);
//
//        simulator.Solve();
//
//        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
//        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
//    }

    // Testing and therefore highly un
//    void TestApical1() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Apical1";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetApicalParameter(0, 0, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

//    void TestApical2() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Apical2";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetApicalParameter(1, 1, 0);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

//    void TestApical3() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Apical3";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetApicalParameter(1, 1, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//
//    void TestBasal1() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Basal1";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetBasalParameter(0, 0, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//    void TestBasal2() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Basal2";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetBasalParameter(1, 1, 0);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//    void TestBasal3() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Basal3";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetBasalParameter(1, 1, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

//    void TestApicalAndBasal1() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/ApicalAndBasal1";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetApicalParameter(1, 1, 1);
//        p_force3->SetBasalParameter(1, 1, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//    void TestApicalAndBasal2() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/ApicalAndBasal2";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetApicalParameter(1, 1, 1);
//        p_force3->SetBasalParameter(2, 2, 2);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

//    void TestLateral() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Lateral";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetLateralParameter(1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

//    void TestVolume1() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Volume1";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetVolumeParameter(0, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//    void TestVolume2() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Volume2";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetVolumeParameter(z_height*target_area, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }
//    void TestVolume3() throw (Exception)
//    {
//        const std::string output_file = "TestUniaxialLoad/Volume3";
//        const double z_height = 0.5, target_area = 1;
//        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
//        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
//        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
//        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
//        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
//
//        p_force3->SetVolumeParameter(z_height*target_area*3, 1);
//        simulator.AddForce(p_force3); simulator.Solve();
//    }

    void TestHexagon() throw (Exception)
    {
        const double z_height = 0.5, target_area = 1;
        HoneycombVertexMeshGenerator generator(1,1, false, 0.01, 0.001, target_area);
        MeshBuilderHelper builder("Hexagon" + CURRENT_TEST);
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(*(generator.GetMesh()), z_height);
        builder.PrintMesh();
        builder.WriteVtk("Hexagon");
        const std::string output_file = "TestUniaxialLoad/Hexagon" + CURRENT_TEST;


        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); simulator.SetEndTime(END_TIME);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        ADDFORCEPARAMETER
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestHexagons() throw (Exception)
        {
        const double z_height = 0.5, target_area = 1;
        HoneycombVertexMeshGenerator generator(2, 2, false, 0.01, 0.001, target_area);
        MeshBuilderHelper builder("Hexagons" + CURRENT_TEST);
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(*(generator.GetMesh()), z_height);
        builder.PrintMesh();
        builder.WriteVtk("Hexagons");
        const std::string output_file = "TestUniaxialLoad/Hexagons" + CURRENT_TEST;


        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); simulator.SetEndTime(END_TIME);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        ADDFORCEPARAMETER
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestSquare() throw (Exception)
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 0.0));
        unsigned node_indices_elem_0[4] = {0, 1, 2, 3};
        MeshBuilderHelper builder(nodes, "Square" + CURRENT_TEST);
        builder.buildElementWith(4, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Square");
        builder.PrintMesh();
        const std::string output_file = "TestUniaxialLoad/Square" + CURRENT_TEST;

        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, vertex_mesh.GetNumElements());
        VertexBasedCellPopulation<3> cell_population(vertex_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); simulator.SetEndTime(END_TIME);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        ADDFORCEPARAMETER
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestSquares() throw (Exception)
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true,  2.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(5, true,  2.0, 1.0, 0.0));
        unsigned node_indices_elem_0[4] = {0, 1, 2, 3};
        unsigned node_indices_elem_1[4] = {1, 4, 5, 2};
        MeshBuilderHelper builder(nodes, "Squares" + CURRENT_TEST);
        builder.buildElementWith(4, node_indices_elem_0);
        builder.buildElementWith(4, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtk("Squares");
        builder.PrintMesh();
        const std::string output_file = "TestUniaxialLoad/Squares" + CURRENT_TEST;

        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, vertex_mesh.GetNumElements());
        VertexBasedCellPopulation<3> cell_population(vertex_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); simulator.SetEndTime(END_TIME);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        ADDFORCEPARAMETER
        simulator.AddForce(p_force3); simulator.Solve();
    }
};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
