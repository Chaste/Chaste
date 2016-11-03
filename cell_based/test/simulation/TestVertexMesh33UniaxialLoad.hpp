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
#include "VoronoiVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "SmartPointers.hpp"
//#include "MisraForce.hpp"
#include "OffLatticeSimulation.hpp"

#include "FakePetscSetup.hpp"

#include "GeneralMonolayerVertexMeshForce.hpp"

#define OUTPUT_NAME "TestUniaxialLoad/InitialMesh"
#define ADDFORCEPARAMETER p_force3->SetVolumeParameter(0.8, 10);
#define CURRENT_TEST std::string("Volume2")
#define END_TIME 0.5

#include "HoneycombVertexMeshGenerator.hpp"
template<unsigned DIM>
class HorizontalStretchForce : public AbstractForce<DIM>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mForceMagnitude;
        archive & mRelativeWidth;
    }

    /**
     * Internal variable to save the magnitude
     */
    double mForceMagnitude;

    /**
     * The relative width from the boundary upon with the force is acting
     */
    double mRelativeWidth;

public:
    /**
     * default constructor
     */
    HorizontalStretchForce(const double ForceMagnitude=1, const double RelativeWidth=0.1)
    : mForceMagnitude(ForceMagnitude),
      mRelativeWidth(RelativeWidth)
    {
    }

    virtual ~HorizontalStretchForce()
    {
    }

    /**
     * Adding force contributions to the Populations
     * @param rCellPopulation
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        // Throw an exception message if not using a VertexBasedCellPopulation
        ///\todo check whether this line influences profiling tests - if so, we should remove it.
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("MisraForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();

        const unsigned num_nodes = p_cell_population->GetNumNodes();
        c_vector<double, DIM> force_in_positive_y = zero_vector<double>(DIM);
        force_in_positive_y[0] = mForceMagnitude;
        c_vector<double, DIM> force_in_negative_y = -force_in_positive_y;

        ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
        const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
        const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
        assert(max_x > min_x);
        const double width_x = max_x - min_x;
        const double width = width_x*mRelativeWidth;
        const double left_bound = min_x + width;
        const double right_bound = max_x - width;
        // could have built two more Cuboids to check if the point is inside,
        // but will have a lot of wrapping which reduce efficiency.
        // Iterate over nodes in the cell population
        for (unsigned node_global_index=0; node_global_index<num_nodes; ++node_global_index)
        {
            Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
 //            if ( ! p_this_node->IsBoundaryNode())
 //                continue;
            const double loc_x = p_this_node->rGetLocation()[0];
            if (loc_x < left_bound)
            {
                p_this_node->ClearAppliedForce();
                p_this_node->AddAppliedForceContribution(force_in_negative_y);
            }
            else if (loc_x > right_bound)
            {
                p_this_node->ClearAppliedForce();
                p_this_node->AddAppliedForceContribution(force_in_positive_y);
            }
        }
    }

    /**
     * For the compiler now
     * @param rParamsFile
     */
    virtual void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
        *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";

        // Call method on direct parent class
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)

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
        const std::string TAB = "    " ;
        std::cout <<"=================================================================================" << std::endl;
        const unsigned num_elems = allElements ? mpMesh->GetNumAllElements() : mpMesh->GetNumElements();
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

        const unsigned num_faces = allElements ? mpMesh->GetNumAllFaces() : mpMesh->GetNumFaces();
        for (unsigned i=0; i<num_faces; ++i)
        {
            VertexElement<2, 3>& face = *(mpMesh->GetFace(i));
            std::cout << "FACE (" << i<< ") : " << face.GetIndex() << (face.IsDeleted()?" (DELETED)": "") << std::endl;
            std::cout << TAB << "Face Attribute : " << face.rGetElementAttributes()[0] << (face.IsElementOnBoundary()?" (BOUNDARY)": "") << std::endl;
            std::cout << TAB << "number of Nodes : " << face.GetNumNodes() << " {  ";
            for (unsigned j=0; j<face.GetNumNodes(); ++j)
            {
                std::cout << face.GetNode(j)->GetIndex() << "  ";
            }
            std::cout << "}" << std::endl << "---------------------------------------------------------" << std::endl;
        }
        std::cout <<"***************************************************************" << std::endl;

        const unsigned num_nodes = allElements ? mpMesh->GetNumAllNodes() : mpMesh->GetNumNodes();
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
        mpWriter->WriteVtkUsingMeshWithCellId(*mpMesh, AdditionalTag, true);
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
                faces_orientation.push_back(true);
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
                faces_orientation.push_back(false);
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


class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
public:

    void TestVertexMesh33UiaxialLoad() throw (Exception)
    {
        // Make a mesh of 10x5
        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 11;
        const unsigned num_cells_y = 5;
        // There seems to be a bug somewhere in voronoiprism3dVertexMeshGenerator....
//        VoronoiPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
//        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MeshBuilderHelper builder("UniaxialLoad");
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        builder.WriteVtk("Before");

        char tmp_name[50];
        sprintf(tmp_name, "TestUniaxialLoad/HoneyTest%dx%d", num_cells_x, num_cells_y);
        VertexMeshWriter<3, 3> vertex_mesh_writer(tmp_name, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(tmp_name);
        simulator.SetSamplingTimestepMultiple(10);
        const double end_time = 4;
        simulator.SetEndTime(end_time);

//        MAKE_PTR(MisraForce<3>, p_force);
//        p_force->SetLateralSurfaceEnergyParameter(0);
//        p_force->SetApicalLineTensionParameter(1);
//        p_force->SetBasalSurfaceEnergyParameter(1);
//        p_force->SetTargetVolume(z_height*target_area*2);
//        p_force->SetVolumeCompressibilityParameter(30);
//        simulator.AddForce(p_force);
        MAKE_PTR(GeneralMonolayerVertexMeshForce<3>, p_force3);
        p_force3->SetApicalParameter(5, 5, 1);
        p_force3->SetBasalParameter(5, 5, 1);
        p_force3->SetLateralParameter(2);
        p_force3->SetVolumeParameter(100, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR_ARGS(HorizontalStretchForce<3>, p_force2, (2));
        simulator.AddForce(p_force2);

            simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
