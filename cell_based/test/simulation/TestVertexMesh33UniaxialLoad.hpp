/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "AbstractCellBasedTestSuite.hpp"

#include "VoronoiVertexMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "GeodesicSphere23Generator.hpp"

#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "HorizontalStretchForce.hpp"
#include "OffLatticeSimulation.hpp"
#include <boost/lexical_cast.hpp>

// Cell writers
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellVolumesWriter.hpp"

#include "Debug.hpp"

#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"


#include "FakePetscSetup.hpp"

class LateralNodeModifier : public AbstractCellBasedSimulationModifier<3, 3>
{
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<3, 3> >(*this);
    }

public:
    /**
     * Default constructor.
     */
    LateralNodeModifier()
    : AbstractCellBasedSimulationModifier<3, 3>()
    {
    }

    /**
     * Destructor.
     */
    virtual ~LateralNodeModifier()
    {
    }

    /**
     * Change the lateral node by its definition at the end of time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<3, 3>& rCellPopulation)
    {
        UpdateCellData(rCellPopulation);
    }

    /**
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<3, 3> &rCellPopulation, std::string outputDirectory)
    {
        UpdateCellData(rCellPopulation);
    }

    void UpdateCellData(AbstractCellPopulation<3, 3>& rCellPopulation)
    {
        AbstractMesh<3, 3>& r_tmp_mesh = rCellPopulation.rGetMesh();

        if(dynamic_cast<MutableVertexMesh<3, 3>*>(&r_tmp_mesh) == NULL)
        {
            EXCEPTION("only to vertex mesh");
        }
        MutableVertexMesh<3, 3>* p_mesh = static_cast<MutableVertexMesh<3, 3>*>(&r_tmp_mesh);

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            Node<3>* p_lateral_node = p_mesh->GetNode(i);
            if (IsLateralNode(p_lateral_node))
            {
                const std::set<unsigned>& containing_faces = p_lateral_node->rGetContainingFaceIndices();
                std::vector<Node<3>*> basal_nodes;
                std::vector<Node<3>*> apical_nodes;

                for (std::set<unsigned>::const_iterator it = containing_faces.begin();
                it != containing_faces.end(); ++it)
                {
                    VertexElement<2, 3>* p_face_tmp = p_mesh->GetFace(*it);
                    if (p_face_tmp->GetNumNodes() != 3)
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

                
                const c_vector<double, 3> mid_apical = (apical_nodes[0]->rGetLocation() + apical_nodes[1]->rGetLocation())/2;
                const c_vector<double, 3> mid_basal = (basal_nodes[0]->rGetLocation() + basal_nodes[1]->rGetLocation())/2;

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
                    if (p_small_triangular_face == NULL || p_small_triangular_face->GetNumNodes() !=3u || p_small_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
                    {
                        NEVER_REACHED;
                    }

                    VertexElement<2, 3>* p_big_triangular_face = GetSharedLateralFace(p_node_x, p_node_y, p_mesh);
                    if (p_big_triangular_face == NULL || p_big_triangular_face->GetNumNodes() !=3u || p_big_triangular_face->GetNodeLocalIndex(p_lateral_node->GetIndex()) == UINT_MAX)
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

    /**
     * Output any simulation modifier parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        AbstractCellBasedSimulationModifier<3>::OutputSimulationModifierParameters(rParamsFile);
    }
};

class BielmeierExternalForce : public AbstractForce<3>
{
    double mEcmSpringConstant;

    double mT_ext;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<3> >(*this);
    }

  public:
    BielmeierExternalForce(const double springConstant = 0.0, const double T_ext = 0.0)
        : AbstractForce<3>(),
          mEcmSpringConstant(springConstant),
          mT_ext(T_ext)
    {
    }

    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
    {
        if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("GeneralMonolayerVertexMeshForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);
        MutableVertexMesh<3, 3>& rMesh = p_cell_population->rGetMesh();

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            Node<3>* p_node = rMesh.GetNode(i);

            if (IsBasalNode(p_node))
            {
                c_vector<double, 3> result = -1 * mEcmSpringConstant * Create_c_vector(0, 0, p_node->rGetLocation()[2]);
                p_node->AddAppliedForceContribution(result);
            }
        }

        for (unsigned i = 0; i < rMesh.GetNumFaces(); ++i)
        {
            VertexElement<2, 3>* p_face = rMesh.GetFace(i);

            if (IsApicalFace(p_face))
            {
                for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
                {
                    c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, j);
                    result *= mT_ext;
                    p_face->GetNode(j)->AddAppliedForceContribution(result);
                }
            }

        }
    }

    void OutputForceParameters(out_stream &rParamsFile)
    {
        *rParamsFile << "\t\t\t<EcmSpringConstant>" << mEcmSpringConstant << "</EcmSpringConstant>\n";
        *rParamsFile << "\t\t\t<T_ext>" << mT_ext << "</T_ext>\n";
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BielmeierExternalForce)
// CHASTE_CLASS_EXPORT(LateralNodeModifier)

#include "SerializationExportWrapperForCpp.hpp"
// CHASTE_CLASS_EXPORT(LateralNodeModifier)
CHASTE_CLASS_EXPORT(BielmeierExternalForce)

class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
private:
    static const double z_height = 1;
    static const double target_area = 1;
    const unsigned num_cells_x = 9;
    const unsigned num_cells_y = 5;
    static const double end_time = 10;

public:
    void TestOnHexagonalMesh() throw(Exception)
    {
        std::string output_filename = "TestUniaxialLoad/HoneyTest_ori" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
        builder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(400, target_area * 1.2);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.5);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);

        // loweredALLApical
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_loweredALLApical";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            // cell_population.AddCellWriter<CellIdWriter>();
            // cell_population.AddCellWriter<CellVolumesWriter>();

            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(5, 5, 0.2);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // LoweredApicalLineTension
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_LoweredApicalLineTension";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            // cell_population.AddCellWriter<CellIdWriter>();
            // cell_population.AddCellWriter<CellVolumesWriter>();

            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(5, 20, 0.7);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // LoweredApicalFaceParam
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_LoweredApicalFaceParam";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            // cell_population.AddCellWriter<CellIdWriter>();
            // cell_population.AddCellWriter<CellVolumesWriter>();

            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(20, 5, 0.2);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // AdjustedBasal (not yet implemented)
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_AdjustedBasal";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            // cell_population.AddCellWriter<CellIdWriter>();
            // cell_population.AddCellWriter<CellVolumesWriter>();

            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(20, 20, 0.7);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }
    }

    void TestOnVoronoiMesh() throw(Exception)
    {
        std::string output_filename = "TestUniaxialLoad/VoronoiTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
        builder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1.0);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestCylindricalMesh()
    {
        const unsigned x = 10;
        const unsigned y = 11;
        const double c_end_time = 5;
        const double a = 2;
        const double length = 3 * sqrt(3) * y + sqrt(3);
        const double radius = a / M_PI / 2 * x;
        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2 * sqrt(3));
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        std::string output_filename = "TestUniaxialLoad/CylinderTest_SmallVolume" + boost::lexical_cast<std::string>(x)
            + "x" + boost::lexical_cast<std::string>(y);
        builder.WriteVtk(output_filename, "InitialMesh");

        builder.ConvertMeshToCylinder(2 * x, 1, radius*0.8, 1.5, 1);

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            c_vector<double, 3>& tmp_loc = p_mesh->GetNode(i)->rGetModifiableLocation();
            double xx = tmp_loc[0];
            tmp_loc[0] = tmp_loc[1];
            tmp_loc[1] = -xx;
        }
        for (unsigned i = 0; i < p_mesh->GetNumElements(); ++i)
        {
            p_mesh->GetElement(i)->MonolayerElementRearrangeFacesNodes();
        }
        builder.WriteVtk(output_filename, "After");

        PRINT_4_VARIABLES(p_mesh->GetVolumeOfElement(0), p_mesh->GetVolumeOfElement(5), p_mesh->CalculateAreaOfFace(p_mesh->GetFace(0)), p_mesh->CalculateAreaOfFace(p_mesh->GetFace(1)))

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        // MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        // CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(c_end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(5, 5, 0.7);
        p_force3->SetBasalParameters(5, 5, 0.7);
        p_force3->SetLateralParameter(7);
        p_force3->SetVolumeParameters(100, 6);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x * y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), c_end_time, 1e-10);
    }

    void TestCellGrowth() throw(Exception)
    {
        // Make a mesh of 10x10
        //        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 6;
        const unsigned num_cells_y = 3;
        std::string output_filename = "TestCellDivision/HoneyTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder("CellGrowth");
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        builder.WriteVtk(output_filename, "Before");

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.2);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 55u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestOnSphere() throw(Exception)
    {
        const double s_end_time = 0.2;
        std::string output_filename = "TestUniaxialLoad/SphereTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        GeodesicSphere23Generator builder;
        builder.SubDivide(); // n=42
        builder.SubDivide(); // n=162

        MutableVertexMesh<2, 3>* p_dual_mesh = builder.GetDual();
        VertexMeshWriter<2, 3> Writer(output_filename, "Geodesic_Dual", false);
        Writer.WriteVtkUsingMesh(*p_dual_mesh);

        const unsigned radius = sqrt(p_dual_mesh->GetNumElements() * target_area / 4 / M_PI);
        MonolayerVertexMeshGenerator sBuilder;
        MutableVertexMesh<3, 3>* p_mesh = sBuilder.MakeSphericalMesh33(p_dual_mesh, 5, 0.5);
        sBuilder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        // CellsGenerator<NoCellCycleModel, 3> cells_generator;
        // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(s_end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 2);
        simulator.AddForce(p_force3);
        // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        // p_force2->SetForceMagnitude(1.0);
        // p_force2->SetRelativeWidth(0.15);
        // simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 176u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), s_end_time, 1e-10);
    }

    void TestSingleCube()
    {
        /*
         * Create a mesh with a square element, as shown below.
         *  _____
         * |     |
         * |     |
         * |_____|
         */
        std::string output_filename = "TestUniaxialLoad/SingleCubeTest";
        const double end_time = 4;
        const double target_volume = 3;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 1.0, 1.0));
        const unsigned node_indices_elem_0[4] = { 0, 1, 3, 2 };

        MonolayerVertexMeshGenerator builder(nodes, "TestGradientDeviation");
        builder.BuildElementWith(4, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>* p_mesh = builder.GenerateMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        // p_force3->SetApicalParameters(20, 20, 0.7);
        // p_force3->SetBasalParameters(20, 20, 0.7);
        // p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, target_volume);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), target_volume, 0.05);

        simulator.SetEndTime(end_time * 2);
        p_force3->SetVolumeParameters(350, target_volume / 2);
        simulator.Solve();
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), target_volume / 2, 0.05);
    }

    void TestAsynchronousT1()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap doesn't occur when max(l1,l2)>mCellRearrangementThreshold.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        // Set the threshold distance between vertices for a T1 swap as follows
        // so that it will not trigger CheckForSwapsFromShortEdges

        std::string output_filename = "TestAsynchronousT1/FirstTest";
        const double end_time = 1;
        const double target_volume = 2.2;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[4] = { 4, 1, 2, 5 };
        unsigned node_indices_elem_2[3] = { 0, 1, 4 };
        unsigned node_indices_elem_3[4] = { 4, 5, 3, 0 };

        MonolayerVertexMeshGenerator builder(nodes, "AsynchronousT1Swap");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        MutableVertexMesh<3, 3>* p_mesh = builder.GenerateMesh();

        p_mesh->GetNode(4)->rGetModifiableLocation()[1] = 0.46;
        p_mesh->GetNode(5)->rGetModifiableLocation()[1] = 0.54;
        // p_mesh->SetCellRearrangementThreshold(0.10);
        // p_mesh->SetCellRearrangementRatio(8);

        // p_mesh->ReMesh();
        // p_mesh->SetCellRearrangementThreshold(0.01);
        // p_mesh->SetCellRearrangementRatio(1.5);


        for (unsigned i = 0; i < p_mesh->GetNumNodes() ; ++i)
        {
            p_mesh->GetNode(i)->rGetModifiableLocation() *= 2;
        }

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(15, 15, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(9.25);
        p_force3->SetVolumeParameters(350, target_volume*1.2);
        simulator.AddForce(p_force3);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.50);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.SetEndTime(end_time);
        simulator.Solve();

        // const double lateral_params[6] = {8, 8.5, 9, 9.5, 10, 10.5};
        // for (unsigned i = 0; i < 6; ++i)
        // {
        //     p_force3->SetLateralParameter(lateral_params[i]);
        //     simulator.SetEndTime(end_time*(i+1));
        //     simulator.Solve();
        // }

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4);
        // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
