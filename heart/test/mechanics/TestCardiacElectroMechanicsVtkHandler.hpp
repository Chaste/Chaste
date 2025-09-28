/*

Copyright (c) 2005-2025, University of Oxford.
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


#ifndef TESTCARDIACELECTROMECHANICSVTKHANDLER_HPP_
#define TESTCARDIACELECTROMECHANICSVTKHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "FineCoarseMeshPair.hpp"
#include "ImplicitCardiacMechanicsSolver.hpp"
#include "SimpleStimulus.hpp"
#include "CardiacElectroMechanicsVtkHandler.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "PetscSetupAndFinalize.hpp"

// cell factory which stimulates everything at once - 2D version
class EntirelyStimulatedTissueCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    EntirelyStimulatedTissueCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-100000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
    }
};


// cell factory which stimulates everything at once - 3D version
class EntirelyStimulatedTissueCellFactory3D : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    EntirelyStimulatedTissueCellFactory3D()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-100000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
    }
};

/**
* The tests below build the VTK writer object
* then tests whether it actually writeds valid VTK files 
* with the expected data 
*/
class TestCardiacElectroMechanicsVtkHandler : public CxxTest::TestSuite
{
public:

    void TestEMVtkOutput2D()
    {
        EntirelyStimulatedTissueCellFactory cell_factory;

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.05, 0.05);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.025, 0.05, 0.05);

        MonodomainProblem<2> mono_problem( &cell_factory );
        mono_problem.SetMesh(&electrics_mesh);
        mono_problem.Initialise();
        std::string output_dir = "TestEMvtkWriter2D";

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > fixed_node_locations;

        // fix the node at the origin so that the solution is well-defined (ie unique)
        fixed_nodes.push_back(0);
        fixed_node_locations.push_back(zero_vector<double>(2));

        // for the rest of the nodes, if they lie on X=0, fix x=0 but leave y free.
        for (unsigned i=1 /*not 0*/; i<mechanics_mesh.GetNumNodes(); i++)
        {
            if (fabs(mechanics_mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                c_vector<double,2> new_position;
                new_position(0) = 0.0;
                new_position(1) = SolidMechanicsProblemDefinition<2>::FREE;
                fixed_nodes.push_back(i);
                fixed_node_locations.push_back(new_position);
            }
        }

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
        problem_defn.SetMechanicsSolveTimestep(1.0);


        ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<2>,2> mechanics_solver(
                                        mechanics_mesh,problem_defn,output_dir);

        FineCoarseMeshPair<2>* p_mesh_pair = new FineCoarseMeshPair<2>(electrics_mesh, mechanics_mesh);
        p_mesh_pair->SetUpBoxesOnFineMesh();
        p_mesh_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(mechanics_solver.GetQuadratureRule()), false);
        p_mesh_pair->DeleteFineBoxCollection();
        
        mechanics_solver.SetFineCoarseMeshPair(p_mesh_pair);
        mechanics_solver.Initialise();
        
        CardiacElectroMechanicsVtkHandler<2> handler(mechanics_solver,
                                      mechanics_mesh,
                                      electrics_mesh,
                                      mono_problem,
                                      output_dir);
        
        // fake a solution
        Vec fake_solution = mono_problem.CreateInitialCondition();
        ReplicatableVector fake_repl(fake_solution);
        fake_repl[0] = 987.5;//put a different number just for variety
        handler.WriteSolution(1u, fake_repl);

        {
            // Check that the reader can see it
            VtkMeshReader<2,2> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir+"/vtk/deformed_mechanics_mesh_1.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mechanics_mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mechanics_mesh.GetNumElements());

            //Check node locations (no deformation, no solve called)
            for (unsigned i = 0u; i < mechanics_mesh.GetNumNodes(); ++i)
            {
                std::vector<double> next_node = vtk_reader.GetNextNode();
                TS_ASSERT_EQUALS(next_node.size(),3);//the GetNextNode always fills up 3 coordinates regardless of SPACE_DIM (not sure why)
                TS_ASSERT_DELTA(next_node[0],mechanics_mesh.GetNode(i)->rGetLocation()[0],1e-9);
                TS_ASSERT_DELTA(next_node[1],mechanics_mesh.GetNode(i)->rGetLocation()[1],1e-9);
            }

            // Check that it has the correct data
            std::vector<double> voltage;
            vtk_reader.GetPointData("V", voltage);
            TS_ASSERT_EQUALS(voltage.size(),mechanics_mesh.GetNumNodes());
            for (unsigned i=0u; i < voltage.size(); ++i)
            {
                TS_ASSERT_DELTA(fake_repl[i] , voltage[i], 1e-9);
            }
            
            // No displacement expected
            std::vector< c_vector<double, 2> > displ;
            vtk_reader.GetPointData("displacements", displ);
            TS_ASSERT_EQUALS(displ.size(),mechanics_mesh.GetNumNodes());
            for (unsigned i=0u; i < displ.size(); ++i)
            {
                TS_ASSERT_DELTA(displ[i](0) , 0.0, 1e-9);
                TS_ASSERT_DELTA(displ[i](1) , 0.0, 1e-9);
            }
        }
        
        delete p_mesh_pair;

    }

    void TestEMVtkOutput3D()
    {
        MechanicsEventHandler::Reset();//prevent some warnings
        HeartEventHandler::Reset();//prevent some warnings
        
        EntirelyStimulatedTissueCellFactory3D cell_factory;

        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.05, 0.05, 0.05);

        QuadraticMesh<3> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.025, 0.05, 0.05, 0.05);

        MonodomainProblem<3> mono_problem( &cell_factory );
        mono_problem.SetMesh(&electrics_mesh);
        mono_problem.Initialise();
        std::string output_dir = "TestEMvtkWriter3D";

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,3> > fixed_node_locations;

        // fix the node at the origin so that the solution is well-defined (ie unique)
        fixed_nodes.push_back(0);
        fixed_node_locations.push_back(zero_vector<double>(3));

        // for the rest of the nodes, if they lie on X=0, fix x=0 but leave y free.
        for (unsigned i=1 /*not 0*/; i<mechanics_mesh.GetNumNodes(); i++)
        {
            if (fabs(mechanics_mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                c_vector<double,3> new_position;
                new_position(0) = 0.0;
                new_position(1) = SolidMechanicsProblemDefinition<2>::FREE;
                fixed_nodes.push_back(i);
                fixed_node_locations.push_back(new_position);
            }
        }

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
        problem_defn.SetMechanicsSolveTimestep(1.0);


        ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<3>,3> mechanics_solver(
                                        mechanics_mesh,problem_defn,output_dir);

        FineCoarseMeshPair<3>* p_mesh_pair = new FineCoarseMeshPair<3>(electrics_mesh, mechanics_mesh);
        p_mesh_pair->SetUpBoxesOnFineMesh();
        p_mesh_pair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(mechanics_solver.GetQuadratureRule()), false);
        p_mesh_pair->DeleteFineBoxCollection();
        
        mechanics_solver.SetFineCoarseMeshPair(p_mesh_pair);
        mechanics_solver.Initialise();
        
        CardiacElectroMechanicsVtkHandler<3> handler(mechanics_solver,
                                      mechanics_mesh,
                                      electrics_mesh,
                                      mono_problem,
                                      output_dir);
        
        // fake a solution
        Vec fake_solution = mono_problem.CreateInitialCondition();
        ReplicatableVector fake_repl(fake_solution);
        fake_repl[0] = 987.5;//put a different number just for variety
        handler.WriteSolution(1u, fake_repl);

        {
            // Check that the reader can see it
            VtkMeshReader<3,3> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir+"/vtk/deformed_mechanics_mesh_1.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mechanics_mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mechanics_mesh.GetNumElements());

            //Check node locations (no deformation, no solve called)
            for (unsigned i = 0u; i < mechanics_mesh.GetNumNodes(); ++i)
            {
                std::vector<double> next_node = vtk_reader.GetNextNode();
                TS_ASSERT_EQUALS(next_node.size(),3);//the GetNextNode always fills up 3 coordinates regardless of SPACE_DIM (not sure why)
                TS_ASSERT_DELTA(next_node[0],mechanics_mesh.GetNode(i)->rGetLocation()[0],1e-9);
                TS_ASSERT_DELTA(next_node[1],mechanics_mesh.GetNode(i)->rGetLocation()[1],1e-9);
                TS_ASSERT_DELTA(next_node[2],mechanics_mesh.GetNode(i)->rGetLocation()[2],1e-9);
            }

            // Check that it has the correct data
            std::vector<double> voltage;
            vtk_reader.GetPointData("V", voltage);
            TS_ASSERT_EQUALS(voltage.size(),mechanics_mesh.GetNumNodes());
            for (unsigned i=0u; i < voltage.size(); ++i)
            {
                TS_ASSERT_DELTA(fake_repl[i] , voltage[i], 1e-9);
            }
            
            // No displacement expected
            std::vector< c_vector<double, 3> > displ;
            vtk_reader.GetPointData("displacements", displ);
            TS_ASSERT_EQUALS(displ.size(),mechanics_mesh.GetNumNodes());
            for (unsigned i=0u; i < displ.size(); ++i)
            {
                TS_ASSERT_DELTA(displ[i](0) , 0.0, 1e-9);
                TS_ASSERT_DELTA(displ[i](1) , 0.0, 1e-9);
                TS_ASSERT_DELTA(displ[i](2) , 0.0, 1e-9);
            }
        }
        delete p_mesh_pair;
    }
  

};
#endif /*TESTCARDIACELECTROMECHANICSVTKHANDLER_HPP_*/
