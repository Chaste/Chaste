/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTMONODOMAINPURKINJEPROBLEM_HPP_
#define TESTMONODOMAINPURKINJEPROBLEM_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "MixedDimensionMesh.hpp"
#include "PetscTools.hpp"
#include "LuoRudy1991.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscMatTools.hpp"
#include "TetrahedralMesh.hpp"
#include "MonodomainProblem.hpp"
#include "MonodomainPurkinjeProblem.hpp"
#include "NumericFileComparison.hpp"
#include "PurkinjeVentricularJunctionStimulus.hpp"
#include "MultiStimulus.hpp"
#include "ChasteSyscalls.hpp"

class NonPurkinjeCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    NonPurkinjeCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};


class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2> //inherits from different base
{
    friend class TestMonodomainPurkinjeProblem;
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    bool mStimulatePurkinje;
    bool mMakeJunction;

public:
    PurkinjeCellFactory(bool stimulatePurkinje = true, bool makeJunction = false)
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5)),
          mStimulatePurkinje(stimulatePurkinje),
          mMakeJunction(makeJunction)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        ChastePoint<2> location = pNode->GetPoint();

        AbstractCardiacCell* p_purkinje_cell;
        if (mStimulatePurkinje && fabs(location[0])<1e-6)
        {
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
        if (mMakeJunction && fabs(location[0] - 0.05)<1e-6)//Junction is in the middle of the square mesh at (0.05, 0.05)
        {
            double resistance = 1000.0; //kOhm
            CreateJunction(pNode, p_purkinje_cell, pCardiacCell, resistance);
        }
        return p_purkinje_cell;
    }
};

class PurkinjeStarCellFactory : public AbstractPurkinjeCellFactory<2>
{
    friend class TestMonodomainPurkinjeProblem;
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    bool mStimulatePurkinje;
    bool mMakeJunction;

public:
    PurkinjeStarCellFactory(bool stimulatePurkinje = true, bool makeJunction = true)
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5)),
          mStimulatePurkinje(stimulatePurkinje),
          mMakeJunction(makeJunction)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        ChastePoint<2> location = pNode->GetPoint();

        AbstractCardiacCell* p_purkinje_cell;
        if (fabs(location[0])<1e-6 && mStimulatePurkinje)
        {
            //Left side of the mesh is stimulated
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
        if(mMakeJunction)
        {
            /// Junctions - needs to have a dumb partition
            unsigned index = pNode->GetIndex();
            ///\todo What about 59 which is on the cable already?
            if (index == 60 ||index == 49 ||index == 50 ||index == 61 ||index == 71 ||index == 70 )//Junctions around the middle of the square mesh at (0.05, 0.05), node 60
            {
                double resistance = 1000.0; //kOhm
                if (index != 60)
                {
                    //The ones on the outside of the star get a different resistance
                    resistance = 1e4; //kOhm
                }
                CreateJunction(pNode, p_purkinje_cell, pCardiacCell, resistance);
            }
        }
        return p_purkinje_cell;
    }
};

class PurkinjeStarCellFactoryFromFile : public AbstractPurkinjeCellFactory<2>
{
    friend class TestMonodomainPurkinjeProblem;
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    bool mStimulatePurkinje;
    bool mMakeJunction;

public:
    PurkinjeStarCellFactoryFromFile(bool stimulatePurkinje = true, bool makeJunction = true)
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5)),
          mStimulatePurkinje(stimulatePurkinje),
          mMakeJunction(makeJunction)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        ChastePoint<2> location = pNode->GetPoint();

        AbstractCardiacCell* p_purkinje_cell;
        if (fabs(location[0])<1e-6 && mStimulatePurkinje)
        {
            //Left side of the mesh is stimulated
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            p_purkinje_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
        if(mMakeJunction)
        {
            CreateJunctionFromFile(pNode, p_purkinje_cell, pCardiacCell);
        }

        return p_purkinje_cell;
    }
};

class TestMonodomainPurkinjeProblem : public CxxTest::TestSuite
{
public:
    // Solve a Purkinje problem and check that the solution for the myocardium nodes is exactly the same as
    // the solution of an equivalent monodomain-only problem, that the purkinje voltage for
    // purkinje nodes doesn't change (no purkinje stimulus is given), and is 0 for non-purkinje nodes.
    // This version uses the problem classes.  It also checks things still work with a permuted mesh.
    void TestMonodomainPurkinjeProblemRunning()
    {
        // Settings common to both problems
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(1.0);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();

        ReplicatableVector soln_repl;
        ReplicatableVector soln_mono_repl;

        {
            // Set up normal monodomain
            NonPurkinjeCellFactory cell_factory_for_just_monodomain;
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_200_elements");
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_normal");
            MonodomainProblem<2,2> monodomain_problem(&cell_factory_for_just_monodomain);
            monodomain_problem.Initialise();

            // Solve
            monodomain_problem.Solve();
            soln_mono_repl.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // Set up Purkinje problem
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_purkinje");
        PurkinjeCellFactory cell_factory;
        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        MixedDimensionMesh<2, 2>& r_mesh = static_cast<MixedDimensionMesh<2, 2>& >(purkinje_problem.rGetMesh());

        //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
        //test
        for (MixedDimensionMesh<2, 2>::CableElementIterator iter = r_mesh.GetCableElementIteratorBegin();
             iter != r_mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            (*iter)->SetAttribute(1.0/std::sqrt(M_PI)); //So that the fibre cross section area is 1.0
        }

        purkinje_problem.SetWriteInfo(); // For coverage

        // Solve
        purkinje_problem.Solve();
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        // Compare solutions at final time
        for (AbstractTetrahedralMesh<2,2>::NodeIterator node_iter = r_mesh.GetNodeIteratorBegin();
             node_iter != r_mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            // Purkinje cable runs along y=0.05
            unsigned i = node_iter->GetIndex();
            if (fabs(node_iter->rGetLocation()[1] - 0.05) < 1e-6)
            {
                // The Purkinje domain is a 1D line embedded within the tissue.
                // It is stimulated in the same way as the tissue domain, therefore
                // the propagation velocity should be the same (within numerical error).
                TS_ASSERT_DELTA(soln_repl[2*i+1], soln_repl[2*i], 1e-1);
            }
            else
            {
                TS_ASSERT_DELTA(soln_repl[2*i+1], 0.0, 1e-4);
            }

            TS_ASSERT_DELTA(soln_repl[2*i], soln_mono_repl[i], 1e-5);
        }

        // Thorough comparison of output results - check that the solution in the 'main' tissue (not Purkinje) is unchanged.
        FileFinder file_finder_monodomain_no_purkinje("TestMonodomainPurkinjeProblem_normal/output/SimulationResults_V.dat",RelativeTo::ChasteTestOutput);
        FileFinder file_finder_monodomain_with_purkinje("TestMonodomainPurkinjeProblem_purkinje/output/SimulationResults_V.dat",RelativeTo::ChasteTestOutput);

        TS_ASSERT_EQUALS(file_finder_monodomain_no_purkinje.Exists(), true);
        TS_ASSERT_EQUALS(file_finder_monodomain_with_purkinje.Exists(), true);

        NumericFileComparison comparison(file_finder_monodomain_no_purkinje.GetAbsolutePath(),
                                         file_finder_monodomain_with_purkinje.GetAbsolutePath());
        comparison.CompareFiles();

        // For coverage, call Solve again to extend the solution
        HeartConfig::Instance()->SetSimulationDuration(1.1);
        purkinje_problem.Solve();
    }

    // Solve a Purkinje problem on a branched domain
    void TestBranchedMonodomainPurkinjeProblem()
    {
        //Sets up two problems with fibres of the same length,
        // * one with a single fibre of radius 0.5
        // * one with branched fibres of radius 0.3 & 0.4
        //By equivalent cylinder theory, propagation in the two branched fibres should be the same
        //and equal to the single fibre. This occurs because the cross sectional area of the fibres
        //is equal (0.3^2 + 0.4^2 = 0.5^2)

        // Settings common to both problems
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(1.0);
        HeartConfig::Instance()->SetVisualizeWithVtk();
        HeartConfig::Instance()->SetVisualizeWithParallelVtk();
        ReplicatableVector soln_repl_branched;
        ReplicatableVector soln_repl;

        //branched
        {
            // Need to load the mesh ourselves at present, as AbstractCardiacProblem::Initialise will assume a DistributedTetrahedralMesh.
            std::string mesh_base("mesh/test/data/mixed_dimension_meshes/branched_cable_2D_0_to_1mm_200_elements");
            TrianglesMeshReader<2,2> reader(mesh_base);
            MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
            mesh.ConstructFromMeshReader(reader);

            //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
            //test
            for (MixedDimensionMesh<2, 2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
                 iter != mesh.GetCableElementIteratorEnd();
                 ++iter)
            {
                if((*iter)->GetIndex() < 5) //parent branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() >= 5 && (*iter)->GetIndex() < 10) //first sub branch
                {
                    (*iter)->SetAttribute(0.3);
                }
                else if((*iter)->GetIndex() >= 10 && (*iter)->GetIndex() < 15) //second sub branch
                {
                    (*iter)->SetAttribute(0.4);
                }
            }

            // Set up Purkinje problem
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_purkinje_branched");
            PurkinjeCellFactory cell_factory;
            MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
            purkinje_problem.SetMesh(&mesh);
            purkinje_problem.Initialise();

            // Solve
            purkinje_problem.Solve();
            soln_repl_branched.ReplicatePetscVector(purkinje_problem.GetSolution());
        }


        // Need to load the mesh ourselves at present, as AbstractCardiacProblem::Initialise will assume a DistributedTetrahedralMesh.
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
        //test
        for (MixedDimensionMesh<2, 2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            if((*iter)->GetIndex() < 5) //First half of branch
            {
                (*iter)->SetAttribute(1.0);
            }
            else if((*iter)->GetIndex() >= 5 && (*iter)->GetIndex() < 10) //Second half of branch
            {
                (*iter)->SetAttribute(0.5);
            }
        }

        {
            // Set up Purkinje problem
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_purkinje_branched_symmetric");
            PurkinjeCellFactory cell_factory;
            MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
            purkinje_problem.SetMesh(&mesh);
            purkinje_problem.Initialise();

            // Solve
            purkinje_problem.Solve();
            soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());
        }

#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        //Check for VTK file
        std::string filepath = OutputFileHandler::GetChasteTestOutputDirectory() + "TestMonodomainPurkinjeProblem_purkinje_branched/vtk_output/";

        FileFinder vtk_file(filepath + "SimulationResults.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        //Check for parallel VTK file
        if (PetscTools::IsParallel())
        {
            std::string filepath2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestMonodomainPurkinjeProblem_purkinje_branched_symmetric/vtk_output/";

            FileFinder pvtk_file(filepath2 + "SimulationResults.pvtu", RelativeTo::Absolute);
            TS_ASSERT(pvtk_file.Exists());
        }
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // purkinje nodes for this mesh
            if((55<=i && i<=65))
            {
                //Compare the solution between the single Purkinje cable, the parent and first sub-branch
                //of the branched Purkinje cables.
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], soln_repl[2*i+1], 1e-5);
            }
            else if(i == 71 || i == 82 || i == 93 || i == 104 || i == 115)
            {
                //Compare the solution between the single Purkinje cable and the second sub-branch of the branched solution
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], soln_repl[(((i - 61) / 11) + 61)*2 + 1], 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], 0.0, 1e-4)
            }
        }
    }

    void TestBranchedMonodomainPurkinjeProblemAntiSymmetric()
    {
        //Sets up two problems with branched fibres of different lengths,
        //The radii of the child branches are anti-symmetric and the results compared.
        //The different lengths are required to avoid the branches having the
        //same transmembrane potential (see previous test)

        // Settings common to both problems
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(1.0);

        ReplicatableVector soln_repl_branched;
        ReplicatableVector soln_repl_anti_branched;

        //branched
        {
            // Need to load the mesh ourselves at present, as AbstractCardiacProblem::Initialise will assume a DistributedTetrahedralMesh.
            std::string mesh_base("mesh/test/data/mixed_dimension_meshes/branched_cable_2D_0_to_1mm_200_elements");
            TrianglesMeshReader<2,2> reader(mesh_base);
            MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
            mesh.ConstructFromMeshReader(reader);

            //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
            //test
            for (MixedDimensionMesh<2, 2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
                 iter != mesh.GetCableElementIteratorEnd();
                 ++iter)
            {
                if((*iter)->GetIndex() < 5) //parent branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() >= 5 && (*iter)->GetIndex() < 10) //first sub branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() >= 10 && (*iter)->GetIndex() < 14) //second sub branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() == 14) //Second sub branch is made shorter to test anti-symmetry
                {
                    // set the radius to, essentially, zero. Can't set it to exactly zero as would then
                    // throw exception
                    (*iter)->SetAttribute(2*DBL_EPSILON);
                }
            }

            // Set up Purkinje problem
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblemPurkinjeBranched");
            PurkinjeCellFactory cell_factory;
            MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
            purkinje_problem.SetMesh(&mesh);
            purkinje_problem.Initialise();

            // Solve
            purkinje_problem.Solve();
            soln_repl_branched.ReplicatePetscVector(purkinje_problem.GetSolution());
        }


        // Need to load the mesh ourselves at present, as AbstractCardiacProblem::Initialise will assume a DistributedTetrahedralMesh.
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/branched_cable_2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
        //test
        for (MixedDimensionMesh<2, 2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
                if((*iter)->GetIndex() < 5) //parent branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() >= 5 && (*iter)->GetIndex() < 9) //first sub branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() >= 10 && (*iter)->GetIndex() < 15) //second sub branch
                {
                    (*iter)->SetAttribute(1.0);
                }
                else if((*iter)->GetIndex() == 9) //First sub branch is made shorter to test anti-symmetry
                {
                    // set the radius to, essentially, zero. Can't set it to exactly zero as would then
                    // throw exception
                    (*iter)->SetAttribute(2*DBL_EPSILON);
                }
        }

        {
            // Set up Purkinje problem
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblemPurkinjeBranched");
            PurkinjeCellFactory cell_factory;
            MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
            purkinje_problem.SetMesh(&mesh);
            purkinje_problem.Initialise();

            // Solve
            purkinje_problem.Solve();
            soln_repl_anti_branched.ReplicatePetscVector(purkinje_problem.GetSolution());
        }

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // Purkinje nodes for the parent branch
            if((55<=i && i<=60))
            {
                //Compare the solution between the parent branches of the two solutions
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], soln_repl_anti_branched[2*i+1], 1e-5);
            }
            else if (61 <= i && i <=65) //note,
            {
                //Test the anti-symmetric child branches
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], soln_repl_anti_branched[(((i - 60) * 11) + 60)*2 + 1], 1e-5);
                TS_ASSERT_DELTA(soln_repl_branched[(((i - 60) * 11) + 60)*2 + 1], soln_repl_anti_branched[2*i+1], 1e-5);
            }
            else if(i != 71 && i != 82 && i != 93 && i != 104 && i != 115) //These nodes are on the sub-branch, but their values are tested in the block above
            {
                TS_ASSERT_DELTA(soln_repl_branched[2*i+1], 0.0, 1e-4)
            }
        }
    }

    //Sets up a PVJ stimulus between two independent (not in tissue) cell models
    void TestPVJStimulusTwoCellModels()
    {
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double magnitude_stimulus = -100;  // uA/cm2
        double duration_stimulus = 1;  // ms
        double start_stimulus = 0.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                magnitude_stimulus,
                duration_stimulus,
                start_stimulus));

        //Create PVJ stimuli & multi stimuli
        double pvj_resistance = 10; //kilo Ohms
        boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_ventricular_stim(new PurkinjeVentricularJunctionStimulus(pvj_resistance));
        boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_purkinje_stim(new PurkinjeVentricularJunctionStimulus(pvj_resistance));
        p_pvj_purkinje_stim->SetAppliedToPurkinjeCellModel();

        boost::shared_ptr<MultiStimulus> p_multi_stim_ventricular(new MultiStimulus);
        p_multi_stim_ventricular->AddStimulus(p_pvj_ventricular_stim);

        boost::shared_ptr<MultiStimulus> p_multi_stim_purkinje(new MultiStimulus);
        p_multi_stim_purkinje->AddStimulus(p_stimulus);
        p_multi_stim_purkinje->AddStimulus(p_pvj_purkinje_stim);

        //Create the cell models
        CellLuoRudy1991FromCellML purkinje_cell(p_solver, p_multi_stim_purkinje);
        CellLuoRudy1991FromCellML myocardial_cell(p_solver, p_multi_stim_ventricular);

        //set cell models on pvj stimuli
        p_pvj_ventricular_stim->SetVentricularCellModel(&myocardial_cell);
        p_pvj_ventricular_stim->SetPurkinjeCellModel(&purkinje_cell);
        p_pvj_purkinje_stim->SetVentricularCellModel(&myocardial_cell);
        p_pvj_purkinje_stim->SetPurkinjeCellModel(&purkinje_cell);

        double dummy_time = -100.0;  //Note, stimulus time is ignored
        //Initially the ventricular & Purkinje cells have the same transmembrane potential, thus
        //we expect zero current across the PVJ
        TS_ASSERT_DELTA(p_pvj_ventricular_stim->GetStimulus(dummy_time), 0.0, 1e-6); //Note, stimulus time is ignored
        TS_ASSERT_DELTA(p_pvj_purkinje_stim->GetStimulus(dummy_time), 0.0, 1e-6);

        //Check that a potential difference creates equal and opposite currents
        double rest_voltage = purkinje_cell.GetVoltage();
        purkinje_cell.SetVoltage(50.0);
        TS_ASSERT_DELTA(p_pvj_ventricular_stim->GetStimulus(dummy_time), -13.3853, 1e-4);
        TS_ASSERT_DELTA(p_pvj_ventricular_stim->GetStimulus(dummy_time), -p_pvj_purkinje_stim->GetStimulus(dummy_time), 1e-6);
        purkinje_cell.SetVoltage(rest_voltage);

        double time_step = 0.01;
        for (unsigned i=0; i<500; i++)
        {
            double time = i*time_step;
            myocardial_cell.Compute(time, time+time_step, time_step);
            purkinje_cell.Compute(time, time+time_step, time_step);
        }

        //Purkinje was stimulated
        TS_ASSERT_LESS_THAN(0.0, purkinje_cell.GetVoltage());

        //Myocardium got some of the action
        TS_ASSERT_LESS_THAN(rest_voltage, myocardial_cell.GetVoltage());
        TS_ASSERT_LESS_THAN(0.0, myocardial_cell.GetVoltage());


        // Now try setting up the junction using the cell factory's helper method, instead of doing it manually
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        boost::shared_ptr<AbstractCardiacCell> p_purkinje_cell2(new CellLuoRudy1991FromCellML(p_solver, p_stimulus));
        boost::shared_ptr<AbstractCardiacCell> p_myocardial_cell2(new CellLuoRudy1991FromCellML(p_solver, p_zero_stimulus));
        PurkinjeCellFactory cell_factory;
        cell_factory.CreateJunction(NULL, p_purkinje_cell2.get(), p_myocardial_cell2.get(), pvj_resistance);

        // As above, should have zero current across the junction initially (with no external stimulus)
        TS_ASSERT_DELTA(p_purkinje_cell2->GetStimulus(dummy_time), 0.0, 1e-6);
        TS_ASSERT_DELTA(p_myocardial_cell2->GetStimulus(dummy_time), 0.0, 1e-6);

        // Now simulate and check the results match those above
        for (unsigned i=0; i<500; i++)
        {
            double time = i*time_step;
            p_myocardial_cell2->Compute(time, time+time_step, time_step);
            p_purkinje_cell2->Compute(time, time+time_step, time_step);
        }
        TS_ASSERT_DELTA(purkinje_cell.GetVoltage(), p_purkinje_cell2->GetVoltage(), 1e-12);
        TS_ASSERT_DELTA(myocardial_cell.GetVoltage(), p_myocardial_cell2->GetVoltage(), 1e-12);
    }

    //Solve a retrograde activation problem. A square of myocardium is coupled through a PMJ to a
    //single Purkinje fibre. The myocardium is stimulated and the activation wave should propagate through
    //the PMJ to activate the Purkinje fibre.
    void TestMonodomainPurkinjeRetrogradeActivation()
    {
        // This test slows dramatically with number of processes.  Not know to deadlock, just to run slowly.
        if (PetscTools::GetNumProcs() > 5u)
        {
            TS_TRACE("This test is not suitable for more than 5 processes.");
            return;
        }
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(6.5);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_half_cable");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_retrograde_purkinje");
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        PurkinjeCellFactory cell_factory(false, true);  // false: Only stimulate the myocardial region, not the Purkinje region
                                                        // true: Make a junction at (0.05, 0.05)


        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        // Solve
        purkinje_problem.Solve();

        //Replicate the solution for easy access
        ReplicatableVector soln_repl;
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        unsigned purkinje_edge_node_index = 55u;
        unsigned purkinje_junction_node_index = 60u;
        double purkinje_edge_voltage = soln_repl[2*purkinje_edge_node_index+1];
        double purkinje_junction_voltage = soln_repl[2*purkinje_junction_node_index+1];
        ///\todo #1900 The junction is doing *something*, but it's only very slight
        TS_ASSERT_LESS_THAN(0.0, purkinje_junction_voltage); //Junction got stimulus
        TS_ASSERT_LESS_THAN(0.0, purkinje_edge_voltage); //Wave travelled back along cable

    }

    //Solve an activation problem. A square of myocardium is coupled through several PVJs to a
    //single Purkinje fibre with a manually drawn forward-star at the end. In this case the Purkinje system
    //isn't stimulated, to check that the PVJ isn't erroneously generating current.
    void TestMonodomainPurkinjeNoActivation()
    {
        // This test slows dramatically with number of processes.  Not know to deadlock, just to run slowly.
        if (PetscTools::GetNumProcs() > 5u)
        {
            TS_TRACE("This test is not suitable for more than 5 processes.");
            return;
        }
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_noactivate_purkinje");
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        PurkinjeStarCellFactory cell_factory(false, true); //no stimulus, junction
        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        //With no stimulus & the same cell model there should be no difference between the Purkinje and the myocardium.
        //However, the PMJ subtly magnifies numerical errors so that the Purkinje & Myocardial voltages drift apart slightly
        //Here we check that they do not drift apart during a longer running simulation
        HeartConfig::Instance()->SetSimulationDuration(50.0);
        purkinje_problem.Solve();

        ReplicatableVector soln_repl;
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        unsigned purkinje_junction_node_index = 60u;
        double purkinje_junction_voltage = soln_repl[2*purkinje_junction_node_index+1];
        double monodomain_junction_voltage = soln_repl[2*purkinje_junction_node_index];

        TS_ASSERT_DELTA(monodomain_junction_voltage, purkinje_junction_voltage, 1e-4);
    }


    //Solve an activation problem. A square of myocardium is coupled through several PVJs to a
    //single Purkinje fibre with a manually drawn forward-star at the end. The Purkinje system is stimulated and the activation wave should
    //propagate through to the tissue.
    void TestMonodomainPurkinjeActivationViaStar()
    {
        // This test slows dramatically with number of processes.  Not know to deadlock, just to run slowly.
        if (PetscTools::GetNumProcs() > 5u)
        {
            TS_TRACE("This test is not suitable for more than 5 processes.");
            return;
        }
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_activate_purkinje");
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        PurkinjeStarCellFactory cell_factory; //Stimulates at the left (only on Purkinje cable) and makes multiple junctions
        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        // Solve to the point where the stimulus enters the PVJ
        HeartConfig::Instance()->SetSimulationDuration(1.8);
        purkinje_problem.Solve();

        ReplicatableVector soln_repl;
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        //Check that during the upstroke the myocardium is behind the Purkinje
        unsigned purkinje_junction_node_index = 60u;
        double purkinje_junction_voltage = soln_repl[2*purkinje_junction_node_index+1];
        double monodomain_junction_voltage = soln_repl[2*purkinje_junction_node_index];

        TS_ASSERT_LESS_THAN(monodomain_junction_voltage, purkinje_junction_voltage);

        //Solve so that the whole myocardium is activated
        HeartConfig::Instance()->SetSimulationDuration(6.5);
        purkinje_problem.Solve();

        //Check that the myocardium has been activated.
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());
        unsigned monodomain_corner_index = 0u;
        double monodomain_corner_voltage = soln_repl[2*monodomain_corner_index];
        TS_ASSERT_LESS_THAN(10.0, monodomain_corner_voltage);
    }

    //Solve an activation problem. A square of myocardium is coupled through several PVJs to a
    //single Purkinje fibre with a manually drawn forward-star at the end. The Purkinje system is stimulated and the activation wave should
    //propagate through to the tissue. The PVJ nodes are defined in a file.
    void TestMonodomainPurkinjeActivationFromFile()
    {
        // This test slows dramatically with number of processes.  Not know to deadlock, just to run slowly.
        if (PetscTools::GetNumProcs() > 5u)
        {
            TS_TRACE("This test is not suitable for more than 5 processes.");
            return;
        }
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_activate_purkinje");
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        PurkinjeStarCellFactoryFromFile cell_factory; //Stimulates at the left (only on Purkinje cable) and makes multiple junctions

        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        // Solve to the point where the stimulus enters the PVJ
        HeartConfig::Instance()->SetSimulationDuration(1.8);
        purkinje_problem.Solve();

        ReplicatableVector soln_repl;
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        //Check that during the upstroke the myocardium is behind the Purkinje
        unsigned purkinje_junction_node_index = 60u;
        double purkinje_junction_voltage = soln_repl[2*purkinje_junction_node_index+1];
        double monodomain_junction_voltage = soln_repl[2*purkinje_junction_node_index];

        TS_ASSERT_LESS_THAN(monodomain_junction_voltage, purkinje_junction_voltage);

        //Solve so that the whole myocardium is activated
        HeartConfig::Instance()->SetSimulationDuration(6.5);
        purkinje_problem.Solve();

        //Check that the myocardium has been activated.
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());
        unsigned monodomain_corner_index = 0u;
        double monodomain_corner_voltage = soln_repl[2*monodomain_corner_index];
        TS_ASSERT_LESS_THAN(10.0, monodomain_corner_voltage);
    }

    //Solve the same activation problem as  TestMonodomainPurkinjeActivationFromFile with a non-dumb partitioning
    void TestMonodomainPurkinjeActivationFromFileWithPartitioning()
    {
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star");
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_activate_purkinje");

        ///Note that if PETSc partitioning is not available then
        // * We drop through and create a extra warning
        //HeartConfig::Instance()->SetMeshPartitioning("petsc");
        HeartConfig::Instance()->SetMeshPartitioning("parmetis");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        PurkinjeStarCellFactoryFromFile cell_factory; //Stimulates at the left (only on Purkinje cable) and makes multiple junctions

        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        purkinje_problem.Initialise();

        // Solve to the point where the stimulus enters the PVJ
        HeartConfig::Instance()->SetSimulationDuration(1.8);
        purkinje_problem.Solve();

        ReplicatableVector soln_repl;
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

        //Check that during the upstroke the myocardium is behind the Purkinje
        unsigned purkinje_junction_node_index = 60u;
        double purkinje_junction_voltage = soln_repl[2*purkinje_junction_node_index+1];
        double monodomain_junction_voltage = soln_repl[2*purkinje_junction_node_index];

        TS_ASSERT_LESS_THAN(monodomain_junction_voltage, purkinje_junction_voltage);

        //Solve so that the whole myocardium is activated
        HeartConfig::Instance()->SetSimulationDuration(6.5);
        purkinje_problem.Solve();

        //Check that the myocardium has been activated.
        soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());
        unsigned monodomain_corner_index = 0u;
        double monodomain_corner_voltage = soln_repl[2*monodomain_corner_index];
        TS_ASSERT_LESS_THAN(10.0, monodomain_corner_voltage);
    }

    void TestPvjFileErrorHandling()
    {
        // Set up Purkinje problem
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinjeProblem_activate_purkinje");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshPartitioning(), DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);

        //Cover warning that a .pvj file doesn't exist
        PurkinjeStarCellFactoryFromFile cell_factory; //Stimulates at the left (only on Purkinje cable) and makes multiple junctions
        MonodomainPurkinjeProblem<2,2> purkinje_problem(&cell_factory);
        {
            Warnings::QuietDestroy();
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
            purkinje_problem.Initialise();
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "No Purkinje-Ventricular junction (.pvj) file found. Junctions must be specified manually.");
            Warnings::QuietDestroy();
        }

        //Cover warning that a .pvj file can't be read
        if (PetscTools::AmMaster())
        {
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star");

            chmod("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star.pvj", 0000);
            TS_ASSERT_THROWS_CONTAINS(purkinje_problem.Initialise(), "Couldn't open data file:");
            chmod("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements_cable_and_star.pvj", CHASTE_READ_WRITE);
        }

        //Covers empty lines in .pvj file
        {
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/mixed_dimension_meshes/bad_pvj_mesh");

            TS_ASSERT_THROWS_CONTAINS(purkinje_problem.Initialise(), "No data found on line in file:");
        }
    }
};

#endif // TESTMONODOMAINPURKINJEPROBLEM_HPP_
