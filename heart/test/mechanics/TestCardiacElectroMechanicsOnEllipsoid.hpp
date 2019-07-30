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

#ifndef TESTONELLIPSOID_HPP_
#define TESTONELLIPSOID_HPP_

#include <cxxtest/TestSuite.h>
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
#include "SimpleStimulus.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"

class ApexStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    double mZvalueToStimulateBelow;

public:
    ApexStimulusCellFactory(double zValueToStimulateBelow)
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1000.0*200, 0.5)),
          mZvalueToStimulateBelow(zValueToStimulateBelow)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        // Stimulate the apex
        if (pNode->rGetLocation()[2] < mZvalueToStimulateBelow)
        {
            return new CellLuoRudy1991FromCellML(mpSolver,mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver,mpZeroStimulus);
        }
    }
};


// Obviously only run this test with an optimised build such as build=GccOpt_ndebug.
//
// To watch progress, run from the command line with -mech_very_verbose -mesh_pair_verbose

class TestCardiacElectroMechanicsOnEllipsoid : public CxxTest::TestSuite
{
public:
    void TestOnEllipsoid()
    {
        EXIT_IF_PARALLEL; // needs investigation, possibly related to #2057

        /////////////////////////////////////////////////////////////
        //
        //   Read the meshes for electrics and mechanics.
        //   Also output mesh resolutions
        //
        //   Note: the electrics mesh is a bit too coarse to expect
        //   accurate propagation velocities
        //
        //
        /////////////////////////////////////////////////////////////
        TetrahedralMesh<3,3> electrics_mesh;
        QuadraticMesh<3> mechanics_mesh;

        {
            TrianglesMeshReader<3,3> reader1("mesh/test/data/ellipsoid_15811_elements");
            electrics_mesh.ConstructFromMeshReader(reader1);

            TrianglesMeshReader<3,3> reader2("mesh/test/data/ellipsoid_8225_elements_quad",
                                             2 /*quadratic elements*/,
                                             2 /*quadratic boundary elements*/);
            mechanics_mesh.ConstructFromMeshReader(reader2);
        }

        c_vector<double, 2> edge_length_electrics = electrics_mesh.CalculateMinMaxEdgeLengths();
        c_vector<double, 2> edge_length_mechanics = mechanics_mesh.CalculateMinMaxEdgeLengths();
        std::cout << "Typical edge length of electrics and mechanics meshes are, respectively: "
                  << (edge_length_electrics[0] + edge_length_electrics[1])/2.0 << " " <<
                     (edge_length_mechanics[0] + edge_length_mechanics[1])/2.0<< "\n";


        ////////////////////////////////////////////////////////////////
        //
        // Dirichlet boundary conditions: fix base in Z-direction,
        // which one point fixed in all directions to remove rotations
        //
        ////////////////////////////////////////////////////////////////
        double base_threshold = 0.0;
        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,3> > locations;

        bool first = true;
        for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double x = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
            double z = mechanics_mesh.GetNode(i)->rGetLocation()[2];
            if (z >= base_threshold)
            {
                fixed_nodes.push_back(i);
                c_vector<double,3> new_location;

                if (first)
                {
                    new_location(0) = x;
                    new_location(1) = y;
                    first = false;
                }
                else
                {
                    new_location(0) = SolidMechanicsProblemDefinition<3>::FREE;
                    new_location(1) = SolidMechanicsProblemDefinition<3>::FREE;
                }

                new_location(2)= z;
                locations.push_back(new_location);
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        //
        // Find the boundary elements on the endocardium (using knowledge on
        // the geometry)
        //
        ////////////////////////////////////////////////////////////////////////////////

        double a_endo = 0.25;
        double b_endo = 0.35;
        double c_endo = 0.85;

        std::vector<BoundaryElement<2,3>*> boundary_elems;
        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
                = mechanics_mesh.GetBoundaryElementIteratorBegin();
              iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
              ++iter)
        {
            ChastePoint<3> centroid = (*iter)->CalculateCentroid();
            if (centroid[2]<0)
            {
                Node<3>* p_node = (*iter)->GetNode(0);
                double x = p_node->rGetLocation()[0];
                double y = p_node->rGetLocation()[1];
                double z = p_node->rGetLocation()[2];

                double v = sqrt(x*x/(a_endo*a_endo) + y*y/(b_endo*b_endo) + (z*z)/(c_endo*c_endo));

                if (fabs(v-1.0)<0.2 && z > -0.9)
                {
                    BoundaryElement<2,3>* p_element = *iter;
                    boundary_elems.push_back(p_element);
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        //
        // Other config
        //
        ////////////////////////////////////////////////////////////////////////////////

        double apex_stim_threshold = -0.95;
        ApexStimulusCellFactory cell_factory(apex_stim_threshold);

        HeartConfig::Instance()->SetOdeTimeStep(0.005);
        HeartConfig::Instance()->SetSimulationDuration(70.0); // see comment below on mechanics timestep

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.1);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);

        problem_defn.SetMechanicsSolveTimestep(1.0);   // runs through contraction, fails during repolarisation at about 300ms
        //problem_defn.SetMechanicsSolveTimestep(0.1); // runs through entire cycle but takes very very long time - todo: adaptive mechanics timestep
        problem_defn.SetSolveUsingSnes();

        problem_defn.SetVerboseDuringSolve();

        problem_defn.SetFixedNodes(fixed_nodes,locations);
        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, -1.0);
        problem_defn.SetNumIncrementsForInitialDeformation(7);

        FileFinder finder("heart/test/data/fibre_tests/ellipsoid_8225_elements.ortho",RelativeTo::ChasteSourceRoot);
        problem_defn.SetVariableFibreSheetDirectionsFile(finder, false);

        ////////////////////////////////////////////////////////////////////////////////
        //
        // Solve
        //
        ////////////////////////////////////////////////////////////////////////////////
        CardiacElectroMechanicsProblem<3,1> problem(COMPRESSIBLE,
                                                  MONODOMAIN,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacElectroMechanicsEllipsoid");

        // problem.SetNoElectricsOutput();

        problem.Solve();
    }
};

#endif // TESTONELLIPSOID_HPP_
