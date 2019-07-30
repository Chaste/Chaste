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

#ifndef TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_
#define TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_


#include <cxxtest/TestSuite.h>

#include "LabelBasedContractionCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "FakeBathContractionModel.hpp"
#include "ElectroMechanicsProblemDefinition.hpp"
#include "AbstractContractionCellFactory.hpp"

template <unsigned DIM>
class ExampleContractionCellFactory : public AbstractContractionCellFactory<DIM>
{
public:
    ExampleContractionCellFactory()
     : AbstractContractionCellFactory<DIM>()
    {};

    AbstractContractionModel* CreateContractionCellForElement(Element<DIM, DIM>* pElement)
    {
        AbstractContractionModel* p_model;

        c_vector<double, 2> centroid = pElement->CalculateCentroid();

        if (centroid[0] >= 0.07)
        {
            p_model = new FakeBathContractionModel;
        }
        else
        {
            p_model = new Kerchoffs2003ContractionModel;
            if (centroid[1] >= 0.05)
            {
                static_cast<AbstractOdeBasedContractionModel*>(p_model)->SetParameter("tr", 25.0);
                static_cast<AbstractOdeBasedContractionModel*>(p_model)->SetParameter("td", 25.0); // ms
                static_cast<AbstractOdeBasedContractionModel*>(p_model)->SetParameter("b", 75.0); // ms/um
            }
        }
        return p_model;
    }
};

class TestAbstractContractionCellFactory : public CxxTest::TestSuite
{
public:
    void TestContractionCellFactory()
    {
        {
            LabelBasedContractionCellFactory<2> factory(CONSTANT);

            ConstantActiveTension* p_model =
                dynamic_cast<ConstantActiveTension*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(NONPHYSIOL1);

            NonPhysiologicalContractionModel* p_model =
                dynamic_cast<NonPhysiologicalContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(NASH2004);

            Nash2004ContractionModel* p_model =
                dynamic_cast<Nash2004ContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(KERCHOFFS2003);

            Kerchoffs2003ContractionModel* p_model =
                dynamic_cast<Kerchoffs2003ContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }
    }

    void TestContractionCellFactoryOnSquare()
    {

        //2D square meshes for electrics and mechanics
        // create electrics mesh
        TetrahedralMesh<2,2>* p_mesh_e = new TetrahedralMesh<2,2>();
        p_mesh_e->ConstructRegularSlabMesh(0.01, 0.1, 0.1);//width/no of ele, width, height

        // create mechanics mesh
        QuadraticMesh<2>* p_mesh_m = new QuadraticMesh<2>(0.02, 0.1, 0.1);//width/no of ele, width, height

        //General simulation setup
        HeartConfig::Instance()->SetSimulationDuration(5.0);
        //Set up the electrophysiology
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000); //stimulates along X=0

        //Set up the mechanics
        ElectroMechanicsProblemDefinition<2> problem_defn(*p_mesh_m);
        //fixed nodes along X=0
        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(*p_mesh_m, 0, 0.0); // all the X=0.0 nodes

        /*
         * HOW_TO_TAG Cardiac/Electro-mechanics
         * Set heterogeneous contraction models by using a contraction cell factory.
         *
         * Here we give a factory to put a contraction model at each element
         * (this is called by AbstractCardiacMechanicsProblem::Initialise()).
         */
        ExampleContractionCellFactory<2> mechanics_factory;
        problem_defn.SetContractionCellFactory(&mechanics_factory);
        problem_defn.SetContractionModelOdeTimestep(0.01); // NB you need to use the finest timestep of any model the factory will create.

        // The following line is an example of the 'old' way of setting a homogeneous contraction model.
        // This method now just creates a LabelBasedCellFactory as per above...
        //problem_defn.SetContractionModel(KERCHOFFS2003,0.01/*contraction model ODE timestep*/);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);
        //We are adding the mechano-electric coupling to ensure stretches are calculated for testing.
        problem_defn.SetDeformationAffectsElectrophysiology(false,true);
        //Solve the problem
        //2 is space dim, 1 is monodomain (2 would be bidomain)
        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    p_mesh_e,
                                                    p_mesh_m,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestContractionCellFactoryOnSquare");

        problem.Solve();

        // We get the stretches directly as this test is a friend of CardiacElectroMechanicsProblem.
        std::vector<double> stretches = problem.mStretchesForEachMechanicsElement;

        for (unsigned i=0; i<stretches.size(); i++)
        {
            c_vector<double, 2> centroid = p_mesh_m->GetElement(i)->CalculateCentroid();

            //std::cout << "Element " << i << ": pos = [" << centroid[0] << "," << centroid[1] << "], stretch = " << stretches[i] << "\n";

            if (centroid[0] > 0.09) // We set a Fake contraction model near x = 0.1.
            {
                // Slight increases in stretch in this region
                TS_ASSERT_LESS_THAN(0.999, stretches[i]); // 0.99 if run for 50ms.
            }
            else if (centroid[0] > 0.02 && centroid[0] < 0.07)
            {
                // contraction in this region.
                TS_ASSERT_LESS_THAN(stretches[i], 0.998); // 0.9 if run for 50ms.
            }
            // In other regions (near x=0 and x=0.8) you can get nothing or contraction
            // depending on whether you are near the sides (y=0, y=0.1).

            // If this ever fails, then run the test for 50 milliseconds.
            // The top tolerance was 0.99, the bottom one was 0.9
        }

        delete p_mesh_e;
        delete p_mesh_m;
    }
};

#endif /* TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_ */
