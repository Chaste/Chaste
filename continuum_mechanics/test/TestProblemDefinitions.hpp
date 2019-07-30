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


#ifndef TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_
#define TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_


#include <cxxtest/TestSuite.h>
#include "SolidMechanicsProblemDefinition.hpp"
#include "StokesFlowProblemDefinition.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"

#include "PetscSetupAndFinalize.hpp"

c_vector<double,2> SomeFunction(c_vector<double,2>& rX, double t)
{
    c_vector<double,2> body_force;
    body_force(0) = rX(0)+t;
    body_force(1) = 2*(rX(1)+t);
    return body_force;
}

c_vector<double,2> AnotherFunction(c_vector<double,2>& rX, double t)
{
    c_vector<double,2> body_force;
    body_force(0) = rX(0)*t;
    body_force(1) = 10*rX(1)*t;
    return body_force;
}

double PressureFunction(double t)
{
    return 3*t;
}

// Helper function for creating std vectors. Ugly, but makes things below much neater (avoids many, many push_backs).
template<class T>
std::vector<T> MakeStdVec(unsigned size, T value0=0, T value1=0, T value2=0, T value3=0, T value4=0)
{
    std::vector<T> ret(size);
    ret[0] = value0;
    if (size > 1)
    {
        ret[1] = value1;
    }
    if (size > 2)
    {
        ret[2] = value2;
    }
    if (size > 3)
    {
        ret[3] = value3;
    }
    if (size > 4)
    {
        ret[4] = value4;
    }
    assert(size <= 5);
    return ret;
}

class TestProblemDefinitions : public CxxTest::TestSuite
{
public:
    // Test all the functionality inside ContinuumMechanicsProblemDefinition,
    // which will be common to other problem definition classes
    void TestContinuumMechanicsProblemDefinition()
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);

        ContinuumMechanicsProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_THROWS_THIS(problem_defn.Validate(), "No Dirichlet boundary conditions (eg fixed displacement or fixed flow) have been set");

        TS_ASSERT_DELTA(problem_defn.GetDensity(), 1.0, 1e-12);

        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), CONSTANT_BODY_FORCE);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(1), 0.0, 1e-12);

        TS_ASSERT_EQUALS(problem_defn.GetTractionBoundaryConditionType(), NO_TRACTIONS);

        problem_defn.SetDensity(2.0);
        TS_ASSERT_DELTA(problem_defn.GetDensity(), 2.0, 1e-12);


        //////////////////////////////////
        // Body force
        //////////////////////////////////

        problem_defn.SetBodyForce(SomeFunction);
        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), FUNCTIONAL_BODY_FORCE);
        c_vector<double,2> X;
        X(0) = 10.0;
        X(1) = 11.0;
        double t = 0.5;
        TS_ASSERT_DELTA(problem_defn.EvaluateBodyForceFunction(X,t)(0), 10.5, 1e-12);
        TS_ASSERT_DELTA(problem_defn.EvaluateBodyForceFunction(X,t)(1), 23.0, 1e-12);

        TS_ASSERT_DELTA(problem_defn.GetBodyForce(X,t)(0), 10.5, 1e-12);
        TS_ASSERT_DELTA(problem_defn.GetBodyForce(X,t)(1), 23.0, 1e-12);

        c_vector<double,2> body_force;
        body_force(0) = -9.81;
        body_force(1) = 0.01;
        problem_defn.SetBodyForce(body_force);
        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), CONSTANT_BODY_FORCE);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(0), -9.81,  1e-12);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(1),  0.01, 1e-12);

        TS_ASSERT_DELTA(problem_defn.GetBodyForce(X,t)(0), -9.81, 1e-12);
        TS_ASSERT_DELTA(problem_defn.GetBodyForce(X,t)(1),  0.01, 1e-12);


        //////////////////////////////////
        // Traction
        //////////////////////////////////

        std::vector<BoundaryElement<1,2>*> boundary_elements;
        std::vector<c_vector<double,2> > tractions;

        TetrahedralMesh<2,2>::BoundaryElementIterator iter
           = mesh.GetBoundaryElementIteratorBegin();

        c_vector<double,2> vec = zero_vector<double>(2);
        vec(0)=1.0;
        boundary_elements.push_back(*iter);
        tractions.push_back(vec);

        ++iter;
        vec(1)=2.0;
        boundary_elements.push_back(*iter);
        tractions.push_back(vec);

        problem_defn.SetTractionBoundaryConditions(boundary_elements, tractions);

        TS_ASSERT_EQUALS(problem_defn.GetTractionBoundaryConditionType(), ELEMENTWISE_TRACTION);

        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 2u);
        TS_ASSERT_EQUALS(problem_defn.rGetElementwiseTractions().size(), 2u);

        // comparing addresses
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[0], boundary_elements[0]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[1], boundary_elements[1]);

        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[0](0), 1.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[1](0), 1.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[1](1), 2.0, 1e-12);

        ++iter;
        boundary_elements.push_back(*iter);
        double pressure = 3423.342;

        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elements, pressure);

        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 3u);

        // comparing addresses
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[0], boundary_elements[0]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[1], boundary_elements[1]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[2], boundary_elements[2]);

        TS_ASSERT_DELTA(problem_defn.GetNormalPressure(), 3423.342, 1e-12);

        problem_defn.SetPressureScaling(10.0/43.0);
        TS_ASSERT_DELTA(problem_defn.GetNormalPressure(), 3423.342*10/43, 1e-12);

        problem_defn.SetPressureScaling(1);


        ++iter;
        boundary_elements.push_back(*iter);
        problem_defn.SetTractionBoundaryConditions(boundary_elements, AnotherFunction);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 4u);

        TS_ASSERT_DELTA(problem_defn.EvaluateTractionFunction(X,t)(0), 5.0,  1e-12);
        TS_ASSERT_DELTA(problem_defn.EvaluateTractionFunction(X,t)(1), 55.0, 1e-12);

        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elements, PressureFunction);
        TS_ASSERT_EQUALS(problem_defn.GetTractionBoundaryConditionType(), FUNCTIONAL_PRESSURE_ON_DEFORMED);
        TS_ASSERT_DELTA(problem_defn.EvaluateNormalPressureFunction(3.64455), 3*3.64455, 1e-12);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        problem_defn.SetZeroDirichletNodes(fixed_nodes); // note, this functionality is all tested properly below

        // should not throw anything
        problem_defn.Validate();
    }

    // Test the functionality specific to SolidMechanicsProblemDefinition
    void TestSolidMechanicsProblemDefinition()
    {
        TS_ASSERT_EQUALS(SolidMechanicsProblemDefinition<2>::FREE, DBL_MAX);
        TS_ASSERT_LESS_THAN(0, SolidMechanicsProblemDefinition<2>::FREE);

        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        //////////////////////////////////
        // Fixed nodes
        //////////////////////////////////

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        fixed_nodes.push_back(4);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes().size(), 2u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodeValues().size(), 2u);

        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](1), 0.0, 1e-12);

        fixed_nodes.push_back(8);
        fixed_nodes.push_back(9);
        fixed_nodes.push_back(10);


        std::vector<c_vector<double,2> > locations;
        c_vector<double,2> location = zero_vector<double>(2);
        // Node 0 is to be placed at (0,0)
        locations.push_back(location);

        // Node 4 is to be placed at (0,0.1)
        location(1)=0.1;
        locations.push_back(location);

        // Node 8 is to be placed at (0.1,0.1)
        location(0)=0.1;
        locations.push_back(location);

        // Node 9 is to be placed at (0.5,FREE)
        location(0) = 0.5;
        location(1) = SolidMechanicsProblemDefinition<2>::FREE;
        locations.push_back(location);

        // Node 9 is to be placed at (FREE,1.5)
        location(0) = SolidMechanicsProblemDefinition<2>::FREE;
        location(1) = 1.5;
        locations.push_back(location);

        problem_defn.SetFixedNodes(fixed_nodes, locations);

        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes().size(), 5u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodeValues().size(), 5u);

        // the fully fixed nodes
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[2], 8u);

        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](0), 0.0 - mesh.GetNode(0)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](1), 0.0 - mesh.GetNode(0)->rGetLocation()[1], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](0), 0.0 - mesh.GetNode(4)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](1), 0.1 - mesh.GetNode(4)->rGetLocation()[1], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[2](0), 0.1 - mesh.GetNode(8)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[2](1), 0.1 - mesh.GetNode(8)->rGetLocation()[1], 1e-12);

        // the partial fixed nodes
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[3], 9u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[4], 10u);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[3](0), 0.5 - mesh.GetNode(9)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[3](1), SolidMechanicsProblemDefinition<2>::FREE, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[4](0), SolidMechanicsProblemDefinition<2>::FREE, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[4](1), 1.5 - mesh.GetNode(10)->rGetLocation()[1], 1e-12);


        ///////////////////////////////////////
        // Set an incompressible material law
        ///////////////////////////////////////
        TS_ASSERT_THROWS_THIS(problem_defn.Validate(), "No material law has been set");

        // set a homogeneous law
        MooneyRivlinMaterialLaw<2> incomp_mooney_rivlin_law(1.0);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&incomp_mooney_rivlin_law);

        TS_ASSERT_EQUALS(problem_defn.IsHomogeneousMaterial(), true);
        TS_ASSERT_EQUALS(problem_defn.GetCompressibilityType(), INCOMPRESSIBLE);
        TS_ASSERT_EQUALS(problem_defn.GetIncompressibleMaterialLaw(0), &incomp_mooney_rivlin_law);

        // set a heterogeneous law
        MooneyRivlinMaterialLaw<2> incomp_mooney_rivlin_law_2(2.0);
        std::vector<AbstractMaterialLaw<2>*> laws;
        for (unsigned i=0; i<mesh.GetNumElements()/2; i++)
        {
            laws.push_back(&incomp_mooney_rivlin_law);
        }
        for (unsigned i=mesh.GetNumElements()/2; i<mesh.GetNumElements(); i++)
        {
            laws.push_back(&incomp_mooney_rivlin_law_2);
        }

        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,laws);

        TS_ASSERT_EQUALS(problem_defn.IsHomogeneousMaterial(), false);
        for (unsigned i=0; i<mesh.GetNumElements()/2; i++)
        {
            TS_ASSERT_EQUALS(problem_defn.GetIncompressibleMaterialLaw(i), &incomp_mooney_rivlin_law);
        }
        for (unsigned i=mesh.GetNumElements()/2; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(problem_defn.GetIncompressibleMaterialLaw(i), &incomp_mooney_rivlin_law_2);
        }

        /////////////////////////////////////////////////////////////////////////
        // Set a compressible material law (clears the incompressible laws)
        /////////////////////////////////////////////////////////////////////////

        CompressibleMooneyRivlinMaterialLaw<2> comp_mooney_rivlin_law(2.0, 1.0);
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&comp_mooney_rivlin_law);

        TS_ASSERT_EQUALS(problem_defn.IsHomogeneousMaterial(), true);
        TS_ASSERT_EQUALS(problem_defn.GetCompressibilityType(), COMPRESSIBLE);
        TS_ASSERT_EQUALS(problem_defn.GetCompressibleMaterialLaw(0), &comp_mooney_rivlin_law);

        // set a heterogeneous law
        CompressibleMooneyRivlinMaterialLaw<2> comp_mooney_rivlin_law_2(4.0, 1.0);
        std::vector<AbstractMaterialLaw<2>*> comp_laws;
        for (unsigned i=0; i<mesh.GetNumElements()/2; i++)
        {
            comp_laws.push_back(&comp_mooney_rivlin_law);
        }
        for (unsigned i=mesh.GetNumElements()/2; i<mesh.GetNumElements(); i++)
        {
            comp_laws.push_back(&comp_mooney_rivlin_law_2);
        }

        problem_defn.SetMaterialLaw(COMPRESSIBLE,comp_laws);

        TS_ASSERT_EQUALS(problem_defn.IsHomogeneousMaterial(), false);
        for (unsigned i=0; i<mesh.GetNumElements()/2; i++)
        {
            TS_ASSERT_EQUALS(problem_defn.GetCompressibleMaterialLaw(i), &comp_mooney_rivlin_law);
        }
        for (unsigned i=mesh.GetNumElements()/2; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(problem_defn.GetCompressibleMaterialLaw(i), &comp_mooney_rivlin_law_2);
        }

        // should not throw anything
        problem_defn.Validate();

        TS_ASSERT_THROWS_THIS(problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&comp_mooney_rivlin_law),"Compressibility type was declared as INCOMPRESSIBLE but a compressible material law was given");
        TS_ASSERT_THROWS_THIS(problem_defn.SetMaterialLaw(COMPRESSIBLE,&incomp_mooney_rivlin_law),"Incompressibility type was declared as COMPRESSIBLE but an incompressible material law was given");

        ///////////////////////////////
        // solver stuff
        ///////////////////////////////
        TS_ASSERT_EQUALS(problem_defn.GetSolveUsingSnes(), false);
        TS_ASSERT_EQUALS(problem_defn.GetVerboseDuringSolve(), false);

        problem_defn.SetSolveUsingSnes();
        problem_defn.SetVerboseDuringSolve();

        TS_ASSERT_EQUALS(problem_defn.GetSolveUsingSnes(), true);
        TS_ASSERT_EQUALS(problem_defn.GetVerboseDuringSolve(), true);

        problem_defn.SetSolveUsingSnes(false);
        problem_defn.SetVerboseDuringSolve(false);

        TS_ASSERT_EQUALS(problem_defn.GetSolveUsingSnes(), false);
        TS_ASSERT_EQUALS(problem_defn.GetVerboseDuringSolve(), false);

    }

    void TestStokesFlowProblemDefinition()
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);

        StokesFlowProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_THROWS_THIS(problem_defn.GetViscosity(), "Viscosity hasn't been set yet (for the Stokes' flow problem)");

        problem_defn.SetViscosity(1.3423423);
        TS_ASSERT_EQUALS(problem_defn.GetViscosity(),1.3423423);

        //////////////////////////////////
        // Fixed nodes
        //////////////////////////////////

        std::vector<unsigned> fixed_flow_nodes;
        fixed_flow_nodes.push_back(0);
        fixed_flow_nodes.push_back(4);
        problem_defn.SetZeroFlowNodes(fixed_flow_nodes);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes().size(), 2u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodeValues().size(), 2u);

        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](1), 0.0, 1e-12);

        fixed_flow_nodes.push_back(8);
        fixed_flow_nodes.push_back(9);
        fixed_flow_nodes.push_back(10);

        std::vector<c_vector<double,2> > fixed_flows;
        c_vector<double,2> flow = zero_vector<double>(2);
        fixed_flows.push_back(flow);
        flow(1)=0.1;
        fixed_flows.push_back(flow);
        flow(0)=0.1;
        fixed_flows.push_back(flow);

        flow(0) = 0.5;
        flow(1) = StokesFlowProblemDefinition<2>::FREE;
        fixed_flows.push_back(flow);

        flow(0) = StokesFlowProblemDefinition<2>::FREE;
        flow(1) = 1.5;
        fixed_flows.push_back(flow);

        problem_defn.SetPrescribedFlowNodes(fixed_flow_nodes, fixed_flows);

        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes().size(), 5u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodeValues().size(), 5u);

        // the fully fixed nodes
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[2], 8u);

        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[1](1), 0.1, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[2](0), 0.1, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[2](1), 0.1, 1e-12);

        // the partial fixed nodes
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[3], 9u);
        TS_ASSERT_EQUALS(problem_defn.rGetDirichletNodes()[4], 10u);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[3](0), 0.5, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[3](1), StokesFlowProblemDefinition<2>::FREE, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[4](0), StokesFlowProblemDefinition<2>::FREE, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetDirichletNodeValues()[4](1), 1.5, 1e-12);

        // should not throw anything
        problem_defn.Validate();
    }
};

#endif /* TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_ */
