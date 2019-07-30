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


#ifndef _TESTELECTROMECHANICSEXACTSOLUTION_HPP_
#define _TESTELECTROMECHANICSEXACTSOLUTION_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/assign.hpp>

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "ExplicitCardiacMechanicsSolver.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "ConstantActiveTension.hpp"


double MATERIAL_PARAM = 2.0;
double ACTIVE_TENSION = 1.5;
double ALPHA = 0.2;

// Body force (when active tension is constant and fibres in X-direction)
// corresponding to the deformation x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
c_vector<double,2> MyBodyForce(c_vector<double,2>& rX, double t)
{
    assert(rX(0)>=0 && rX(0)<=1 && rX(1)>=0 && rX(1)<=1);

    c_vector<double,2> body_force;

    double lam = 1+ALPHA*rX(0);
    double mu = rX(1)*ALPHA/(lam*lam);
    double denom_sqd = pow(lam*lam + mu*mu,2);
    double aY = rX(1)*ALPHA;

    body_force(0) = -2*MATERIAL_PARAM * ALPHA                                - ACTIVE_TENSION*ALPHA*(5*mu*mu -lam*lam)/denom_sqd;
    body_force(1) = -2*MATERIAL_PARAM * 2*ALPHA*ALPHA*rX(1)/(lam*lam*lam)    - ACTIVE_TENSION*ALPHA*aY*(4.0/lam - 2*aY*aY*pow(lam,-7))/denom_sqd;

    return body_force;
}

// Surface traction (when active tension is constant and fibres in X-direction)
// on three sides of a cube, corresponding to x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
c_vector<double,2> MyTraction(c_vector<double,2>& location, double t)
{
    c_vector<double,2> traction = zero_vector<double>(2);

    double lam = 1+ALPHA*location(0);
    double mu = location(1)*ALPHA/(lam*lam);
    double denom = lam*lam + mu*mu;

    if (fabs(location(0)-1.0) <= DBL_EPSILON) //Right edge
    {
        traction(0) =  2*MATERIAL_PARAM * (lam - 1.0/lam)              +  ACTIVE_TENSION*lam/denom;
        traction(1) = -2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam)  -  ACTIVE_TENSION*mu/denom;
    }
    else if (fabs(location(1))  <= DBL_EPSILON) //Bottom edge
    {
        traction(0) =  2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam);
        traction(1) = -2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else if (fabs(location(1) - 1.0) <= DBL_EPSILON)//Top edge
    {
        traction(0) = -2*MATERIAL_PARAM * location(1)*ALPHA/(lam*lam);
        traction(1) =  2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else
    {
        NEVER_REACHED;
    }

    return traction;
}

/**
 *  Test against an exact solution.
 *
 *  Choose x=X+0.5*alpha*X*X, y=Y/(1+alpha*X), p=2c, then F has determinant 1. Take
 *  the active tension to be constant.
 *
 *  Then S can be computed, and the boundary tractions and body needed to give this
 *  deformation as a solution can be computed, they are implemented above.
 *
 *  As active tension is constant, we really just testing cardiac mechanics
 *  here, not electrics solution, but electrics is tested elsewhere and mechanics does
 *  not affect electrics by default. Class is called TestElectroMechanicsExactSolution
 *  because we use the CardiacElectroMechanicsProblem interface)
 *
 *  See also equations and finite element implementation document.
 *
 */

class TestElectroMechanicsExactSolution : public CxxTest::TestSuite
{
public:
    void TestIncompressibleSolveWithExactSolution()
    {
        MechanicsEventHandler::Reset();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(0.0);
        HeartConfig::Instance()->SetSimulationDuration(1.0);


        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.02, 1.0, 1.0);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.05, 1.0, 1.0);

        MooneyRivlinMaterialLaw<2> law(MATERIAL_PARAM);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh,0,0);

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mechanics_mesh.GetBoundaryElementIteratorBegin();
            iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            // get all boundary elems except those on X=0
            if (fabs((*iter)->CalculateCentroid()[0])>1e-6)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
            }
        }

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetBodyForce(MyBodyForce);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, MyTraction);
        problem_defn.SetContractionModel(CONSTANT,1.0/*ODE timestep, unused*/);
        problem_defn.SetSolverType(EXPLICIT);
        problem_defn.SetMechanicsSolveTimestep(1.0);


        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestElectroMechanicsExactSolution");

        problem.Initialise();


        /////////////////////////////////////////////////////////////////
        //
        // Get access to the contraction models, to set the constant
        // to something other than 1.0. Difficult to do at the moment
        //
        //////////////////////////////////////////////////////////////////
        ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>* p_solver
            = dynamic_cast<ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>*>(problem.mpCardiacMechSolver);

        for (std::map<unsigned,DataAtQuadraturePoint>::iterator iter = p_solver->rGetQuadPointToDataAtQuadPointMap().begin();
             iter != p_solver->rGetQuadPointToDataAtQuadPointMap().end();
             iter++)
        {
            ConstantActiveTension* p_contraction_model = dynamic_cast<ConstantActiveTension*>(iter->second.ContractionModel);
            p_contraction_model->SetActiveTensionValue(ACTIVE_TENSION);
        }


        //////////////////////////////////////////////////////////////////
        //
        //  Solve and test
        //
        /////////////////////////////////////////////////////////////////
        problem.Solve();

        std::vector<c_vector<double,2> >& r_solution = problem.rGetDeformedPosition();

        for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];

            double exact_x = X + 0.5*ALPHA*X*X;
            double exact_y = Y/(1+ALPHA*X);

            TS_ASSERT_DELTA(r_solution[i](0), exact_x, 1e-3);
            TS_ASSERT_DELTA(r_solution[i](1), exact_y, 1e-3);
        }

        std::vector<double>& r_pressures = p_solver->rGetPressures();
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA( r_pressures[i]/(2*MATERIAL_PARAM), 1.0, 5e-3);
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
};

#endif //_TESTELECTROMECHANICSEXACTSOLUTION_HPP_
