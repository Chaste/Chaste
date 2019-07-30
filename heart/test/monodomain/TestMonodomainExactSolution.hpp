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

#ifndef TESTMONODOMAINEXACTSOLUTION_HPP_
#define TESTMONODOMAINEXACTSOLUTION_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"

//////////////////////////////////////////////////////////////////////////////////
//
//  Test the monodomain solver against exact solutions by using fake cell model
//  that returns an ionic current that is a chosen function of position, and
//  ditto the stimulus and initial condition
//
//  The equation solved is:
//
//  Am Cm dV/dt   +   Am Iion  =   div (sigma gradV) - Istim,   with Neumann bcs
//
//  In 2d, say, on the unit square, let
//
//   ---  V(t,x) = exp(-t)*cos(pi*x)cos(pi*y)  ---
//
//  (which satisfies the bcs). This is the solution to the problem if we choose
//  initial conditions V0 = cos(pi*x)cos(pi*y), and
//
//    Iion  = alpha * exp(-t)*cos(pi*x)cos(pi*y)
//    Istim =  beta * exp(-t)*cos(pi*x)cos(pi*y)
//
//  as long as the parameters are chosen to satisfy:
//
//    -Am Cm   +  Am alpha  =  -(sigma00+sigma11) pi^2  - beta
//
//  And similarly in 1d and 3d
//////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
class ToyCellModel : public AbstractCardiacCell
{
private:
    c_vector<double,DIM> mX;

public :
    ToyCellModel(boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus,
                 c_vector<double,DIM> x)
        : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                              1, 0, pIntracellularStimulus)
    {
        mX = x;
        mStateVariables.resize(1);

        ///////////////////////
        // initial condition
        ///////////////////////
        if (DIM==1)
        {
            mStateVariables[0] = cos(M_PI*mX[0]);
        }
        else if (DIM==2)
        {
            mStateVariables[0] = cos(M_PI*mX[0])*cos(M_PI*mX[1]);
        }
        else
        {
            mStateVariables[0] = cos(M_PI*mX[0])*cos(M_PI*mX[1])*cos(M_PI*mX[2]);
        }
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
    }

    ///////////////////////
    // ionic current
    ///////////////////////
    double GetIIonic(const std::vector<double>* pStateVariables=NULL)
    {
        // We use the next timestep to be consistent with how GetStimulus(time) is called.

        double t = PdeSimulationTime::GetNextTime();
        if (DIM==1)
        {
            return (-0.5*M_PI*M_PI + 1.5)*exp(-t)*cos(M_PI*mX[0]);
        }
        else if (DIM==2)
        {
            return (-0.5*M_PI*M_PI + 1.5)*exp(-t)*cos(M_PI*mX[0])*cos(M_PI*mX[1]);
        }
        else
        {
            return (-0.5*M_PI*M_PI + 1.5)*exp(-t)*cos(M_PI*mX[0])*cos(M_PI*mX[1])*cos(M_PI*mX[2]);
        }
    }

    void ComputeExceptVoltage(double tStart, double tEnd)
    {
    }
};


template<unsigned DIM>
class PositionDependentStimulus : public AbstractStimulusFunction
{
private:
    c_vector<double,DIM> mX;

public:
    PositionDependentStimulus(c_vector<double,DIM> x)
    {
        mX = x;
    }

    ///////////////////////
    // stimulus current
    ///////////////////////
    double GetStimulus(double time)
    {
        if (DIM==1)
        {
            return  -M_PI*M_PI *exp(-time)*cos(M_PI*mX[0]);
        }
        else if (DIM==2)
        {

            return  -M_PI*M_PI *exp(-time)*cos(M_PI*mX[0])*cos(M_PI*mX[1]);
        }
        else
        {
            return  -M_PI*M_PI *exp(-time)*cos(M_PI*mX[0])*cos(M_PI*mX[1])*cos(M_PI*mX[2]);
        }
    }
};

template<unsigned DIM>
class MyCellFactory : public AbstractCardiacCellFactory<DIM>
{
public:
    MyCellFactory()
        : AbstractCardiacCellFactory<DIM>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        c_vector<double,DIM> x = pNode->rGetLocation();

        boost::shared_ptr<PositionDependentStimulus<DIM> >
            p_stimulus(new PositionDependentStimulus<DIM>(x));

        ToyCellModel<DIM>* p_cell_model = new ToyCellModel<DIM>(p_stimulus, x);
        return p_cell_model;
    }
};

class TestMonodomainExactSolution : public CxxTest::TestSuite
{
public:
    void TestMonodomainExactSolution1d()
    {
        TetrahedralMesh<1,1> mesh;
        double h=0.01; //cm

        mesh.ConstructRegularSlabMesh(h, 1.0);

        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainExactSolution1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetSimulationDuration(1); //ms


        // Note pde_timestep=0.01 is hardcoded above in GetIonic, and printing
        // timestep hardcoding below in test!
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2);
        HeartConfig::Instance()->SetCapacitance(1.5);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(2));

        MyCellFactory<1> cell_factory;

        MonodomainProblem<1> problem( &cell_factory );
        problem.SetMesh(&mesh);
        problem.Initialise();

        problem.SetWriteInfo();

        problem.Solve();


        Hdf5DataReader reader("TestMonodomainExactSolution1d","results");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        DistributedVectorFactory factory(mesh.GetNumNodes());
        Vec voltage = factory.CreateVec();

        for (unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(voltage, "V", timestep);
            ReplicatableVector voltage_repl(voltage);

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];

                double V = voltage_repl[i];
                double t = 0.01*timestep;
                double exact_solution = exp(-t)*cos(M_PI*x);

                // for testing suitable tols:
                if (fabs(V-exact_solution)>0.002*exp(-t))
                {
                    std::cout << "t, i = " << t << " " << i << "\n";
                    std::cout << "V, exact_solution = " << V << " " << exact_solution << "\n";
                    assert(0);
                }

                // 0.1% error tolerance. Note that the tol is scaled by exp(-t)=max(V)
                TS_ASSERT_DELTA(V, exact_solution, 0.001*exp(-t));
            }
        }
    }

    void TestMonodomainExactSolution2d()
    {
        TetrahedralMesh<2,2> mesh;
        double h=0.02; //cm

        mesh.ConstructRegularSlabMesh(h, 1.0, 1.0);

        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainExactSolution2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetSimulationDuration(1); //ms


        // Note pde_timestep=0.01 is hardcoded above in GetIonic, and printing
        // timestep hardcoding below in test!
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2);
        HeartConfig::Instance()->SetCapacitance(1.5);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.2,0.8));

        MyCellFactory<2> cell_factory;

        MonodomainProblem<2> problem( &cell_factory );
        problem.SetMesh(&mesh);
        problem.Initialise();

        problem.SetWriteInfo();

        problem.Solve();

        Hdf5DataReader reader("TestMonodomainExactSolution2d","results");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        DistributedVectorFactory factory(mesh.GetNumNodes());
        Vec voltage = factory.CreateVec();

        for (unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(voltage, "V", timestep);
            ReplicatableVector voltage_repl(voltage);

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x=mesh.GetNode(i)->rGetLocation()[0];
                double y=mesh.GetNode(i)->rGetLocation()[1];

                double V = voltage_repl[i];
                double t = 0.01*timestep;
                double exact_solution = exp(-t)*cos(M_PI*x)*cos(M_PI*y);

                // for testing suitable tols:
                if (fabs(V-exact_solution)>0.002*exp(-t))
                {
                    std::cout << "t, i = " << t << " " << i << "\n";
                    std::cout << "V, exact_solution = " << V << " " << exact_solution << "\n";
                    assert(0);
                }

                //// With UNIT parameters (Am=Cm=1, sigma=(1,0; 0,1)):
                // 10 elements each dir: passes tol = 4%, fails tol = 3%
                // 20 elements each dir: passes tol = 1%, fails tol = 0.5%
                // 50 elements each dir: passes tol = 0.2%, fails tol = 0.1%

                // 0.2% error tolerance. Note that the tol is scaled by exp(-t)=max(V)
                TS_ASSERT_DELTA(V, exact_solution, 0.002*exp(-t));

            }
        }
    }


    void TestMonodomainExactSolution3d()
    {
        TetrahedralMesh<3,3> mesh;
        double h=0.05; //cm

        mesh.ConstructRegularSlabMesh(h, 1.0, 1.0, 1.0);

        HeartConfig::Instance()->SetOutputDirectory("TestMonodomainExactSolution3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetSimulationDuration(1); //ms

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2);
        HeartConfig::Instance()->SetCapacitance(1.5);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.8, 0.7, 0.5));

        MyCellFactory<3> cell_factory;

        MonodomainProblem<3> problem( &cell_factory );
        problem.SetMesh(&mesh);
        problem.Initialise();

        problem.SetWriteInfo();

        problem.Solve();


        Hdf5DataReader reader("TestMonodomainExactSolution3d","results");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        DistributedVectorFactory factory(mesh.GetNumNodes());
        Vec voltage = factory.CreateVec();

        for (unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(voltage, "V", timestep);
            ReplicatableVector voltage_repl(voltage);

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x=mesh.GetNode(i)->rGetLocation()[0];
                double y=mesh.GetNode(i)->rGetLocation()[1];
                double z=mesh.GetNode(i)->rGetLocation()[2];

                double V = voltage_repl[i];
                double t = 0.01*timestep;
                double exact_solution = exp(-t)*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);

                // for testing suitable tols:
                if (fabs(V-exact_solution)>0.03*exp(-t))
                {
                    std::cout << "t, i = " << t << " " << i << "\n";
                    std::cout << "V, exact_solution = " << V << " " << exact_solution << "\n";
                    assert(0);
                }

                //// With UNIT parameters (Am=Cm=1, sigma=(1,1,1)):
                // 10 elements each dir: passes tol = 7%, fails tol = 6%
                // 20 elements each dir: passes tol = 2%, fails tol = 1%
                //  --look reasonable compared to 2d results above, since would expect error
                //    constant to scale with dimension

                // 3% error tolerance. Note that the tol is scaled by exp(-t)=max(V)
                TS_ASSERT_DELTA(V, exact_solution, 0.03*exp(-t));
            }
        }
    }
};

#endif /*TESTMONODOMAINEXACTSOLUTION_HPP_*/
