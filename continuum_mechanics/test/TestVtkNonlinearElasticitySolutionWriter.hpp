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


#ifndef _TESTVTKNONLINEARELASTICITYWRITER_HPP_
#define _TESTVTKNONLINEARELASTICITYWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/assign.hpp>

#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "NonlinearElasticityTools.hpp"
#include "SolidMechanicsProblemDefinition.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "VtkNonlinearElasticitySolutionWriter.hpp"

class TestVtkNonlinearElasticitySolutionWriter : public CxxTest::TestSuite
{
public:
    void TestException()
    {
#ifdef CHASTE_VTK
        QuadraticMesh<3> mesh(1.0, 1.0, 1.0, 1.0);

        // set up a solver object
        MooneyRivlinMaterialLaw<3> law(1.0, 0.0);
        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh, 0, 0.0);
        SolidMechanicsProblemDefinition<3> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        IncompressibleNonlinearElasticitySolver<3> solver(mesh,problem_defn,"");

        VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
        TS_ASSERT_THROWS_THIS(vtk_writer.Write(),"No output directory was given to the mechanics solver");
#endif //CHASTE_VTK
    }

    void TestVtuFile()
    {
#ifdef CHASTE_VTK
        for (unsigned run=0; run<3; run++)
        {
            std::stringstream dir;
            dir << "TestVtkNonlinearElasticityWriter_" << run;

            QuadraticMesh<3> bar_mesh(0.5, 10.0, 1.0, 1.0);

            // set up a solver object
            MooneyRivlinMaterialLaw<3> law(1.0, 0.0);
            std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(bar_mesh, 0, 0.0);
            SolidMechanicsProblemDefinition<3> problem_defn(bar_mesh);
            problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
            problem_defn.SetZeroDisplacementNodes(fixed_nodes);
            IncompressibleNonlinearElasticitySolver<3> solver(bar_mesh,problem_defn,dir.str());
            solver.Solve();

            // solution is currently no deformation. Hack into the solution and set it to be something known.
            for (unsigned i=0; i<bar_mesh.GetNumNodes(); i++)
            {
                double X = bar_mesh.GetNode(i)->rGetLocation()[0];
                double Y = bar_mesh.GetNode(i)->rGetLocation()[1];
                //double Z = bar_mesh.GetNode(i)->rGetLocation()[2];

                double u,v,w,p;
                if (X<5.01)
                {
                    u = 1.5*X - X;
                    v = 0.0;
                    w = 0.0;
                    p = 23.0;
                }
                else
                {
                    u = 1.5*5 + 2*(X-5.0) - X;
                    v = Y + 0.1*(X-5.0) - Y;
                    w = 0;
                    p = 46.0;
                }

                solver.rGetCurrentSolution()[4*i] = u;
                solver.rGetCurrentSolution()[4*i+1] = v;
                solver.rGetCurrentSolution()[4*i+2] = w;
                solver.rGetCurrentSolution()[4*i+3] = p;
            }

            StrainType strain_type = (run==0 ? DEFORMATION_GRADIENT_F : (run==1 ? DEFORMATION_TENSOR_C : LAGRANGE_STRAIN_E) );

            VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
            vtk_writer.SetWriteElementWiseStrains(strain_type);
            vtk_writer.Write();

            // .vtu files have been visualised, everything looks good.
            FileFinder vtk_file(dir.str() + "/vtk/solution.vtu", RelativeTo::ChasteTestOutput);
            TS_ASSERT(vtk_file.Exists());


            // we can't really test the vtu file easily, but the strain data is stored in a member variable so we can
            // test it was computed correctly.
            for (unsigned i=0; i<bar_mesh.GetNumElements(); i++)
            {
                double X = bar_mesh.GetElement(i)->CalculateCentroid()[0];
                if (X < 5)
                {
                    double F[3][3] = { {1.5, 0, 0}, {0, 1, 0}, {0, 0, 1} };
                    double C[3][3] = { {2.25, 0, 0}, {0, 1, 0}, {0, 0, 1} };
                    double E[3][3] = { {0.625, 0, 0}, {0, 0, 0}, {0, 0, 0} };

                    for (unsigned M=0; M<3; M++)
                    {
                        for (unsigned N=0; N<3; N++)
                        {
                            double value = (run==0 ? F[M][N] : (run==1 ? C[M][N] : E[M][N]) );
                            TS_ASSERT_DELTA(vtk_writer.mTensorData[i](M,N), value, 1e-12);
                        }
                    }
                }
                else
                {
                    double F[3][3] = { {2.0, 0, 0}, {0.1, 1, 0}, {0, 0, 1} };
                    double C[3][3] = { {4.01, 0.1, 0}, {0.1, 1, 0}, {0, 0, 1} };
                    double E[3][3] = { {1.505, 0.05, 0}, {0.05, 0, 0}, {0, 0, 0} };

                    for (unsigned M=0; M<3; M++)
                    {
                        for (unsigned N=0; N<3; N++)
                        {
                            double value = (run==0 ? F[M][N] : (run==1 ? C[M][N] : E[M][N]) );
                            TS_ASSERT_DELTA(vtk_writer.mTensorData[i](M,N), value, 1e-12);
                        }
                    }
                }
            }
        }
#endif //CHASTE_VTK
    }
};


#endif //_TESTVTKNONLINEARELASTICITYWRITER_HPP_
