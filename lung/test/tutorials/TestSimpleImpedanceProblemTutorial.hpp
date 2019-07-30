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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTSIMPLEIMPEDANCEPROBLEMTUTORIAL_HPP_
#define TESTSIMPLEIMPEDANCEPROBLEMTUTORIAL_HPP_

/* HOW_TO_TAG Lung/Simulation
 * Calculate transfer impedance of an airway tree using a simple impedance model
 */

/*
 * = An example showing how to calculate transfer impedance of an airway tree using a simple impedance model =
 *
 * In this tutorial we demonstrate the use of !SimpleImpedanceProblem to calculate transfer impedance on an
 * airway tree model. We further demonstrate post-processing of the output using !ImpedancePostProcessor to
 * calculate a number of clinically relevant measures.
 *
 * Note that !SimpleImpedanceProblem uses Poiseuille formulas to calculate impedance
 * rather than more accurate acoustic impedance equations. For the more accurate version see !ImpedanceProblem.
 */

/* The usual headers are included */
#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"

/* !SimpleImpedanceProblem does most of the work in calculating impedance. */
#include "SimpleImpedanceProblem.hpp"

/* !ImpedancePostProcessor allows easy calculation of a number of clinically relevant measures. */
#include "ImpedancePostProcessor.hpp"


/* Define the test */
class TestSimpleImpedanceProblemTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestCalculateImpedance()
    {
        EXIT_IF_PARALLEL;

        /* First, we load up a mesh containing the centre lines and radii of the a complete conducting airway tree.
         * The mesh will typically have been developed using a combination of computed tomography (CT) image segmentation
         * and algorithmic airway generation.
         */
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("lung/test/data/TestSubject002");
        mesh.ConstructFromMeshReader(mesh_reader);

        /* Note that the mesh defined above was developed using a CT scan taken at full inspiration. Impedance is more
         * commonly recorded during tidal breathing. Here we use a simple scaling to bring the airway radii down into the
         * tidal breathing range.
         */
        for (TetrahedralMesh<1,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->rGetNodeAttributes()[0] *= 0.7;
        }

        /* Setup a !SimpleImpedanceProblem and tell it that the given mesh is defined in millimetres
         */
        SimpleImpedanceProblem problem(mesh, 0u);
        problem.SetMeshInMilliMetres();


        /* This vector lists the input frequencies at which to calculate impedance. They must be
         * monotonically increasing.
         */
        std::vector<double> test_frequencies;
        test_frequencies.push_back(1.0);
        test_frequencies.push_back(2.0);
        test_frequencies.push_back(3.0);
        test_frequencies.push_back(5.0);
        test_frequencies.push_back(10.0);
        test_frequencies.push_back(20.0);
        test_frequencies.push_back(30.0);
        problem.SetFrequencies(test_frequencies);               //Set & get frequencies for coverage

        /* The simple impedance model defines a linear spring at each terminal of the airway tree.
         * This method allows us to set the elastance of the whole lung (in Pa/m^3). This elastance
         * is then evenly distributed over the terminals.
         */
        problem.SetElastance(5.8*98.0665*1e3);

        /* Calculates the impedance at the given frequencies*/
        problem.Solve();

        /* Get the calculated impedances. The impedance at each frequency is
         * made up of a real component (the resistance) and a complex component
         * (the elastance).
         */
        std::vector<std::complex<double> > impedances = problem.rGetImpedances();

        /* The impedances calculated above could at this stage be written to a file
         * and plotted. Instead, we make use of !ImpedancePostProcessor to calculate
         * a number of common clinical summary statistics from the data.
         */
        ImpedancePostProcessor processor(test_frequencies, impedances);

        std::cout << "\n";
        std::cout << "R5       = " << real(impedances[3])*1e-6 << " kPa.s.L^-1" << std::endl;
        std::cout << "R20      = " << real(impedances[5])*1e-6 << " kPa.s.L^-1" << std::endl;
        std::cout << "R5 - R20 = " << processor.GetR5MinusR20()*1e-6 << " kPa.s.L^-1" << std::endl;
        std::cout << "Rrs      = " << processor.GetRrs()*1e-6 << " kPa.s.L^-1" << std::endl;
        std::cout << "\n";
        std::cout << "X5       = " << imag(impedances[5])*1e-6 << " kPa.s.L^-1" << std::endl;
        std::cout << "Fres     = " << processor.GetResonantFrequency() << " Hz " << std::endl;
        std::cout << "Ax       = " << processor.GetAx()*1e-6 << " kPa.L^-1" << std::endl;
        std::cout << "\n";
    }
};

#endif /*TESTSIMPLEIMPEDANCEPROBLEMTUTORIAL_HPP_*/
