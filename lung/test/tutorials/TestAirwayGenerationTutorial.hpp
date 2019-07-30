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
#ifndef TESTAIRWAYGENERATIONTUTORIAL_HPP_
#define TESTAIRWAYGENERATIONTUTORIAL_HPP_

/* HOW_TO_TAG Lung/Anatomy definition
 * Generate a complete conducting airway model given segmentations of CT airways and lobes.
 */

/*
 * = An example showing how generate a complete conducting airway model given segmentations of CT airways and lobes =
 *
 * In this tutorial we demonstrate using Chaste's airway generation algorithm to create a complete model of
 * the conducting airways (mean generation ~16) from computed tomography segmentations. Note that the execution
 * time for this tutorial on a standard desktop PC is several minutes.
 */

/*
 * Note that the airway generation code is dependent of having VTK installed.  However, we cannot put a guard around the
 * whole file since that gives compiler errors if VTK is not installed.  Instead we guard the internals of each test, and
 * any includes that will be missing if VTK is not present.
 */
#ifdef CHASTE_VTK

/*
 * We include some VTK classes to allow STL files to be read
 */
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkSTLReader.h"

#endif // CHASTE_VTK

/* The usual headers are included */
#include <cxxtest/TestSuite.h>

/* {{{MultiLobeAirwayGenerator}}} is the class that does most of the work in generating a complete airway tree */
#include "MultiLobeAirwayGenerator.hpp"

/* All test suites should include either {{{PetscSetupAndFinalize}}} or {{{FakePetscSetup}}}.  This code does not
 * currently use any parallel functionality so it might include either.
 */
#include "PetscSetupAndFinalize.hpp"

/* Define the test */
class TestAirwayGenerationTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestGenerateAirways()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        EXIT_IF_PARALLEL;

        /* First, we load up a mesh containing the centre lines and radii of the central airways extracted from
         * a CT image. The mesh needs to be of the type {{{SPACE_DIM=3}}} and {{{ELEMENT_DIM=1}}}; that is it defines a mesh that exists
         * in 3D space and is made up of 1D line elements. Each node in the mesh is expected to have two attributes
         * associated with it. The first attribute specifies the radius of the airways that node. Thus the mesh defines a series
         * of cylinders that represent the airways. The second attribute specifies whether the node is a terminal node or not.
         */
        TetrahedralMesh<1,3> airways_mesh;
        VtkMeshReader<1,3> airways_mesh_reader("lung/test/data/TestSubject002MajorAirways.vtu");
        airways_mesh.ConstructFromMeshReader(airways_mesh_reader);

        /* Note that the central airways mesh used here is defined in VTK unstructured grid format,
         * for this format we have to manually copy over the node attribute information. This step would
         * be unnecessary if using a mesh in !Triangles/Tetgen format.
         */
        std::vector<double> node_radii;
        airways_mesh_reader.GetPointData("radius", node_radii);
        std::vector<double> terminal_marker;
        airways_mesh_reader.GetPointData("start_id", terminal_marker);
        for (TetrahedralMesh<1,3>::NodeIterator iter = airways_mesh.GetNodeIteratorBegin();
             iter != airways_mesh.GetNodeIteratorEnd();
             ++iter)
        {
            iter->AddNodeAttribute(node_radii[iter->GetIndex()]);
            iter->AddNodeAttribute(fmod(terminal_marker[iter->GetIndex()],2));
        }

        /* We now define a {{{MultiLobeAirwayGenerator}}} to allow us to generate the distal airways to form a complete
         * conducting airway tree. {{{MultiLobeAirwayGenerator}}} provides an easy to use interface to {{{AirwayGenerator}}}
         * and facilitates the generation of airways into a complete lung, rather than the user having to do
         * each lobe separately.
         */
        MultiLobeAirwayGenerator generator(airways_mesh);

        /* We need to set a number of parameters to ensure the resulting airway tree is consistent with known
         * human morphometric data. The values given here can be considered standard for human lungs.
         * The most important of these is the {{{NumberOfPointsPerLung}}}, which (approximately)
         * specifies the number of terminals in the tree. The next is the {{{BranchingFraction}}}, which controls how long
         * the generated airways will be. The diameter ratio is used to control the rate at which airway diameters
         * decrease between airway orders.
         */
        generator.SetNumberOfPointsPerLung(15000);
        generator.SetBranchingFraction(0.4);
        generator.SetDiameterRatio(1.15);

        /* These parameters are less important for producing a consistent airway tree, but are useful for
         * debugging etc.
         */
        generator.SetMinimumBranchLength(0.00001);
        generator.SetPointLimit(1);
        generator.SetAngleLimit(180.0);

        /* We now add lobar surface definitions for the five human lung lobes. Less 'lobes' can be
         * added if full lobar segmentation data is not available. Lobes must be tagged as 'left'
         * or 'right' to enable the correct number of acini to be created. Lobes are represented
         * by triangle surface definitions defined in STL files.
         */
        vtkSmartPointer<vtkSTLReader> lll_reader = vtkSmartPointer<vtkSTLReader>::New();
        lll_reader->SetFileName("lung/test/data/lll.stl");
        lll_reader->Update();
        generator.AddLobe(lll_reader->GetOutput(), LEFT);

        vtkSmartPointer<vtkSTLReader> lul_reader = vtkSmartPointer<vtkSTLReader>::New();
        lul_reader->SetFileName("lung/test/data/lul.stl");
        lul_reader->Update();
        generator.AddLobe(lul_reader->GetOutput(), LEFT);

        vtkSmartPointer<vtkSTLReader> rll_reader = vtkSmartPointer<vtkSTLReader>::New();
        rll_reader->SetFileName("lung/test/data/rll.stl");
        rll_reader->Update();
        generator.AddLobe(rll_reader->GetOutput(), RIGHT);

        vtkSmartPointer<vtkSTLReader> rml_reader = vtkSmartPointer<vtkSTLReader>::New();
        rml_reader->SetFileName("lung/test/data/rml.stl");
        rml_reader->Update();
        generator.AddLobe(rml_reader->GetOutput(), RIGHT);

        vtkSmartPointer<vtkSTLReader> rul_reader = vtkSmartPointer<vtkSTLReader>::New();
        rul_reader->SetFileName("lung/test/data/rul.stl");
        rul_reader->Update();
        generator.AddLobe(rul_reader->GetOutput(), RIGHT);

        /* We now perform two preprocessing steps prior to generation. {{{AssignGrowthApices}}} determine
         * which lobe each of the terminal ends of the central airways segmentation are in.
         */
        generator.AssignGrowthApices();

        /* Distribute points creates the target acinar points within the lung volume. The number
         * created is as specified previously using {{{SetNumberOfPointsPerLung}}}.
         */
        generator.DistributePoints();

        /* We now generate the distal airways. The output is automatically written as a mesh in
         * both tetgen format and VTK unstructured grid format to
         * {{{$CHASTE_TEST_OUTPUT/TestAirwayGenerationTutorial/}}}
         *
         * The resulting geometry can most easily be viewed in Paraview by loading the unstructured grid
         * file. Application of a 'Extract Surface' filter followed by a 'Tube' filter allows the centreline
         * and radius information to be view as a series of tubes.
         */
        generator.Generate("TestAirwayGenerationTutorial", "example_complete_conducting_airway");

#endif // VTK >= 5.6
    }
};

#endif /*TESTAIRWAYGENERATIONTUTORIAL_HPP_*/
