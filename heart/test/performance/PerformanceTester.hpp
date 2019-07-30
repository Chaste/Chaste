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


#ifndef PERFORMANCETESTER_HPP_
#define PERFORMANCETESTER_HPP_

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class PerformanceTester
{
protected:
    static const double mMeshWidth;//=0.2;//cm
public:
    PerformanceTester(const std::string& rTestName)
    : OdeTimeStep(0.0025),
      PdeTimeStep(0.0025),
      MeshNum(2u),
      KspAtol(5e-5),
      //RelativeConvergenceCriterion(1e-4),
      SimTime(8.0),
      PrintingTimeStep(0.04),
      mTestName(rTestName)
    {
    }

    void Run()
    {
        HeartConfig::Instance()->SetOdeTimeStep(this->OdeTimeStep);
        HeartConfig::Instance()->SetSimulationDuration(SimTime);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(KspAtol);

        // Create the meshes on which the test will be based
        const std::string mesh_dir = mTestName+"Mesh/";
        OutputFileHandler output_file_handler(mesh_dir);

        unsigned prev_mesh_num=9999;
        std::string mesh_filename="temp_mesh";

        CuboidMeshConstructor<DIM> constructor;

        if (MeshNum!=prev_mesh_num)
        {
            //Here we temporarily construct a mesh and write files out.
            //This is so that performance metrics include the time to read in the mesh from disk
            TetrahedralMesh<DIM,DIM> mesh;
            constructor.Construct(mesh, MeshNum, mMeshWidth);
            TrianglesMeshWriter<DIM,DIM> writer(mesh_dir, mesh_filename);
            writer.WriteFilesUsingMesh(mesh);
            mNumElements = constructor.GetNumElements();
            mNumNodes = constructor.GetNumNodes();
            prev_mesh_num = MeshNum;
        }

        unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2);
        GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(num_ele_across, mMeshWidth);
        CARDIAC_PROBLEM cardiac_problem(&cell_factory);

        HeartConfig::Instance()->SetMeshFileName(output_file_handler.GetOutputDirectoryFullPath()+mesh_filename);
        HeartConfig::Instance()->SetOutputDirectory (mTestName);
        HeartConfig::Instance()->SetOutputFilenamePrefix ("Results");

//        cardiac_problem.SetLinearSolverRelativeTolerance(KspRtol);


        assert(fabs(0.04/PdeTimeStep - round(0.04/PdeTimeStep)) <1e-15 );

        HeartConfig::Instance()->SetPdeTimeStep(PdeTimeStep);
        HeartConfig::Instance()->SetPrintingTimeStep(PrintingTimeStep);

        cardiac_problem.Initialise();

        //// use this to get some info printed out
        //cardiac_problem.SetWriteInfo();

        try
        {
            cardiac_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout<<"Warning - this run threw an exception.  Check performance results\n";
            throw(e);
        }

        DisplayRun();
    }

    static void DisplayHeadings()
    {
        const unsigned NUM_HEADINGS=7;
        const char* heading[NUM_HEADINGS]={"Dimen", "Elts", "Nodes", "PdeStp", "OdeStp", "PriStp", "SimTim"};
        for (unsigned i=0; i<NUM_HEADINGS; i++)
        {
            printf("%6s\t", heading[i]);
        }
    }

    void DisplayRun()
    {
        printf("%6u\t%6u\t%6u\t%2.1e\t%2.1e\t%2.1e\t%2.1e\t",
               DIM, mNumElements, mNumNodes, PdeTimeStep, OdeTimeStep, PrintingTimeStep, SimTime);
    }


public:
    double OdeTimeStep;
    double PdeTimeStep;
    unsigned MeshNum;
    double KspAtol;
    //double RelativeConvergenceCriterion;
    double SimTime;
    double PrintingTimeStep;

private:
    std::string mTestName;
    unsigned mNumNodes;
    unsigned mNumElements;
};

///Set static const on instantiation
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
const double PerformanceTester<CELL, CARDIAC_PROBLEM, DIM>::mMeshWidth=0.2;//cm

#endif /*PERFORMANCETESTER_HPP_*/
