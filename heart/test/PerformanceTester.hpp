/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
    static const double mMeshWidth=0.2;//cm
public:
    PerformanceTester()
    : OdeTimeStep(0.0025),
      PdeTimeStep(0.0025),
      MeshNum(2u),
      KspAtol(5e-5),
      //RelativeConvergenceCriterion(1e-4),
      SimTime(8.0),
      PrintingTimeStep(0.04)
    {
    }

    void Run()
    {
        HeartConfig::Instance()->SetOdeTimeStep(this->OdeTimeStep);
        HeartConfig::Instance()->SetSimulationDuration(SimTime);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(KspAtol);

        // Create the meshes on which the test will be based
        const std::string mesh_dir = "PerformanceMesh/";
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
        HeartConfig::Instance()->SetOutputDirectory ("Performance");
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
    unsigned mNumNodes;
    unsigned mNumElements;

};
#endif /*PERFORMANCETESTER_HPP_*/
