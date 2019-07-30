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


#ifndef TESTFIBREWRITER_HPP_
#define TESTFIBREWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>

#include "FibreReader.hpp"
#include "FibreWriter.hpp"
#include "PetscTools.hpp"
#include "UblasIncludes.hpp"
#include "FileComparison.hpp"
#include "FibreConverter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestFibreWriter : public CxxTest::TestSuite
{
public:
    void TestAxiWriterAscii()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);

        // Write ascii file
        FibreWriter<3> fibre_writer("TestFibreWriter", "SimpleAxisymmetric2", true);
        fibre_writer.WriteAllAxi(fibre_vector);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        FileComparison comparer(results_dir + "/SimpleAxisymmetric2.axi","heart/test/data/fibre_tests/SimpleAxisymmetric2.axi");
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestAxiWriterBinary()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);

        //Write binary file
        FibreWriter<3> fibre_writer("TestFibreWriter", "SimpleAxisymmetric2Bin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllAxi(fibre_vector);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        FileComparison comparer(results_dir + "/SimpleAxisymmetric2Bin.axi","heart/test/data/fibre_tests/SimpleAxisymmetric2Bin.axi");
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestOrthoWriterAscii()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);

        //Write ascii file
        FibreWriter<3> fibre_writer("TestFibreWriter", "Orthotropic3D", false);
        fibre_writer.WriteAllOrtho(fibres, second, third);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";

        FileComparison comparer(results_dir + "/Orthotropic3D.ortho","heart/test/data/fibre_tests/Orthotropic3D.ortho");
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestOrthoWriterBinary()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);

        //Write binary file
        FibreWriter<3> fibre_writer("TestFibreWriter", "Orthotropic3DBin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllOrtho(fibres, second, third);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        FileComparison comparer(results_dir + "/Orthotropic3DBin.ortho","heart/test/data/fibre_tests/Orthotropic3DBin.ortho");
        TS_ASSERT(comparer.CompareFiles());
    }


    //Following are convenience methods
private:
    void ConvertToBinaryOrtho(std::string fullPath, std::string baseName, bool copyOutput)
    {
        FileFinder file_finder(fullPath+"/"+baseName+".ortho", RelativeTo::ChasteSourceRoot);

        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
        TS_ASSERT_EQUALS(fibres.size(), second.size());
        TS_ASSERT_EQUALS(fibres.size(), third.size());

        std::string dir = "FibreConverter";
        FibreConverter ortho_converter;
        ortho_converter.Convert(file_finder, dir);

        //Read it back
        FileFinder file_finder_bin(dir+"/"+baseName+"_bin.ortho", RelativeTo::ChasteTestOutput);
        FibreReader<3> fibre_reader_bin(file_finder_bin, ORTHO);
        std::vector< c_vector<double, 3> > fibres_bin;
        std::vector< c_vector<double, 3> > second_bin;
        std::vector< c_vector<double, 3> > third_bin;
        fibre_reader_bin.GetAllOrtho(fibres_bin, second_bin, third_bin);

        TS_ASSERT_EQUALS(fibres.size(), fibres_bin.size());
        TS_ASSERT_EQUALS(second.size(), second_bin.size());
        TS_ASSERT_EQUALS(third.size(), third_bin.size());

        for (unsigned i = 0; i< fibres.size(); i++)
        {
            for (unsigned j = 0; j< 3u ; j++)
            {
                TS_ASSERT_DELTA(fibres[i][j], fibres_bin[i][j], 1e-16);
                TS_ASSERT_DELTA(second[i][j], second_bin[i][j], 1e-16);
                TS_ASSERT_DELTA(third[i][j], third_bin[i][j], 1e-16);
            }
        }
        if (copyOutput)
        {
            // This code is not called by default and only exists to copy test output into the trunk.
            FileFinder file_finder_test_folder(fullPath, RelativeTo::ChasteSourceRoot);
            file_finder_bin.CopyTo(file_finder_test_folder);
        }
    }

    void ConvertToBinaryAxi(std::string fullPath, std::string baseName, bool copyOutput)
    {
        FileFinder file_finder(fullPath+"/"+baseName+".axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibres;
        fibre_reader.GetAllAxi(fibres);

        std::string dir = "FibreConverter";
        FibreConverter axi_converter;
        axi_converter.Convert(file_finder, dir);

        //Read it back
        FileFinder file_finder_bin(dir+"/"+baseName+"_bin.axi", RelativeTo::ChasteTestOutput);
        FibreReader<3> fibre_reader_bin(file_finder_bin, AXISYM);
        std::vector< c_vector<double, 3> > fibres_bin;
        fibre_reader_bin.GetAllAxi(fibres_bin);

        TS_ASSERT_EQUALS(fibres.size(), fibres_bin.size());

        for (unsigned i = 0; i< fibres.size(); i++)
        {
            for (unsigned j = 0; j< 3u ; j++)
            {
                TS_ASSERT_DELTA(fibres[i][j], fibres_bin[i][j], 1e-16);
            }
        }
        if (copyOutput)
        {
            // This code is not called by default and only exists to copy test output into the trunk.
            FileFinder file_finder_test_folder(fullPath, RelativeTo::ChasteSourceRoot);
            file_finder_bin.CopyTo(file_finder_test_folder);
        }
    }

public:
    void TestConvertFiles()
    {
          ConvertToBinaryOrtho("heart/test/data/box_shaped_heart/", "box_heart", false);
          ConvertToBinaryAxi("heart/test/data/box_shaped_heart/", "box_heart", false);
    }

    void doNotTestReallyConvertFiles()
    {
          ConvertToBinaryOrtho("heart/test/data/fibre_tests/", "downsampled", true);
          ConvertToBinaryOrtho("heart/test/data/point50_heart_mesh/", "point50", true);
          ConvertToBinaryAxi("apps/texttest/weekly/Propagation3d/", "OxfordRabbitHeart_482um", true);
          //This one is dodgy... ConvertToBinaryAxi("notforrelease/test/data/simplified_very_low_res/", "heart_chaste2_renum_e_triangles", true);
    }
};

#endif /*TESTFIBREREADER_HPP_*/
