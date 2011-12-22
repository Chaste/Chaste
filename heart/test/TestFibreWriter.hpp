/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef TESTFIBREWRITER_HPP_
#define TESTFIBREWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>

#include "FibreReader.hpp"
#include "FibreWriter.hpp"
#include "PetscTools.hpp"
#include "UblasIncludes.hpp"


class TestFibreWriter : public CxxTest::TestSuite
{
public:
    void TestAxiWriterAscii()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);

        //Write ascii file
        FibreWriter<3> fibre_writer("TestFibreWriter", "SimpleAxisymmetric2", true);
        fibre_writer.WriteAllAxi(fibre_vector);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_dir + "/SimpleAxisymmetric2.axi heart/test/data/fibre_tests/SimpleAxisymmetric2.axi").c_str()), 0);

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
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/SimpleAxisymmetric2Bin.axi heart/test/data/fibre_tests/SimpleAxisymmetric2Bin.axi").c_str()), 0);
    }

    void TestOrthoWriterAscii() throw (Exception)
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
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_dir + "/Orthotropic3D.ortho heart/test/data/fibre_tests/Orthotropic3D.ortho").c_str()), 0);
    }

    void TestOrthoWriterBinary() throw (Exception)
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
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/Orthotropic3DBin.ortho heart/test/data/fibre_tests/Orthotropic3DBin.ortho").c_str()), 0);
    }


    //Following are convenience methods
private:
    void ConvertToBinaryOrtho(std::string fullPath, std::string baseName, bool output) throw (Exception)
    {
        FileFinder file_finder(fullPath+"/"+baseName+".ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
        TS_ASSERT_EQUALS(fibres.size(), second.size());
        TS_ASSERT_EQUALS(fibres.size(), third.size());
        std::string dir =   "FibreConverter";
        //Write binary file
        FibreWriter<3> fibre_writer(dir, baseName+"_bin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllOrtho(fibres, second, third);

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
        if (output)
        {
            std::cout<<"cp "<<
                OutputFileHandler::GetChasteTestOutputDirectory()
                << dir <<"/"<<baseName<<"_bin.ortho " <<
                fullPath << "/" << baseName << ".ortho\n";
            system(("ls -lh " + OutputFileHandler::GetChasteTestOutputDirectory()
                + dir + "/" + baseName + "_bin.ortho " +
                fullPath + "/" + baseName + ".ortho\n").c_str());
            }
        }
    void ConvertToBinaryAxi(std::string fullPath, std::string baseName, bool output) throw (Exception)
    {
        FileFinder file_finder(fullPath+"/"+baseName+".axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibres;
        fibre_reader.GetAllAxi(fibres);
        std::string dir =   "FibreConverter";
        //Write binary file
        FibreWriter<3> fibre_writer(dir, baseName+"_bin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllAxi(fibres);

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
        if (output)
        {
            std::cout<<"cp "<<
                OutputFileHandler::GetChasteTestOutputDirectory()
                << dir <<"/"<<baseName<<"_bin.axi " <<
                fullPath << "/" << baseName << ".axi\n";
            system(("ls -lh " + OutputFileHandler::GetChasteTestOutputDirectory()
                + dir + "/" + baseName + "_bin.axi " +
                fullPath + "/" + baseName + ".axi\n").c_str());
        }
    }
public:
    void TestConvertFiles() throw (Exception)
    {
          ConvertToBinaryOrtho("heart/test/data/box_shaped_heart/", "box_heart", false);
          ConvertToBinaryAxi("heart/test/data/box_shaped_heart/", "box_heart", false);
    }
    void doNotTestReallyConvertFiles() throw (Exception)
    {
          ConvertToBinaryOrtho("heart/test/data/fibre_tests/", "downsampled", true);
          ConvertToBinaryOrtho("heart/test/data/point50_heart_mesh/", "point50", true);
          ConvertToBinaryAxi("apps/texttest/weekly/Propagation3d/", "heart_chaste2_renum_i_triangles", true);
          //This one is dodgy... ConvertToBinaryAxi("notforrelease/test/data/simplified_very_low_res/", "heart_chaste2_renum_e_triangles", true);
    }

};


#endif /*TESTFIBREREADER_HPP_*/
