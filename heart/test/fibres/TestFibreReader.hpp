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


#ifndef TESTFIBREREADER_HPP_
#define TESTFIBREREADER_HPP_

#include <cxxtest/TestSuite.h>

#include "FibreReader.hpp"
#include "HeartFileFinder.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"

// simple helper function
template<unsigned DIM>
double UblasMatrixInfinityNorm(c_matrix<double,DIM,DIM> mat)
{
    double ret = fabs(mat(0,0));
    for (unsigned i=0; i<DIM; i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            if (fabs(mat(i,j)) > ret)
            {
                ret = fabs(mat(i,j));
            }
        }
    }
    return ret;
}

class TestFibreReader : public CxxTest::TestSuite
{
public:
    void TestOrthoReaderSetup()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, ORTHO);

        TS_ASSERT_EQUALS(fibre_reader.IsBinary(), false);

        TS_ASSERT_EQUALS(fibre_reader.GetNumLinesOfData(), 5u);

        c_matrix<double, 2, 2> fibre_matrix;

        fibre_reader.GetFibreSheetAndNormalMatrix(0u, fibre_matrix);
        c_matrix<double, 2, 2> correct_matrix = identity_matrix<double>(2,2);
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetFibreSheetAndNormalMatrix(1u, fibre_matrix);
        correct_matrix(1,1) = -1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // this one isn't orthogonal - the false prevents this being checked
        fibre_reader.GetFibreSheetAndNormalMatrix(2u, fibre_matrix, false);
        correct_matrix(0,1) = 1.0;
        correct_matrix(1,0) = 1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // Non-symmetrical test case, standard rotation matrix
        // [cos(theta)  sin(theta)]fibre_matrix(0,0),
        // [-sin(theta) cos(theta)]
        correct_matrix(0,0) =  1.0/sqrt(2.0);  // fibre0
        correct_matrix(1,0) = -1.0/sqrt(2.0);  // fibre1
        correct_matrix(0,1) =  1.0/sqrt(2.0);  // sheet0
        correct_matrix(1,1) =  1.0/sqrt(2.0);  // sheet1
        fibre_reader.GetFibreSheetAndNormalMatrix(3u, fibre_matrix);
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // next matrix is not orthogonal, here we make sure this is checked
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreSheetAndNormalMatrix(4u, fibre_matrix), "not orthogonal")

        // called out of order
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreSheetAndNormalMatrix(4u, fibre_matrix), "Fibre reads must be monotonically increasing")

        // called too many times
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreSheetAndNormalMatrix(5u, fibre_matrix), "End of file")
   }

    void TestOrthoReaderSkipping()  // Cf above test
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, ORTHO);

        c_matrix<double, 2, 2> fibre_matrix;

        fibre_reader.GetFibreSheetAndNormalMatrix(1u, fibre_matrix);
        c_matrix<double, 2, 2> correct_matrix = identity_matrix<double>(2,2);
        correct_matrix(1,1) = -1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // this one isn't orthogonal - the false prevents this being checked
        fibre_reader.GetFibreSheetAndNormalMatrix(4u, fibre_matrix, false);
        correct_matrix(0,0) = 2.0;
        correct_matrix(1,0) = 0.0;
        correct_matrix(0,1) = 1.0;
        correct_matrix(1,1) = 0.0;

        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);
    }

    void TestAxiReaderSetup()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, AXISYM);

        TS_ASSERT_EQUALS(fibre_reader.GetNumLinesOfData(), 4u);

        c_vector<double, 2> fibre_vector;

        fibre_reader.GetFibreVector(0u, fibre_vector);
        TS_ASSERT_DELTA(fibre_vector(0), 0, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), 1, 1e-9);

        fibre_reader.GetFibreVector(1u, fibre_vector);
        TS_ASSERT_DELTA(fibre_vector(0), 0.6, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), -0.8, 1e-9);

        // next vector is not normalised, the false below prevents this being checked
        fibre_reader.GetFibreVector(2u, fibre_vector, false);
        TS_ASSERT_DELTA(fibre_vector(0), 2, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), 6, 1e-9);

        // next vector is not normalised, here we make sure this is checked
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreVector(3u, fibre_vector), "not normalised")

        // called out of order
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreVector(3u, fibre_vector), "Fibre reads must be monotonically increasing")

        // called too many times
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetFibreVector(4u, fibre_vector), "End of file")
    }

    void TestAxiReaderSkipping()  // Cf above test
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, AXISYM);

        c_vector<double, 2> fibre_vector;
        fibre_reader.GetFibreVector(1u, fibre_vector);
        TS_ASSERT_DELTA(fibre_vector(0), 0.6, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), -0.8, 1e-9);

        fibre_reader.GetFibreVector(3u, fibre_vector, false /* don't check normalised */);
        TS_ASSERT_DELTA(fibre_vector(0), 1.0, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), 0.1, 1e-9);
    }


    void TestFibreConvenienceMethodsForVtk()
    {

        {
            FileFinder file("heart/test/data/fibre_tests/SimpleAxisymmetric.axi", RelativeTo::ChasteSourceRoot);
            FibreReader<3> fibre_reader(file, AXISYM);
            std::vector< c_vector<double, 3> > fibres;
            fibre_reader.GetAllAxi(fibres);
            TS_ASSERT_EQUALS(fibres.size(), 6u);
            for (unsigned i=0; i<fibres.size()-1; i++)
            {
                TS_ASSERT_DELTA(fibres[i][0], 1.0, 1e-10);//x
                TS_ASSERT_DELTA(fibres[i][1], 0.0, 1e-10);//y
                TS_ASSERT_DELTA(fibres[i][2], 0.0, 1e-10);//z
            }
            //Last element differs
            TS_ASSERT_DELTA(fibres[5][0], 0.0, 1e-10);//x
            TS_ASSERT_DELTA(fibres[5][1], 0.0, 1e-10);//y
            TS_ASSERT_DELTA(fibres[5][2], 1.0, 1e-10);//z
        }
        {
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
            FibreReader<3> fibre_reader(file, ORTHO);
            std::vector< c_vector<double, 3> > fibres;
            std::vector< c_vector<double, 3> > second;
            std::vector< c_vector<double, 3> > third;
            fibre_reader.GetAllOrtho(fibres, second, third);
            TS_ASSERT_EQUALS(fibres.size(), 6u);
            TS_ASSERT_EQUALS(second.size(), 6u);
            TS_ASSERT_EQUALS(third.size(), 6u);
            for (unsigned i=0; i<fibres.size()-1; i++)
            {
                TS_ASSERT_DELTA(fibres[i][0], 1.0, 1e-10);//x
                TS_ASSERT_DELTA(fibres[i][1], 0.0, 1e-10);//y
                TS_ASSERT_DELTA(fibres[i][2], 0.0, 1e-10);//z
            }
            //Last element differs
            TS_ASSERT_DELTA(fibres[5][0], 0.0, 1e-10);//x
            TS_ASSERT_DELTA(fibres[5][1], 1.0, 1e-10);//y
            TS_ASSERT_DELTA(fibres[5][2], 0.0, 1e-10);//z
            TS_ASSERT_DELTA(second[5][0], 0.0, 1e-10);//x
            TS_ASSERT_DELTA(second[5][1], 0.0, 1e-10);//y
            TS_ASSERT_DELTA(second[5][2], 1.0, 1e-10);//z
            TS_ASSERT_DELTA(third[5][0], 1.0, 1e-10);//x
            TS_ASSERT_DELTA(third[5][1], 0.0, 1e-10);//y
            TS_ASSERT_DELTA(third[5][2], 0.0, 1e-10);//z
        }
    }

    void TestFibretoVtk()
    {
#ifdef CHASTE_VTK
        //See TestConductivityTensors
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "simple_fibres", false);

        {
            FileFinder file("heart/test/data/fibre_tests/SimpleAxisymmetric.axi", RelativeTo::ChasteSourceRoot);
            FibreReader<3> fibre_reader(file, AXISYM);
            std::vector< c_vector<double, 3> > fibres;
            fibre_reader.GetAllAxi(fibres);
            TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
            writer.AddCellData("AxiFibres", fibres);
        }
        {
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
            FibreReader<3> fibre_reader(file, ORTHO);
            std::vector< c_vector<double, 3> > fibres;
            std::vector< c_vector<double, 3> > second;
            std::vector< c_vector<double, 3> > third;
            fibre_reader.GetAllOrtho(fibres, second, third);
            TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
            TS_ASSERT_EQUALS(second.size(), mesh.GetNumElements());
            TS_ASSERT_EQUALS(third.size(), mesh.GetNumElements());
            writer.AddCellData("OrthoFibres", fibres);
            writer.AddCellData("OrthoSecond", second);
            writer.AddCellData("OrthoThird", third);
        }
        writer.WriteFilesUsingMesh(mesh);
        //Check that it has been written
        OutputFileHandler handler("TestVtkMeshWriter", false);
        std::ifstream vtk_file;
        std::string command = handler.GetOutputDirectoryFullPath()+"/simple_fibres.vtu";
        vtk_file.open(command.c_str());
        TS_ASSERT(vtk_file.is_open());
        vtk_file.close();
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestFibreReaderExceptions()
    {
        c_matrix<double, 2, 2> fibre_matrix;

        // file doesn't exist
        FileFinder finder0("heart/test/data/fibre_tests/dgfsdgjdf.ortho", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(finder0,ORTHO), "Failed to open fibre file");

        // line for first element is incomplete
        FileFinder finder1("heart/test/data/fibre_tests/bad_ortho1.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader1(finder1, ORTHO);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader1.GetFibreSheetAndNormalMatrix(0u, fibre_matrix), "A line is incomplete in");

        // line for third element is missing
        FileFinder finder2("heart/test/data/fibre_tests/bad_ortho2.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader2(finder2, ORTHO);
        fibre_reader2.GetFibreSheetAndNormalMatrix(0u, fibre_matrix);
        fibre_reader2.GetFibreSheetAndNormalMatrix(1u, fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader2.GetFibreSheetAndNormalMatrix(2u, fibre_matrix), "End of file");

        // line for second element has too many entries
        FileFinder finder3("heart/test/data/fibre_tests/bad_ortho3.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader3(finder3, ORTHO);
        fibre_reader3.GetFibreSheetAndNormalMatrix(0u, fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader3.GetFibreSheetAndNormalMatrix(1u, fibre_matrix), "Too many entries in a line in");

        // first line doesn't give the number of lines of data
        FileFinder finder4("heart/test/data/fibre_tests/bad_ortho4.ortho", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(finder4,ORTHO), "First (non comment) line of the fibre orientation file should contain the number of lines");

        // Wrong method call, can't read an 'orthotropic vector'
        c_vector<double, 2> fibre_vector;
        FileFinder finder5("heart/test/data/fibre_tests/random_fibres.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader5(finder5, ORTHO);
        TS_ASSERT_THROWS_THIS(fibre_reader5.GetFibreVector(0u, fibre_vector), "Use GetFibreSheetAndNormalMatrix when reading orthotropic fibres");
        std::vector<c_vector<double,2> > v1;
        TS_ASSERT_THROWS_THIS(fibre_reader5.GetAllAxi(v1), "Use GetAllOrtho when reading orthotropic fibres");
        // wrong method call, can't read an 'axisymmetric matrix'
        FileFinder finder6("heart/test/data/fibre_tests/random_fibres.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader6(finder6, AXISYM);
        TS_ASSERT_THROWS_THIS(fibre_reader6.GetFibreSheetAndNormalMatrix(0u, fibre_matrix), "Use GetFibreVector when reading axisymmetric fibres");
        std::vector<c_vector<double,2> > v2;
        std::vector<c_vector<double,2> > v3;
        TS_ASSERT_THROWS_THIS(fibre_reader6.GetAllOrtho(v1, v2, v3), "Use GetAllAxi when reading axisymmetric fibres");

        // Incomplete axi data
        FileFinder finder7("heart/test/data/fibre_tests/bad_axi.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader7(finder7, AXISYM);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader7.GetFibreVector(0u, fibre_vector), "A line is incomplete in");
    }

    void TestAxiBinaryFileReader()
    {
        // Read in a binary fibres file.
        FileFinder file_finder_bin("heart/test/data/fibre_tests/SimpleAxisymmetric2Bin.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader_bin(file_finder_bin, AXISYM);

        TS_ASSERT_EQUALS(fibre_reader_bin.IsBinary(), true);

        std::vector< c_vector<double, 3> > fibre_vector_bin;
        fibre_reader_bin.GetAllAxi(fibre_vector_bin);

        // Read in the equivalent ascii fibres file.
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);

        TS_ASSERT_EQUALS(fibre_vector_bin.size(), fibre_vector.size());
        for (unsigned i=0; i<fibre_vector_bin.size(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(fibre_vector_bin[i][j], fibre_vector[i][j], 1e-9);
            }
        }
    }

    void TestAxiBinaryFileReaderWithSkipping()
    {
        c_vector<double, 3> binary_vector;
        c_vector<double, 3> ascii_vector;
        // Read in a binary fibres file.
        FileFinder file_finder_bin("heart/test/data/fibre_tests/SimpleAxisymmetric2Bin.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader_bin(file_finder_bin, AXISYM);

        // Read in the equivalent ascii fibres file.
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);

        for (unsigned i=0; i<6u; i += 2u)
        {
            fibre_reader.GetFibreVector(i, ascii_vector);
            fibre_reader_bin.GetFibreVector(i, binary_vector);
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(binary_vector[j], ascii_vector[j], 1e-9);
            }
        }
    }

    void TestOrthoBinaryFileReader()
    {
        // Read in a binary fibres file.
        FileFinder file_finder_bin("heart/test/data/fibre_tests/Orthotropic3DBin.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader_bin(file_finder_bin, ORTHO);
        std::vector< c_vector<double, 3> > fibre_vector_bin;
        std::vector< c_vector<double, 3> > second_vector_bin;
        std::vector< c_vector<double, 3> > third_vector_bin;
        fibre_reader_bin.GetAllOrtho(fibre_vector_bin, second_vector_bin, third_vector_bin);

        // Read in the equivalent ascii fibres file.
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibre_vector;
        std::vector< c_vector<double, 3> > second_vector;
        std::vector< c_vector<double, 3> > third_vector;
        fibre_reader.GetAllOrtho(fibre_vector, second_vector, third_vector);

        TS_ASSERT_EQUALS(fibre_vector_bin.size(),  fibre_vector.size());
        TS_ASSERT_EQUALS(second_vector_bin.size(), fibre_vector.size());
        TS_ASSERT_EQUALS(third_vector_bin.size(),  fibre_vector.size());
        TS_ASSERT_EQUALS(second_vector.size(),     fibre_vector.size());
        TS_ASSERT_EQUALS(third_vector.size(),      fibre_vector.size());
        for (unsigned i=0; i<fibre_vector.size(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(fibre_vector_bin[i][j], fibre_vector[i][j], 1e-9);
                TS_ASSERT_DELTA(second_vector_bin[i][j], second_vector[i][j], 1e-9);
                TS_ASSERT_DELTA(third_vector_bin[i][j], third_vector[i][j], 1e-9);
            }
        }
    }

    void TestOrthoBinaryFileReaderWithSkipping()
    {
        // Read in a binary fibres file.
        FileFinder file_finder_bin("heart/test/data/fibre_tests/Orthotropic3DBin.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader_bin(file_finder_bin, ORTHO);
        std::vector< c_vector<double, 3> > fibre_vector_bin;
        std::vector< c_vector<double, 3> > second_vector_bin;
        std::vector< c_vector<double, 3> > third_vector_bin;
        fibre_reader_bin.GetAllOrtho(fibre_vector_bin, second_vector_bin, third_vector_bin);

        // Read in the equivalent ascii fibres file.
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibre_vector;
        std::vector< c_vector<double, 3> > second_vector;
        std::vector< c_vector<double, 3> > third_vector;
        fibre_reader.GetAllOrtho(fibre_vector, second_vector, third_vector);

        TS_ASSERT_EQUALS(fibre_vector_bin.size(),  fibre_vector.size());
        TS_ASSERT_EQUALS(second_vector_bin.size(), fibre_vector.size());
        TS_ASSERT_EQUALS(third_vector_bin.size(),  fibre_vector.size());
        TS_ASSERT_EQUALS(second_vector.size(),     fibre_vector.size());
        TS_ASSERT_EQUALS(third_vector.size(),      fibre_vector.size());
        for (unsigned i=0; i<fibre_vector.size(); i+=2)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(fibre_vector_bin[i][j], fibre_vector[i][j], 1e-9);
                TS_ASSERT_DELTA(second_vector_bin[i][j], second_vector[i][j], 1e-9);
                TS_ASSERT_DELTA(third_vector_bin[i][j], third_vector[i][j], 1e-9);
            }
        }
    }
};

#endif /*TESTFIBREREADER_HPP_*/
