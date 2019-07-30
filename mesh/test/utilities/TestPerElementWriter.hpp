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


#ifndef _TESTPERELEMENTWRITER_HPP_
#define _TESTPERELEMENTWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "AbstractPerElementWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "NumericFileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"


class CentroidWriter : public AbstractPerElementWriter<3,3,3>
{
public:
    CentroidWriter(AbstractTetrahedralMesh<3, 3>* pMesh)
     : AbstractPerElementWriter<3,3,3>(pMesh)
    {
    }

private:
    void Visit(Element<3,3>* pElement, unsigned localElementIndex, c_vector<double, 3>& rData)
    {
        rData = pElement->CalculateCentroid();
        //Check for a data invariant
        TS_ASSERT_EQUALS(mpMesh->mElements[localElementIndex], pElement);
    }
};

class CentroidWithIndexWriter : public AbstractPerElementWriter<3,3,4>
{
public:
    CentroidWithIndexWriter(AbstractTetrahedralMesh<3, 3>* pMesh)
     : AbstractPerElementWriter<3,3,4>(pMesh)
    {
    }

private:
    void Visit(Element<3,3>* pElement, unsigned localElementIndex, c_vector<double, 4>& rData)
    {
        c_vector<double,3> centroid = pElement->CalculateCentroid();
        rData[0] = pElement->GetIndex();
        for (unsigned i=0; i<3; i++)
        {
            rData[i+1] = centroid[i];
        }
    }

    void WriteElementOnMaster(const c_vector<double, 4>& rData)
    {
        //Put square bracket on the first piece of data
        (*mpMasterFile)<<"["<<rData[0]<<"]\t";
        for (unsigned i=1; i<4; i++)
        {
            (*mpMasterFile)<<rData[i]<<"\t";
        }
        (*mpMasterFile)<<"\n";
    }

    void WriteHeaderOnMaster()
    {
        *(this->mpMasterFile)<<"This one has indices";
    }
};

class TestPerElementWriter : public CxxTest::TestSuite
{

public:

    void TestPerElement()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 3, 4);
        OutputFileHandler handler("TestPerElementWriter");

        CentroidWriter writer(&mesh);
        writer.WriteData(handler, "centroid.dat");
        NumericFileComparison(handler.GetOutputDirectoryFullPath() + "/centroid.dat",
                              "mesh/test/data/TestUtilities/centroid.dat").CompareFiles();

        CentroidWithIndexWriter indexed_writer(&mesh);
        indexed_writer.WriteData(handler, "centroid_indexed.dat");
        NumericFileComparison(handler.GetOutputDirectoryFullPath() + "/centroid_indexed.dat",
                              "mesh/test/data/TestUtilities/centroid_indexed.dat").CompareFiles();
    }
};

#endif //_TESTPERELEMENTWRITER_HPP_
