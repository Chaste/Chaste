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

#include "CmguiDeformedSolutionsWriter.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CmguiDeformedSolutionsWriter<DIM>::CmguiDeformedSolutionsWriter(std::string outputDirectory,
                                                                std::string baseName,
                                                                AbstractTetrahedralMesh<DIM,DIM>& rQuadraticMesh,
                                                                CmguiMeshWriteType writeType)
    : CmguiMeshWriter<DIM, DIM>(outputDirectory, baseName),
      mpQuadraticMesh(&rQuadraticMesh),
      mFinalCounter(0)
{
    QuadraticMesh<DIM>* p_quad_mesh = dynamic_cast<QuadraticMesh<DIM>* >(&rQuadraticMesh);

    if (p_quad_mesh == nullptr)
    {
        EXCEPTION("CmguiDeformedSolutionsWriter only supports use of a QuadraticMesh");
    }

    mNumNodesToUse = mpQuadraticMesh->GetNumVertices();

    if (writeType == WRITE_QUADRATIC_MESH)
    {
        mNumNodesToUse = mpQuadraticMesh->GetNumNodes();

        switch (DIM)
        {
            //WriteCmguiScript Commented as CmguiDeformedSolutionsWriter is meant to correspond to
            // output of nonlinear elasticity problems - 2d or 3d only, and there is
            // no explicit instantiation of this class in 1d.
            //case 1:
            //{
            //   this->mElementFileHeader = CmguiElementFileHeader1DQuadratic;
            //    this->mCoordinatesFileHeader = CmguiCoordinatesFileHeader1DQuadratic;
            //    this->mAdditionalFieldHeader = CmguiAdditionalFieldHeader1DQuadratic;
            //    this->mNumNodesPerElement = 3;
            //    this->mReordering.resize(this->mNumNodesPerElement);
            //    unsigned reordering[6] = {0,2,1};
            //    for (unsigned i=0; i<3; i++)
            //    {
            //        this->mReordering[i] = reordering[i];
            //    }
            //    break;
            //};

            case 2:
            {
                this->mElementFileHeader = CmguiElementFileHeader2DQuadratic;
                this->mCoordinatesFileHeader = CmguiCoordinatesFileHeader2DQuadratic;
                this->mAdditionalFieldHeader = CmguiAdditionalFieldHeader2DQuadratic;
                this->mNumNodesPerElement = 6;
                this->mReordering.resize(this->mNumNodesPerElement);

                // Go from Chaste(tetgen ordering) (see example comments in
                // QuadraticBasisFunction::ComputeBasisFunction() to CMGUI ordering
                // ("psi1 increasing, then psi1 increasing")
                unsigned reordering[6] = {0,5,1,4,3,2};
                for (unsigned i=0; i<6; i++)
                {
                    this->mReordering[i] = reordering[i];
                }
                break;
            }
            case 3:
            {
                this->mElementFileHeader = CmguiElementFileHeader3DQuadratic;
                this->mCoordinatesFileHeader = CmguiCoordinatesFileHeader3DQuadratic;
                this->mAdditionalFieldHeader = CmguiAdditionalFieldHeader3DQuadratic;
                this->mNumNodesPerElement = 10;
                this->mReordering.resize(this->mNumNodesPerElement);

                // Go from Chaste(tetgen ordering) (see example comments in
                // QuadraticBasisFunction::ComputeBasisFunction() to CMGUI ordering
                // ("psi1 increasing, then psi2 increasing, then psi3 increasing")
                unsigned reordering[10] = {0,4,1,6,5,2,7,8,9,3};
                for (unsigned i=0; i<10; i++)
                {
                    this->mReordering[i] = reordering[i];
                }
                break;
            }
            default:
            {
                NEVER_REACHED;
            }
        }
    }
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteInitialMesh(std::string fileName)
{
    std::string saved_base_name = this->mBaseName;
    if (fileName == "")
    {
        this->mBaseName = this->mBaseName + "_0";
    }
    else
    {
        this->mBaseName = fileName;
    }
    this->WriteFilesUsingMesh(*mpQuadraticMesh);
    this->mBaseName = saved_base_name;
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteDeformationPositions(std::vector<c_vector<double,DIM> >& rDeformedPositions,
                                                                  unsigned counter)
{
    if (mpQuadraticMesh->GetNumNodes() != rDeformedPositions.size() )
    {
        EXCEPTION("The size of rDeformedPositions does not match the number of nodes in the mesh");
    }

    mFinalCounter = counter;
    std::stringstream node_file_name_stringstream;
    node_file_name_stringstream <<  this->mBaseName << "_" << counter << ".exnode";

    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name_stringstream.str());

    this->WriteNodeFileHeader(p_node_file);

    // Write each node's data
    for (unsigned index=0; index<this->GetNumNodes(); index++)
    {
        *p_node_file << "Node:\t" << index+1 << "\t";

        for (unsigned i=0; i<DIM; i++)
        {
            *p_node_file << rDeformedPositions[index](i) << "\t";
        }
        *p_node_file << "\n";
    }
    p_node_file->close();
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteCmguiScript(std::string fieldBaseName, std::string undeformedBaseName)
{
    std::string field_string = "";
    if (fieldBaseName != "")
    {
        field_string = " gfx read node " + fieldBaseName + "_$i time $i\n";
    }

    out_stream p_script_file = this->mpOutputFileHandler->OpenOutputFile("LoadSolutions.com");
    *p_script_file << "#\n# Cmgui script automatically generated by Chaste\n#\n";

    if (undeformedBaseName != "")
    {
        *p_script_file << "gfx read node " << undeformedBaseName << " time -1\n";
    }

    *p_script_file << "for ($i=0; $i<=" << mFinalCounter << "; $i++) { \n"
                   << "  gfx read node " << this->mBaseName << "_$i time $i\n"
                   << field_string
                   << "}\n";

    if (undeformedBaseName != "")
    {
        for (unsigned region_index=0; region_index<this->mRegionNames.size(); region_index++)
        {
            *p_script_file << "gfx read ele " << this->mRegionNames[region_index] << "\n";
        }
    }
    else
    {
        *p_script_file << "gfx read ele " << this->mBaseName << "_0\n";
    }
    *p_script_file << "gfx define faces egroup "<<this->mBaseName<<"\n";
    *p_script_file << "gfx modify g_element "<<this->mBaseName<<" lines select_on material default spectrum default selected_material default_selected;\n";
    *p_script_file << "gfx cr win\n\n";
    p_script_file->close();
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::ConvertOutput(std::string inputDirectory,
                                                      std::string inputFileBaseName,
                                                      unsigned finalCounter)
{
    // write the mesh to <inputFileBaseName>_0.exnode and <inputFileBaseName>_0.exelem
    WriteInitialMesh();

    std::vector<c_vector<double,DIM> > deformed_position(mpQuadraticMesh->GetNumNodes(), zero_vector<double>(DIM));

    for (unsigned i=1; i<=finalCounter; i++) //not i=0
    {
        std::stringstream in_file_stream;
        in_file_stream << inputDirectory << "/" << inputFileBaseName << "_" << i << ".nodes";

        std::ifstream ifs(in_file_stream.str().c_str());
        if (!ifs.is_open())
        {
            EXCEPTION("Could not open file: " + in_file_stream.str());
        }

        // the file into deformed_position
        double data;
        for (unsigned index=0; index<mpQuadraticMesh->GetNumNodes(); index++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                ifs >> data;
                if (ifs.fail())
                {
                    EXCEPTION("Error occurred when reading file " << in_file_stream.str()
                                  << ". Expected " << mpQuadraticMesh->GetNumNodes() << " rows and "
                                  << DIM << " columns");
                }
                deformed_position[index](j) = data;
            }
        }

        ifs.close();

        // convert
        WriteDeformationPositions(deformed_position, i);
    }

    WriteCmguiScript();
}


template class CmguiDeformedSolutionsWriter<2>;
template class CmguiDeformedSolutionsWriter<3>;
