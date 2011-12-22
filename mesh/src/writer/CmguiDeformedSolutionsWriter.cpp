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

#include "CmguiDeformedSolutionsWriter.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CmguiDeformedSolutionsWriter<DIM>::CmguiDeformedSolutionsWriter(std::string outputDirectory,
                                                                std::string baseName,
                                                                QuadraticMesh<DIM>& rQuadraticMesh,
                                                                CmguiMeshWriteType writeType)
    : CmguiMeshWriter<DIM, DIM>(outputDirectory, baseName),
      mpQuadraticMesh(&rQuadraticMesh),
      mFinalCounter(0)
{

    mNumNodesToUse = mpQuadraticMesh->GetNumVertices();

    if (writeType==WRITE_QUADRATIC_MESH)
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
void CmguiDeformedSolutionsWriter<DIM>::WriteInitialMesh()
{
    std::string saved_base_name = this->mBaseName;
    this->mBaseName = this->mBaseName + "_0";
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
void CmguiDeformedSolutionsWriter<DIM>::WriteCmguiScript(std::string fieldBaseName)
{
    std::string field_string = "";
    if (fieldBaseName != "")
    {
        field_string = " gfx read node " + fieldBaseName + "_$i time $i\n";
    }

    out_stream p_script_file = this->mpOutputFileHandler->OpenOutputFile("LoadSolutions.com");
    *p_script_file << "#\n# Cmgui script automatically generated by Chaste\n#\n"
                   << "for ($i=0; $i<=" << mFinalCounter << "; $i++) { \n"
                   << "  gfx read node " << this->mBaseName << "_$i time $i\n"
                   << field_string
                   << "}\n"
                   << "gfx read ele " << this->mBaseName << "_0 generate_faces_and_lines\n"
                   << "gfx cr win\n\n";
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
