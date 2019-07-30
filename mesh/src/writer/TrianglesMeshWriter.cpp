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

#include "TrianglesMeshWriter.hpp"

#include "AbstractTetrahedralMesh.hpp"
#include "Version.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshWriter(
    const std::string& rDirectory,
    const std::string& rBaseName,
    const bool clearOutputDir)
        : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::~TrianglesMeshWriter()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetWriteFilesAsBinary()
{
    this->mFilesAreBinary=true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "#\n# " + ChasteBuildInfo::GetProvenanceString();
//    PetscTools::Barrier("DodgyBarrierBeforeNODE");
    MeshEventHandler::BeginEvent(MeshEventHandler::NODE);
    // Write node file
    std::string node_file_name = this->mBaseName + ".node";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name, std::ios::binary | std::ios::trunc);

    // Write the node header
    unsigned num_attr = 0;

    if (this->mpMesh && !this->mFilesAreBinary) ///\todo #1949 Readers do not currently support reading of node attributes, so we cannot yet write them from a reader
    {
        // Assumes that all nodes have the same number of attributes as the first node in the mesh.
        num_attr = this->mpMesh->GetNumNodeAttributes();
    }
    ///\todo #1949
    unsigned max_bdy_marker = 0;
    unsigned num_nodes = this->GetNumNodes();

    *p_node_file << num_nodes << "\t";
    *p_node_file << SPACE_DIM << "\t";
    *p_node_file << num_attr << "\t";
    *p_node_file << max_bdy_marker;
    if (this->mFilesAreBinary)
    {
        *p_node_file << "\tBIN\n";
    }
    else
    {
        *p_node_file << "\n";
    }

    *p_node_file << std::setprecision(20);

    // Write each node's data
    for (unsigned item_num=0; item_num<num_nodes; item_num++)
    {
        if (this->mpMesh && !this->mFilesAreBinary && num_attr!=0) ///\todo #1949 Readers do not currently support reading of node attributes, so we cannot yet write them from a reader
        {

            ///\todo #1949 Will deadlock on GetNode(global ID) in parallel since this code is run on the master process
            WriteItem(p_node_file, item_num, this->GetNextNode(), this->mpMesh->GetNode(item_num)->rGetNodeAttributes());
        }
        else
        {
            WriteItem(p_node_file, item_num, this->GetNextNode());
        }
    }
    *p_node_file << comment << "\n";
    p_node_file->close();
//    PetscTools::Barrier("DodgyBarrierAfterNODE");
    MeshEventHandler::EndEvent(MeshEventHandler::NODE);
    MeshEventHandler::BeginEvent(MeshEventHandler::ELE);
    if (ELEMENT_DIM < SPACE_DIM)
    {
        WriteElementsAsFaces();
        WriteFacesAsEdges();
        return;
    }

    // Write element file
    std::string element_file_name = this->mBaseName + ".ele";
    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name, std::ios::binary | std::ios::trunc);

    // Write the element header
    unsigned num_elements = this->GetNumElements();
    num_attr = 1u; // We have a single region code

    // The condition below allows the writer to cope with a NodesOnlyMesh
    if (num_elements == 0)
    {
        *p_element_file << 0 << "\t";
        *p_element_file << 0 << "\t";
        *p_element_file << 0;
        if (this->mFilesAreBinary)
        {
            *p_element_file << "\tBIN\n";
        }
        else
        {
            *p_element_file << "\n";
        }
        p_element_file->close();
    }
    else
    {
        ElementData element_data = this->GetNextElement();

        unsigned nodes_per_element = element_data.NodeIndices.size();
        if (nodes_per_element != ELEMENT_DIM+1)
        {
            // Check that this is a quadratic mesh
            assert(ELEMENT_DIM == SPACE_DIM);
            assert(nodes_per_element == (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2);
         }

        *p_element_file << num_elements << "\t";
        *p_element_file << nodes_per_element << "\t";
        *p_element_file << num_attr;
        if (this->mFilesAreBinary)
        {
            *p_element_file << "\tBIN\n";
        }
        else
        {
            *p_element_file << "\n";
        }

        // Write each element's data
        std::vector<double> attribute_values(1);
        for (unsigned item_num=0; item_num<num_elements; item_num++)
        {
            /*
             * If item_num==0 we will already have got the element above (in order to
             * get the number of nodes per element).
             */
            if (item_num > 0)
            {
                element_data = this->GetNextElement();
            }
            attribute_values[0] =  element_data.AttributeValue;
            WriteItem(p_element_file, item_num, element_data.NodeIndices, attribute_values);
        }
        *p_element_file << comment << "\n";
        p_element_file->close();
    }
//    PetscTools::Barrier("DodgyBarrierAfterELE");
    MeshEventHandler::EndEvent(MeshEventHandler::ELE);
    MeshEventHandler::BeginEvent(MeshEventHandler::FACE);
    // Write boundary face file
    std::string face_file_name = this->mBaseName;

    if (ELEMENT_DIM == 1)
    {
        // In 1-D there is no boundary file: it's trivial to calculate
        return;
    }
    else if (ELEMENT_DIM == 2)
    {
        face_file_name = face_file_name + ".edge";
    }
    else
    {
        face_file_name = face_file_name + ".face";
    }
    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name, std::ios::binary | std::ios::trunc);

    // Write the boundary face header
    if (num_elements != 0)
    {
        unsigned num_faces = this->GetNumBoundaryFaces();

        *p_face_file << num_faces << "\t";
        ///\todo #1949
        *p_face_file << max_bdy_marker;
        if (this->mFilesAreBinary)
        {
            *p_face_file << "\tBIN\n";
        }
        else
        {
            *p_face_file << "\n";
        }

        // Write each face's data
        std::vector<double> default_marker(0);
        for (unsigned item_num=0; item_num<num_faces; item_num++)
        {
            ElementData face_data = this->GetNextBoundaryElement();
            WriteItem(p_face_file, item_num, face_data.NodeIndices, default_marker);
        }
        *p_face_file << comment << "\n";
        p_face_file->close();

        if (this->GetNumCableElements() > 0)
        {
            // Write cable element file
            std::string cable_element_file_name = this->mBaseName + ".cable";
            out_stream p_cable_element_file = this->mpOutputFileHandler->OpenOutputFile(cable_element_file_name, std::ios::binary | std::ios::trunc);

            // Write the cable element header
            unsigned num_cable_elements = this->GetNumCableElements();
            ///\todo #1949
            num_attr = 1u; // We have a single region code - which is actually a radius

            *p_cable_element_file << num_cable_elements << "\t";
            *p_cable_element_file << 2 << "\t";
            *p_cable_element_file << num_attr;
            if (this->mFilesAreBinary)
            {
                *p_cable_element_file << "\tBIN\n";
            }
            else
            {
                *p_cable_element_file << "\n";
            }

            // Write each element's data
            std::vector<double> attribute_values(1);
            for (unsigned item_num=0; item_num<num_cable_elements; item_num++)
            {
                ElementData cable_element_data = this->GetNextCableElement();
                attribute_values[0] = cable_element_data.AttributeValue;
                WriteItem(p_cable_element_file, item_num, cable_element_data.NodeIndices, attribute_values);
            }
            *p_cable_element_file << comment;
            p_cable_element_file->close();
        }
    }
//    PetscTools::Barrier("DodgyBarrierAfterFACE");
    MeshEventHandler::EndEvent(MeshEventHandler::FACE);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteElementsAsFaces()
{
    std::string comment = "#\n# " + ChasteBuildInfo::GetProvenanceString();

    std::string element_file_name = this->mBaseName;
    if (ELEMENT_DIM == 1 && (SPACE_DIM == 2 || SPACE_DIM == 3))
    {
        element_file_name = element_file_name + ".edge";
    }
    else if (ELEMENT_DIM == 2 && SPACE_DIM == 3)
    {
        element_file_name = element_file_name + ".face";
    }

    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name, std::ios::binary | std::ios::trunc);

    // Write the element header
    unsigned num_elements = this->GetNumElements();
    assert(SPACE_DIM != ELEMENT_DIM);    // LCOV_EXCL_LINE
    unsigned num_attr = 1;

    *p_element_file << num_elements << "\t";
    //*p_element_file << nodes_per_element << "\t";
    *p_element_file << num_attr;
    if (this->mFilesAreBinary)
    {
        *p_element_file << "\tBIN\n";
    }
    else
    {
        *p_element_file << "\n";
    }

    // Write each element's data
    std::vector<double> attribute_values(1);
    for (unsigned item_num=0; item_num<num_elements; item_num++)
    {
        ElementData element_data = this->GetNextElement();
        attribute_values[0] =  element_data.AttributeValue;
        WriteItem(p_element_file, item_num, element_data.NodeIndices, attribute_values);
    }

    *p_element_file << comment << "\n";
    p_element_file->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFacesAsEdges()
{
    std::string comment = "#\n# " + ChasteBuildInfo::GetProvenanceString();

    if (ELEMENT_DIM == 1 && (SPACE_DIM == 2 || SPACE_DIM == 3))
    {
        return;
    }

    assert(SPACE_DIM == 3 && ELEMENT_DIM == 2);    // LCOV_EXCL_LINE

    std::string face_file_name = this->mBaseName;
    face_file_name = face_file_name + ".edge";

    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name, std::ios::binary | std::ios::trunc);

    // Write the boundary face header
    unsigned num_faces = this->GetNumBoundaryFaces();

    unsigned max_bdy_marker = 0;
    std::vector<double> default_marker(0);

    *p_face_file << num_faces << "\t";
    *p_face_file << max_bdy_marker;
    if (this->mFilesAreBinary)
    {
        *p_face_file << "\tBIN\n";
    }
    else
    {
        *p_face_file << "\n";
    }

    // Write each face's data
    for (unsigned item_num=0; item_num<num_faces; item_num++)
    {
        ElementData face_data = this->GetNextBoundaryElement();
        WriteItem(p_face_file, item_num, face_data.NodeIndices, default_marker);
    }
    *p_face_file << comment << "\n";
    p_face_file->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template<class T_DATA>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteItem(out_stream &pFile, unsigned itemNumber,
                                                            const std::vector<T_DATA> &dataPacket)
{
    //Writing with no attribute
    //Instantiates the attribute variety with the attributes empty
    WriteItem(pFile, itemNumber, dataPacket, std::vector<double>());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template<class T_DATA>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteItem(out_stream &pFile, unsigned itemNumber,
                                                            const std::vector<T_DATA> &dataPacket,
                                                            const std::vector<double> &rAttributes)
{
    if (this->mFilesAreBinary)
    {
        // No item numbers
        // Write raw data out of std::vector into the file
        pFile->write((char*)&dataPacket[0], dataPacket.size()*sizeof(T_DATA));

        // Write raw attribute(s)
        // MSVC 10 will trip an assertion on accessing the 0th element if the vector is empty.
        // Note that C++11 gives vectors a .data() method which would be a better way of doing this!
        if (!rAttributes.empty())
        {
            pFile->write((char*)&rAttributes[0], rAttributes.size()*sizeof(double));
        }
    }
    else
    {
        *pFile << itemNumber;
        for (unsigned i=0; i<dataPacket.size(); i++)
        {
            *pFile << "\t" << dataPacket[i];
        }
        for (unsigned i=0; i<rAttributes.size(); i++)
        {
            *pFile << "\t" << rAttributes[i];
        }
        *pFile << "\n";
    }
}

// Explicit instantiation
template class TrianglesMeshWriter<1,1>;
template class TrianglesMeshWriter<1,2>;
template class TrianglesMeshWriter<1,3>;
template class TrianglesMeshWriter<2,2>;
template class TrianglesMeshWriter<2,3>;
template class TrianglesMeshWriter<3,3>;

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void TrianglesMeshWriter<1, 1>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<1, 1>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
template void TrianglesMeshWriter<1, 2>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<1, 2>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
template void TrianglesMeshWriter<1, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<1, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
template void TrianglesMeshWriter<2, 2>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<2, 2>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
template void TrianglesMeshWriter<2, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<2, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
template void TrianglesMeshWriter<3, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, const std::vector<double>   & );
template void TrianglesMeshWriter<3, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, const std::vector<double>   & );
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
