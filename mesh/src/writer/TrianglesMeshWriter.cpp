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

#include "TrianglesMeshWriter.hpp"
#include "Version.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshWriter(
    const std::string &rDirectory,
    const std::string &rBaseName,
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

    // Write node file
    std::string node_file_name = this->mBaseName + ".node";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    // Write the node header
    unsigned num_attr = 0;
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
    unsigned default_marker = UINT_MAX; ///\todo #1899 or #1949 Is this necessary?
    for (unsigned item_num=0; item_num<num_nodes; item_num++)
    {
        WriteItem(p_node_file, item_num, this->GetNextNode(), default_marker);
    }
    *p_node_file << comment << "\n";
    p_node_file->close();

    if (ELEMENT_DIM < SPACE_DIM)
    {
        WriteElementsAsFaces();
        WriteFacesAsEdges();
        return;
    }

    // Write element file
    std::string element_file_name = this->mBaseName + ".ele";
    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    unsigned num_elements = this->GetNumElements();
    ///\todo #1949
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
            //Note: Cable element has a double attribute for radius but regular element have always had an unsigned region marker
            WriteItem(p_element_file, item_num, element_data.NodeIndices, (unsigned) element_data.AttributeValue);
        }
        *p_element_file << comment << "\n";
        p_element_file->close();
    }

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
    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name);

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
        default_marker = UINT_MAX;
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
            out_stream p_cable_element_file = this->mpOutputFileHandler->OpenOutputFile(cable_element_file_name);

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
            for (unsigned item_num=0; item_num<num_cable_elements; item_num++)
            {
                ElementData cable_element_data = this->GetNextCableElement();
                //Cable element has a double attribute for radius
                WriteItem(p_cable_element_file, item_num, cable_element_data.NodeIndices, cable_element_data.AttributeValue);
            }
            *p_cable_element_file << comment;
            p_cable_element_file->close();
        }
    }
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

    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    unsigned num_elements = this->GetNumElements();
    assert(SPACE_DIM != ELEMENT_DIM);
    unsigned num_attr = 0;

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
    for (unsigned item_num=0; item_num<num_elements; item_num++)
    {
         WriteItem(p_element_file, item_num, this->GetNextElement().NodeIndices);
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

    assert(SPACE_DIM == 3 && ELEMENT_DIM == 2);

    std::string face_file_name = this->mBaseName;
    face_file_name = face_file_name + ".edge";

    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name);

    // Write the boundary face header
    unsigned num_faces = this->GetNumBoundaryFaces();

    unsigned max_bdy_marker = 0;
    unsigned default_marker = UINT_MAX;

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
    //Instantiates the attribute variety with the attribute type set to unsigned
    WriteItem(pFile, itemNumber, dataPacket, UINT_MAX);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template<class T_DATA, class T_ATTR>
void TrianglesMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteItem(out_stream &pFile, unsigned itemNumber,
                                                            const std::vector<T_DATA> &dataPacket, T_ATTR attribute)
{
    if (this->mFilesAreBinary)
    {
        // No item numbers
        // Write raw data out of std::vector into the file
        pFile->write((char*)&dataPacket[0], dataPacket.size()*sizeof(T_DATA));

        // Write raw attribute
        if (attribute != (std::numeric_limits<T_ATTR>::max)())  //or attribute != UINT_MAX
        {
            pFile->write((char*) &attribute, sizeof(attribute));
        }
    }
    else
    {
        *pFile << itemNumber;
        for (unsigned i=0; i<dataPacket.size(); i++)
        {
            *pFile << "\t" << dataPacket[i];
        }
        if (attribute != (std::numeric_limits<T_ATTR>::max)())  //or attribute != UINT_MAX
        {
            *pFile << "\t" << attribute;
        }
        *pFile << "\n";
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

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
template void TrianglesMeshWriter<1, 1>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<1, 1>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
template void TrianglesMeshWriter<1, 2>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<1, 2>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
template void TrianglesMeshWriter<1, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<1, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
template void TrianglesMeshWriter<2, 2>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<2, 2>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
template void TrianglesMeshWriter<2, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<2, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
template void TrianglesMeshWriter<3, 3>::WriteItem(out_stream &, unsigned, const std::vector<unsigned> &, double );
template void TrianglesMeshWriter<3, 3>::WriteItem(out_stream &, unsigned, const std::vector<double>   &, double );
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
