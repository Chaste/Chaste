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


#ifndef _TRIANGLESMESHWRITER_HPP_
#define _TRIANGLESMESHWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"
#include "OutputFileHandler.hpp"


/**
 * A concrete mesh writer class that writes Triangle output files.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    TrianglesMeshWriter(const std::string& rDirectory,
                        const std::string& rBaseName,
                        const bool clearOutputDir=true);
    /**
     * Switch this mesh write to write binary files
     *
     * (set to write ascii files in the constructor)
     */
     void SetWriteFilesAsBinary();

    /**
     * Write mesh data to files.
     */
    void WriteFiles();

    /**
     * Write elements as faces (used in the case ELEMENT_DIM== SPACE_DIM-1)
     */
    void WriteElementsAsFaces();

    /**
     * Write faces as edges (used in the case ELEMENT_DIM==2, SPACE_DIM==3)
     */
    void WriteFacesAsEdges();

    /**
     * Write a line (ascii format) to a specific file stream
     * Templated over std::vector dataPacket contents of unsigned or doubles.
     * Templated over type of attribute.
     *
     * @param pFile Pointer to file stream
     * @param itemNumber Index of the element, node or face
     * @param dataPacket List of unsigneds (for node indices) or doubles (for node locations)
     * @param attribute An attribute.  Usually unsigned, but double for cable elements
     */

     template<class  T_DATA, class T_ATTR>
     void WriteItem(out_stream &pFile, unsigned itemNumber, const std::vector<T_DATA> &dataPacket, T_ATTR attribute);

     /**
      * Write a line (ascii format) to a specific file stream
      * Templated over std::vector dataPacket contents of unsigned or doubles.
      *
      * @param pFile Pointer to file stream
      * @param itemNumber Index of the element, node or face
      * @param dataPacket List of unsigneds (for node indices) or doubles (for node locations)
      */
      template<class  T_DATA>
      void WriteItem(out_stream &pFile, unsigned itemNumber, const std::vector<T_DATA> &dataPacket);


    /**
     * Destructor.
     */
    virtual ~TrianglesMeshWriter();
};

#endif //_TRIANGLESMESHWRITER_HPP_
