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


#ifndef CMGUIDEFORMEDSOLUTIONSWRITER_HPP_
#define CMGUIDEFORMEDSOLUTIONSWRITER_HPP_

#include "CmguiMeshWriter.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "DistributedQuadraticMesh.hpp"

/**
 *  Small enumeration for representing whether we want linear
 *  visualisation of the quadratic mesh (just vertices output)
 *  or full quadratic visualisation
 */
typedef enum CmguiMeshWriteType_
{
    WRITE_LINEAR_MESH = 0,
    WRITE_QUADRATIC_MESH
} CmguiMeshWriteType;

/**
 *  CmguiDeformedSolutionsWriter
 *
 *  A class for writing a mesh, and solutions from a solid mechanics (ie deformed meshes) problem,
 *  in Cmgui output format.
 *
 *  Inherits from CmguiMeshWriter
 */
template<unsigned DIM>
class CmguiDeformedSolutionsWriter : public CmguiMeshWriter<DIM, DIM>
{
private:
    /**
     *  The quadratic mesh used in the mechanics simulation. solution_0.exnode and solution_0.exelem
     *  (the only exelem file written) will be written using this mesh
     */
    AbstractTetrahedralMesh<DIM,DIM>* mpQuadraticMesh;

    /**
     *  A counter is given whenever WriteDeformationPositions() is called, this variable
     *  stores the last one used
     */
    unsigned mFinalCounter;

    /**
     *  Number of nodes to output - either mpQuadraticMesh->GetNumVertices() (linear visualisation)
     *  mpQuadraticMesh->GetNumNodes() (quadratic visualisation)
     */
    unsigned mNumNodesToUse;

    /**
     *  Overloaded GetNumNodes().
     *  @return either mpQuadraticMesh->GetNumVertices() for linear
     *  visualisation or mpQuadraticMesh->GetNumNodes() for quadratic visualisation
     */
    unsigned GetNumNodes()
    {
        return mNumNodesToUse;
    }

public:
    /**
     *  Constructor
     *  @param outputDirectory The output directory for the Cmgui files
     *  @param baseName The base name for the Cmgui output files - the files written will be
     *   [basename_0.exnode, [basename]_0.exelem; [basename]_1.exnode, [basename]_2.exnode, ..
     *  @param rQuadraticMesh The quadratic mesh used in the mechanics simulation
     *  @param writeType Should be equal to either WRITE_LINEAR_MESH or WRITE_QUADRATIC_MESH,
     *   depending on whether linear visualisation of the quadratic mesh (just vertices output)
     *   or full quadratic visualisation is required.
     */
    CmguiDeformedSolutionsWriter(std::string outputDirectory,
                                 std::string baseName,
                                 AbstractTetrahedralMesh<DIM,DIM>& rQuadraticMesh,
                                 CmguiMeshWriteType writeType);

    /**
     *  Write [basename]_0.exnode, [basename]_0.exelem using the quadratic mesh
     *  If the optional argument fileName is given, writes [fileName].exnode
     *  and [fileName].exelem instead
     *
     *  @param fileName Optional file name (stem).
     */
    void WriteInitialMesh(std::string fileName = "");

    /**
     *  Write [basename]_i.exnode using the given deformed positions
     *  @param rDeformedPositions std::vector of deformed positions to be used, must have size equal to number
     *  of nodes in the mesh
     *  @param counter the value "i" in "[basename]_i.exnode" to be used.
     */
    void WriteDeformationPositions(std::vector<c_vector<double,DIM> >& rDeformedPositions,
                                   unsigned counter);

    /**
     *  Writes a small cmgui script called LoadSolutions.com, for loading the output that has been written.
     *  Assumes the output was solution_0.exnode .. solution_N.exnode, where N is the counter that was
     *  given in the last call to WriteDeformationPositions()
     *
     *  @param fieldBaseName If there is a field to visualise on top of the deforming mesh, give it's
     *    path (relative to the cmgui deformation directory), and filename prefix. Leave empty
     *    if no field to visualise.
     *    For example,
     *      WriteCmguiScript("../../electrics/cmgui_output/voltage_mechanics_mesh");
     *    for the script to read files voltage_mechanics_mesh_0.exnode .. voltage_mechanics_mesh_N.exnode.
     *
     *  @param undeformedBaseName is assumed to be "solution_0.exnode" and .exelem unless this optional
     *    parameter is given. Depends on what parameters were given to WriteInitialMesh(). If
     *    undeformedBaseName the time given to plot the undeformed shape is set to -1, otherwise
     *    it is set to 0.
     */
    void WriteCmguiScript(std::string fieldBaseName="", std::string undeformedBaseName="");

    /**
     *  For a simulation that has already been run, convert the chaste output to cmgui format.
     *  @param inputDirectory The directory the chaste output is in
     *  @param inputFileBaseName The base name for the chaste output
     *  @param finalCounter   The final counter, ie the value N for the final file [basename]_N.nodes.
     *  The first file is assumed to be [basename]_0.nodes. The files [basename]_0.nodes up to
     *  [basename]_N.nodes will be converted to cmgui format and put in the output directory
     *  given in the constructor.
     */
    void ConvertOutput(std::string inputDirectory,
                       std::string inputFileBaseName,
                       unsigned finalCounter);
};

#endif /*CMGUIDEFORMEDSOLUTIONSWRITER_HPP_*/

