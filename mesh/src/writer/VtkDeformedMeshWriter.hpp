/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef VTKDEFORMEDMESHWRITER_HPP_
#define VTKDEFORMEDMESHWRITER_HPP_

#ifdef CHASTE_VTK

#include <cstring> 
#include "VtkMeshWriter.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * A class to write meshes in VTK format, allowing for specification of 
 * a new set of node positions ("deformed mesh") and a new file name.
 * 
 * This class was designed with time-dependent simulations in mind, where,
 * at every time step, one may want to write out a mesh that is being
 * deformed in the time loop.
 * 
 * Typical usage:
 *
 *  VtkDeformedMeshWriter<DIM> writer(...);
 *  //write undeformed mesh
 *  writer.AddPointData("name of node-wise data", some_data);
 *  writer.AddCellData("name of element-wise data", other_data);
 *  writer.WriteDeformedFiles();
 *
 *  //apply deformation - lines below could be in a time loop
 *  writer.ApplyDeformation(new_node_locations);
 *  writer.SetWriteMeshCells(false); //cells have been alreday written
 *  writer.SetOutputBaseFileName(new_name);//do not overwrite the previous file
 *  writer.WriteDeformedFiles();
 * 
 *  Note that the ApplyDeformation method actually changes
 *  the mesh node coordinates. Use with caution as no other
 *  quantities of the mesh (e.g., Jacobian) are actually modified!      
 */
template<unsigned DIM>
class VtkDeformedMeshWriter : public VtkMeshWriter<DIM,DIM>
{

private:

    /** 
     * Pointer to the deformed mesh, used only for output here. 
     * Initialized to a new mesh object and constructed upon construction of this writer object
     */
    AbstractTetrahedralMesh<DIM,DIM>* mpDeformedMesh;

public:

    /**
     * Constructor
     * 
     * @param pDeformedOutputMesh a pointer to the mesh objectto be used.
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param rCleanDirectory  whether to clean the directory (defaults to true)
     */
    VtkDeformedMeshWriter(AbstractTetrahedralMesh<DIM,DIM>* pDeformedOutputMesh, const std::string& rDirectory, const std::string& rBaseName, const bool& rCleanDirectory=true);

    /**
     * Destructor.
     */
    ~VtkDeformedMeshWriter();

    /**
     * Set the new file name to be used the next time WriteDeformedFiles() is called.
     * The output file will then be rOutputBaseFileName.vtu
     * 
     * @param rOutputBaseFileName the base name of the next output file (".vtu" will be added)
     */
    void SetOutputBaseFileName(const std::string& rOutputBaseFileName);


    /**
     * Write the VTK files using the internal mesh, with any deformation that has been applied.
     * Internally, it calls WriteFilesUsingMesh().
     * The output file name will be based on the latest call to SetOutputBaseFileName,
     * or, if SetOutputBaseFileName was not called, the constructor parameter rBaseName
     */
    void WriteDeformedFiles();

    /**
     * Modify the node locations of the mesh by assigning the coordinates
     * according to rPositions. Note that only node coordinates are changed 
     * (no Jacobians or other related quantities are modified). Use with caution!
     * 
     * @param rPositions the node positions to be applied to the mesh
     */
    void ApplyDeformation(const std::vector<c_vector<double,DIM> >& rPositions);

};

#endif //CHASTE_VTK

#endif /*VTKDEFORMEDMESHWRITER_HPP_*/

