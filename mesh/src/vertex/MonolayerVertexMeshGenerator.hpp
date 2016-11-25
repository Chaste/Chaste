/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef MONOLAYERVERTEXMESHGENERATOR_HPP_
#define MONOLAYERVERTEXMESHGENERATOR_HPP_

#include "MutableVertexMesh.hpp"
#include "VertexMeshWriter.hpp"
#include <set>

/**
 * A class that generates a monolayer vertex mesh with given 2D mesh or sufficient
 * information about the mesh (elements and their nodes).
 */
class MonolayerVertexMeshGenerator
{
protected:

    /**
     * Name of this instance. When writing vtk file, the name will be added at
     * the very beginning.
     */
    std::string mName;

    /**
     * A vector storing all the basal nodes of the mesh.
     */
    std::vector<Node<3>*> mBasalNodes;

    /**
     * A vector storing all the apical nodes of the mesh.
     */
    std::vector<Node<3>*> mApicalNodes;

    /**
     * A helper variables that tell to which lateral faces does a node belongs
     * so that face will not be created twice during mesh-generation.
     * Since the number of nodes is known a-priori, we can use a vector of sets
     * for this map.
     * This couldn't not be a local variable of a function as element might be
     * generated seperately.
     */
    std::vector<std::set<unsigned> > mNodeLateralFaceMap;

    /**
     * A vector storing generated faces.
     */
    std::vector<VertexElement<2, 3>*> mFaces;

    /**
     * A vector storing generated elements.
     */
    std::vector<VertexElement<3, 3>*> mElements;

    /**
     * A pointer for generated mesh.
     */
    MutableVertexMesh<3, 3>* mpMesh;

    /**
     * A pointer for mesh writer.
     */
    VertexMeshWriter<3, 3>* mpWriter;

public:

    /**
     * Default constructor of this class. It does basically nothing.
     *
     * @param name  the name given to this instance.
     */
    MonolayerVertexMeshGenerator(const std::string& name = "mesh");

    /**
     * Constructor that will take basal nodes and create corresponding
     * apical nodes.
     *
     * @param rBasalNodes  basal nodes of the to-be-generated-mesh.
     * @param name  the name given to this instance.
     * @param zHeight  vertical distance (in z-direction) between apical
     *                 and basal nodes. Default is set to 1.
     */
    MonolayerVertexMeshGenerator(const std::vector<Node<3>*>& rBasalNodes,
                                 const std::string& name = "mesh",
                                 const unsigned zHeight = 1);

    /**
     * Destructor.
     */
    ~MonolayerVertexMeshGenerator();

    /**
     * Creates a monolayer vertex mesh from a 2D mesh.
     *
     * @param mesh2d  the given 2D vertex mesh
     * @param zHeight  vertical distance (in z-direction) between apical
     *                 and basal nodes. Default is set to 1.
     * @return  a pointer of the monolayer vertex mesh.
     */
    MutableVertexMesh<3, 3>* MakeMeshUsing2dMesh(const MutableVertexMesh<2, 2>& mesh2d,
                                                 const double zHeight=1);

    /**
     * Builds individual element with given basal node indices.
     * Note that basal nodes should have been created by the constructor.
     *
     * @param basalNodeIndices  indices of basal nodes in CCW.
     */
    void BuildElementWith(const std::vector<unsigned>& basalNodeIndices);

    /**
     * A wrapper function serves as an alternative because c++98 doesn't support
     * the use of initializer list for vector constructor.
     *
     * @param numBasalNodes  number of basal nodes of this element.
     * @param basalNodeIndices  indices of basal nodes in CCW.
     */
    void BuildElementWith(const unsigned numBasalNodes, const unsigned basalNodeIndices[]);

    /**
     * Generates the complete mesh after all the elements are generated.
     *
     * @return  a pointer of the monolayer vertex mesh.
     */
    MutableVertexMesh<3, 3>* GenerateMesh();

    /**
     * Clear all mesh objects stored within the class.
     */
    void ClearStoredMeshObjects();

    /**
     * Output a .vtu data for visualization using Paraview.
     * @param outputFile  the directory of the output data
     * @param additionalTag  optional tag which will be append at the end of the data name
     * @param usingFaceId  whether to write mesh with face ID on it. By default writing with
     *                     element ID.
     */
    void WriteVtk(const std::string& outputFile, const std::string& additionalTag = "",
                  const bool usingFaceId = false);

    /**
     * Alternative to output a .vtu data for visualization using Paraview. A subfolder which
     * has the same name as this class will be created to have a tidier file system.
     * @param outputFile  the directory of the output data
     * @param additionalTag  optional tag which will be append at the end of the data name
     * @param usingFaceId  whether to write mesh with face ID on it. By default writing with
     *                     element ID.
     */
    void WriteVtkWithSubfolder(const std::string& outputFile, const std::string& additionalTag = "",
                               const bool usingFaceId = false);

    /**
     * Prints out the elements, nodes and faces in verbose mode.
     * Mainly for debugging purpose.
     *
     * @param printDeletedObjects  whether to include deleted objects (false by default)
     */
    void PrintMesh(const bool printDeletedObjects = false) const;
};

#endif /*MONOLAYERVERTEXMESHGENERATOR_HPP_*/
