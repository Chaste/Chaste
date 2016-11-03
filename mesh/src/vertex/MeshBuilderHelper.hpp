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

#ifndef MESHBUILDERHELPER_HPP_
#define MESHBUILDERHELPER_HPP_

#include "MutableVertexMesh.hpp"
#include "VertexMeshWriter.hpp"

/*
 * nodeIndicesThisElem  array of global node indices which belongs to this to-be-created element in CCW
 * rLowerNodes  all the lower nodes
 * rUpperNodes  maybe all the upper nodes, if nothing here, it will be populated
 * rExistingFaces  all the created lateral faces so that no repeating faces, non-const as newly created face will be push back
 */
class MeshBuilderHelper
{
private:
    std::string mName;
    std::string mAdditionalPath;
    unsigned mNumLowerNodes;
    std::vector<Node<3>*> mLowerNodes;
    std::vector<Node<3>*> mUpperNodes;
    // since number of nodes is known a priori, vector is good enough
    std::vector<std::vector<unsigned> > mNodeToLateralFaceIndices;
    std::vector<VertexElement<2, 3>*> mFaces;
    std::vector<VertexElement<3, 3>*> mElements;
    MutableVertexMesh<3, 3>* mpMesh;
    VertexMeshWriter<3, 3>* mpWriter;

public:
    MeshBuilderHelper(const std::vector<Node<3>*>& rLowerNodes,
                      const std::string& additionalPath = "",
                      const std::string& name = "mesh",
                      const unsigned zHeight = 1);

    MeshBuilderHelper(const std::string& additionalPath = "",
                      const std::string& name = "mesh");

    ~MeshBuilderHelper();

    MutableVertexMesh<3, 3>* MakeMeshUsing2dMesh(const MutableVertexMesh<2, 2>& mesh2, const double zHeight=1);

    MutableVertexMesh<3, 3>* GenerateMesh();

    void PrintMesh(const bool allElements = false) const;

    void WriteVtk(const std::string& outputName, const std::string& additionalTag = "");

    void BuildElementWith(const unsigned numNodesThis, const unsigned nodeIndicesThis[]);

    void BuildElementWith(const std::vector<unsigned>& nodeIndicesThisElem);
};

#endif /*MESHBUILDERHELPER_HPP_*/
