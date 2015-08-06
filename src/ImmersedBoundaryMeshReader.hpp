/*

Copyright (c) 2005-2015, University of Oxford.
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
#ifndef IMMERSEDBOUNDARYMESHREADER2D_HPP_
#define IMMERSEDBOUNDARYMESHREADER2D_HPP_

#include <string>
#include <fstream>
#include <cassert>
#include <vector>

#include "Exception.hpp"
#include "AbstractMeshReader.hpp"

/**
 * Helper structure that stores the nodes and any attribute value
 * associated with an ImmersedBoundaryElement.
 */
struct ImmersedBoundaryElementData
{
	std::vector<unsigned> NodeIndices; /**< Vector of Node indices owned by the element. */
	unsigned AttributeValue; /**< Attribute value associated with the element. */
	double SpringConstant; /**< Spring constant associated with the element. */
	double RestLength; /**< Spring rest length associated with the element. */
};

/**
 * A mesh reader class for immersed boundary meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:

	/** The base name for mesh files. */
	std::string mFilesBaseName;

	/** The nodes file for the mesh. */
	std::ifstream mNodesFile;

	/** The elements file for the mesh. */
	std::ifstream mElementsFile;

	/** The grid file for the mesh. */
	std::ifstream mGridFile;

	/** True if input data are numbered from zero, false otherwise. */
	bool mIndexFromZero;

	/** Number of nodes in the mesh. */
	unsigned mNumNodes;

	/** Number of elements in the mesh. */
	unsigned mNumElements;

	/** Number of nodes read in by the reader. */
	unsigned mNodesRead;

	/** Number of elements read in by the reader. */
	unsigned mElementsRead;

	/** Is the number of attributes stored at each node. */
	unsigned mNumNodeAttributes;

	/** Is the number of attributes stored for each element. */
	unsigned mNumElementAttributes;

	/** Is the number of grid points in the x direction. */
	unsigned mNumGridPtsX;

	/** Is the number of grid points in the y direction. */
	unsigned mNumGridPtsY;

	/** The characteristic node spacing. */
	double mCharacteristicNodeSpacing;

	/**
	 * Open node and element files.
	 */
	void OpenFiles();

	/**
	 * Open node file.
	 */
	void OpenNodeFile();

	/**
	 * Open element file.
	 */
	void OpenElementsFile();

	/**
	 * Open element file.
	 */
	void OpenGridFile();

	/**
	 * Read the file headers to determine node and element numbers and attributes.
	 */
	void ReadHeaders();

	/**
	 * Close node and element files.
	 */
	void CloseFiles();

	/**
	 * Get the next line from a given file stream.
	 *
	 * @param fileStream the file stream
	 * @param rawLine the raw line (may contain comments)
	 */
	void GetNextLineFromStream(std::ifstream& fileStream, std::string& rawLine);

public:

	/**
	 * Constructor.
	 *
	 * @param pathBaseName the base name for results files
	 */
	ImmersedBoundaryMeshReader(std::string pathBaseName);

	/**
	 * Destructor.
	 */
	~ImmersedBoundaryMeshReader()
	{}

	/**
	 * @return the number of elements in the mesh.
	 */
	unsigned GetNumElements() const;

	/**
	 * @return the number of nodes in the mesh.
	 */
	unsigned GetNumNodes() const;

	/**
	 * @return the number of grid points in the x direction
	 */
	unsigned GetNumGridPtsX() const;

	/**
	 * @return the number of grid points in the y direction
	 */
	unsigned GetNumGridPtsY() const;

	/**
	 * @return the number of attributes in the mesh
	 */
	unsigned GetNumElementAttributes() const;

	/**
	 * @return the characteristic node spacing in the mesh
	 */
	double GetCharacteristicNodeSpacing();

	/**
	 * Reset pointers to beginning.
	 */
	void Reset();

	/**
	 * @return the coordinates of each node in turn.
	 */
	std::vector<double> GetNextNode();

	/**
	 * @return the next row of the fluid velocity grids.
	 */
	std::vector<double> GetNextGridRow();

	/**
	 * @return the nodes of each element (and any attribute information, if there is any) in turn.
	 */
	ImmersedBoundaryElementData GetNextImmersedBoundaryElementData();

	/**
	 * @return the number of faces in the mesh (synonym of GetNumEdges()).
	 */
	unsigned GetNumFaces() const;

	/**
	 * @return the nodes of each element (and any attribute information, if there is any) in turn.
	 */
	ElementData GetNextElementData();

	/**
	 * @return a vector of the nodes of each face in turn.
	 */
	ElementData GetNextFaceData();

};


#endif /*IMMERSEDBOUNDARYMESHREADER2D_HPP_*/
