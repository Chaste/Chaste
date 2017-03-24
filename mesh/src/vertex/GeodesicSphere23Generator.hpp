/*
 * GeodesicSphere23MeshGenerator.hpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Weijie
 */

#ifndef GEODESICSPHERE23GENERATOR_HPP_
#define GEODESICSPHERE23GENERATOR_HPP_

#include <vector>
#include "MutableVertexMesh.hpp"
#include "Node.hpp"
#include "VertexElement.hpp"

class GeodesicSphere23Generator
{
    friend class TestMonolayerVertexMeshGenerator;

protected:
    std::vector<Node<3>*> mNodes;
    std::vector<VertexElement<1, 3>*> mEdges;
    std::vector<VertexElement<2, 3>*> mFaces;

public:
    GeodesicSphere23Generator();

    void SubDivide();

    MutableVertexMesh<2, 3>* GetDual();
};

#endif /* GEODESICSPHERE23GENERATOR_HPP_ */
