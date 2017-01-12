/*
 * GeodesicSphere23MeshGenerator.hpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Weijie
 */

#ifndef GEODESICSPHERE23GENERATOR_HPP_
#define GEODESICSPHERE23GENERATOR_HPP_

#include <vector>
#include "Node.hpp"
#include "VertexElement.hpp"
#include "MutableVertexMesh.hpp"
#include "UblasCustomFunctions.hpp"


class GeodesicSphere23Generator
{
    friend class TestMonolayerVertexMeshGenerator;
    
private:
    unsigned mDepth;
    std::vector<Node<3>*> mNodes;
    std::vector<VertexElement<1,3>*> mEdges;
    std::vector<VertexElement<2,3>*> mFaces;

public:
    
    c_vector<double,3> normalise(const c_vector<double,3>& v)
    {
        return v/norm_2(v);
    }

    GeodesicSphere23Generator(const unsigned numDivision=0u);

    void SubDivide();

    MutableVertexMesh<2, 3>* GetDual();
};

#endif /* GEODESICSPHERE23GENERATOR_HPP_ */
