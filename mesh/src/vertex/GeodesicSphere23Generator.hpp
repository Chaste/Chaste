/*
 * GeodesicSphere23MeshGenerator.hpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Weijie
 */

#ifndef GEODESICSPHERE23GENERATOR_HPP_
#define GEODESICSPHERE23GENERATOR_HPP_

#include <set>
#include <vector>
#include <utility>              // for pair
#include "Node.hpp"
#include "VertexElement.hpp"
#include "MutableVertexMesh.hpp"


class GeodesicSphere23Generator
{
private:
    unsigned mDepth;
    static const double X = 0.525731112119133606;
    static const double Z = 0.850650808352039932;

    class SortWithIndex
    {
    public:
        bool operator()(const std::pair<double, unsigned>& a, const std::pair<double, unsigned>& b) const
        {
            return a.first < b.first;
        }
    };

public:
    std::vector<Node<3>*> mNodes;
    std::vector<VertexElement<1,3>*> mEdges;
    std::vector<VertexElement<2,3>*> mFaces;

    c_vector<double,3> normalise(const c_vector<double,3>& v)
    {
        return v/norm_2(v);
    }

    GeodesicSphere23Generator(const unsigned numDivision=0u);

    void SubDivide();

    MutableVertexMesh<2, 3>* GetDual();
};

#endif /* GEODESICSPHERE23GENERATOR_HPP_ */
