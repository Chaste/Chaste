//
// Created by bartmanski on 02/05/17.
//

#ifndef CHASTE_IMMERSEDBOUNDARYTESTMESHGENERATOR_HPP
#define CHASTE_IMMERSEDBOUNDARYTESTMESHGENERATOR_HPP

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
#include "SuperellipseGenerator.hpp"
#include "UblasCustomFunctions.hpp"

class ImmersedBoundaryTestMeshGenerator
{
protected:
    /** A pointer to the mesh this class creates. */
    ImmersedBoundaryMesh<2,2>* mpMesh;

public:
    /** Default constructor. */
    ImmersedBoundaryTestMeshGenerator(double cellHeight=0.5, double cellWidth=0.5, unsigned numNodes=100, double ellipseExponent=0.5);

    /** Default destructor */
    virtual ~ImmersedBoundaryTestMeshGenerator();

    /** @return a 2D honeycomb mesh based on a 2D plane. */
    ImmersedBoundaryMesh<2, 2>* GetMesh();

};


#endif //CHASTE_IMMERSEDBOUNDARYTESTMESHGENERATOR_HPP
