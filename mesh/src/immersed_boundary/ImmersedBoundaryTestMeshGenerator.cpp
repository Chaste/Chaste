//
// Created by bartmanski on 02/05/17.
//

#include "ImmersedBoundaryTestMeshGenerator.hpp"

ImmersedBoundaryTestMeshGenerator::ImmersedBoundaryTestMeshGenerator(double cellHeight,
                                                                     double cellWidth,
                                                                     unsigned numNodes,
                                                                     double ellipseExponent) : mpMesh(NULL)
{
    // Check for sensible input
    assert(cellHeight > 0.0);
    assert(cellWidth > 0.0);
    assert(numNodes > 3);
    assert(ellipseExponent > 0.0);

    // Helper vectors
    unit_vector<double> x_unit(2, 0);
    unit_vector<double> y_unit(2, 1);

    // Set up a reference cell
    SuperellipseGenerator* p_gen = new SuperellipseGenerator(numNodes, ellipseExponent, cellWidth, cellHeight,
                                                             0.5*(1.0 - cellWidth), 0.5*(1.0 - cellHeight));
    std::vector<c_vector<double, 2> > locations =  p_gen->GetPointsAsVectors();

    // Set up the containers for holding the nodes and ib_elements
    std::vector<Node<2>*> nodes;
    std::vector<ImmersedBoundaryElement<2, 2>*> ib_elements;
    std::vector<ImmersedBoundaryElement<1, 2>*> ib_laminas;

    // Shift the reference cell to an appropriate place
    for (unsigned location_index = 0; location_index < locations.size(); location_index++)
    {
        Node<2>* p_node = new Node<2>(nodes.size(), locations[location_index], true);
        nodes.push_back(p_node);
    }
    ib_elements.push_back(new ImmersedBoundaryElement<2, 2>(ib_elements.size(), nodes));

    mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, ib_laminas);
}

ImmersedBoundaryTestMeshGenerator::~ImmersedBoundaryTestMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2, 2>* ImmersedBoundaryTestMeshGenerator::GetMesh()
{
    return mpMesh;
}