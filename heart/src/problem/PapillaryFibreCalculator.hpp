/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef PAPILLARYFIBRECALCULATOR_HPP_
#define PAPILLARYFIBRECALCULATOR_HPP_

#include "UblasCustomFunctions.hpp"
#include "TetrahedralMesh.hpp"

/**
 * Assigns fibre orientation vectors for papillary muscle structures. Vectors are
 * assigned to be parallel to the direction of the axis of the muscle. To do this
 * we use the "Structure Tensor" method.
 */
class PapillaryFibreCalculator
{
// Allow the test class to use the private functions.
friend class TestPapillaryFibreCalculator;

private:
    /** Reference to the ventricular/papillary mesh*/
    TetrahedralMesh<3,3>& mrMesh;
    /** vectors from the centre of each element to the nearest boundary node */
    std::vector< c_vector<double, 3> > mRadiusVectors;
    /** Tensors created from mRadiusVectors */
    std::vector< c_matrix<double,3,3> > mStructureTensors;
    /** Smoothed tensors created from mStructureTensors */
    std::vector< c_matrix<double,3,3> > mSmoothedStructureTensors;

   /**
     * This method calculates the vector from the centroid of an element to all of
     * the boundary nodes. It returns the shortest of the vectors.
     *
     * @param elementIndex  The index of the element we are calculating radial vectors for
     * @return The shortest radial vector
     */
    c_vector<double,3> GetRadiusVectorForOneElement(unsigned elementIndex);

    /**
     * This method calls GetRadiusVectorForOneElement() for each of the elements sets
     * the radial vector for each element of the mesh.
     */
    void GetRadiusVectors();


    /**
     * This generates structure tensors from the radial vectors by taking
     *
     * T = r.r'
     */
    void ConstructStructureTensors();

    /**
     * Smoothes the structure tensor components for each papillary element by looping
     * over all other papillary elements, calculating
     * distance geometric distance between the two elements;
     * if it is within a certain limit, include this in the Gaussian kernel
     *
     * Here mStructureTensors[i] is the 'rough' tensor for each element
     * and mSmoothedStructureTensors[i] is the smoothed tensor for each element
     */
    void SmoothStructureTensors();

public:
    /**
     * Constructor saves mesh and allocates memory
     *
     * @param rMesh  The mesh to calculate fibres on
     */
    PapillaryFibreCalculator(TetrahedralMesh<3,3>& rMesh);

    /**
     *  Main method - calculate the fibre orientations
     *
     *  @return A fibre vector for each element
     */
     std::vector<c_vector<double,3> > CalculateFibreOrientations();

};

// PUBLIC METHODS
PapillaryFibreCalculator::PapillaryFibreCalculator(TetrahedralMesh<3,3>& rMesh)
    : mrMesh(rMesh)
{
    mRadiusVectors.resize(mrMesh.GetNumElements());
    mStructureTensors.resize(mrMesh.GetNumElements());
    mSmoothedStructureTensors.resize(mrMesh.GetNumElements());
}

std::vector<c_vector<double,3> > PapillaryFibreCalculator::CalculateFibreOrientations()
{
   GetRadiusVectors();

   ConstructStructureTensors();

   SmoothStructureTensors();

   // Calculate eigenvalues
   std::vector<c_vector<double,3> > fibre_orientations(mrMesh.GetNumElements());
   for(unsigned i=0; i<fibre_orientations.size(); i++)
   {
       fibre_orientations[i] = CalculateEigenvectorForSmallestNonzeroEigenvalue(mSmoothedStructureTensors[i]);
   }

   return fibre_orientations;
}

// PRIVATE METHODS
c_vector<double,3> PapillaryFibreCalculator::GetRadiusVectorForOneElement(unsigned elementIndex)
{
    c_vector<double, 3> centroid = (mrMesh.GetElement(elementIndex))->CalculateCentroid();
    // Loops over all papillary face nodes
    c_vector<double,3> coordinates;

    double nearest_r_squared=DBL_MAX;
    unsigned nearest_face_node = 0;

    TetrahedralMesh<3,3>::BoundaryNodeIterator bound_node_iter = mrMesh.GetBoundaryNodeIteratorBegin();
    while (bound_node_iter != mrMesh.GetBoundaryNodeIteratorEnd())
    {
        unsigned bound_node_index =  (*bound_node_iter)->GetIndex();
        coordinates=mrMesh.GetNode(bound_node_index)->rGetLocation();

        // Calculates the distance between the papillary face node and the centroid
        double r_squared =  norm_2(centroid-coordinates);
        // Checks to see if it is the smallest so far - if it is, update the current smallest distance
        if (r_squared < nearest_r_squared)
        {
            nearest_r_squared = r_squared;
            nearest_face_node = bound_node_index;
        }
        ++bound_node_iter;
    }

    coordinates = mrMesh.GetNode(nearest_face_node)->rGetLocation();
    c_vector<double,3> radial_vector = coordinates-centroid;
    return radial_vector;
}

void PapillaryFibreCalculator::GetRadiusVectors()
{
    // Loops over all elements finding radius vector
    for (AbstractTetrahedralMesh<3,3>::ElementIterator iter = mrMesh.GetElementIteratorBegin();
         iter != mrMesh.GetElementIteratorEnd();
         ++iter)
    {
        unsigned element_index = iter->GetIndex();
        mRadiusVectors[element_index] = GetRadiusVectorForOneElement(element_index);
    }
}

void PapillaryFibreCalculator::ConstructStructureTensors()
{
    for(unsigned i=0;i<mRadiusVectors.size();i++)
    {
        mStructureTensors[i] = outer_prod(mRadiusVectors[i],mRadiusVectors[i]);
    }
}

void PapillaryFibreCalculator::SmoothStructureTensors()
{
    const double sigma = 0.05; //cm
    const double r_max = 0.1; //cm
    double g_factor_sum = 0;
    double g_factor = 0;

    for (AbstractTetrahedralMesh<3,3>::ElementIterator elem_iter = mrMesh.GetElementIteratorBegin();
         elem_iter != mrMesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        mSmoothedStructureTensors[ elem_iter->GetIndex()] = zero_matrix<double>(3,3);

        c_vector<double, 3> centroid = elem_iter->CalculateCentroid();
        g_factor_sum = 0;

        for (AbstractTetrahedralMesh<3,3>::ElementIterator iter_2 = mrMesh.GetElementIteratorBegin();
             iter_2 != mrMesh.GetElementIteratorEnd();
             ++iter_2)
        {
            c_vector<double, 3> centroid_2 = iter_2->CalculateCentroid();
            double r = norm_2(centroid-centroid_2);
            if (r < r_max)
            {
                g_factor = exp(-r/(2*sigma*sigma));

                g_factor_sum += g_factor;

                mSmoothedStructureTensors[elem_iter->GetIndex()] += g_factor*mStructureTensors[iter_2->GetIndex()];
            }
        }

        mSmoothedStructureTensors[elem_iter->GetIndex()] /= g_factor_sum;
    }
}

#endif /*PAPILLARYFIBRECALCULATOR_HPP_*/

