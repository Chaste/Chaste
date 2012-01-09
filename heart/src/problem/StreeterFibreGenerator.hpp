/*

Copyright (C) University of Oxford, 2005-2012

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
#ifndef STREETERFIBREGENERATOR_HPP_
#define STREETERFIBREGENERATOR_HPP_

#include <vector>
#include <string>
#include <set>
#include "DistanceMapCalculator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "HeartGeometryInformation.hpp"

/**
 * Generate fibre in a ventricular mesh using the description in
 * Streeter DD, Jr, Spotnitz HM, Patel DP, Ross J, Jr, Sonnenblick EH.
 * Fiber orientation in the canine left ventricle during diastole and systole.
 * Circ Res. 1969 Mar;24(3):339–347.
 *
 * Formulae used to generate the fibre orientations are in
 * A Comparison of Monodomain and Bidomain Reaction-Diffusion Models for Action Potential Propagation in the Human Heart
 * Mark Potse, Bruno Dubé, Jacques Richer, Alain Vinet, and Ramesh M. Gulrajani
 * IEEE Trans. Biomed. Eng. 53(12):2425-2435, 2006.
 */
template<unsigned SPACE_DIM>
class StreeterFibreGenerator
{
private:
    AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh; /**< Reference to the mesh (used for calculating distances to epi and endo surfaces)*/

    HeartGeometryInformation<SPACE_DIM>* mpGeometryInfo; /**< Provides a method to calculate the relative position of a node with respect to two (or three) given surfaces*/

    c_vector <double, SPACE_DIM> mApexToBase; /**< Normalised direction from apex to base */

    /**
     * Compute the wallthickness of a given node based on a
     * neighbourhood average of its thickness and of those in the forward star.
     *
     * @param nodeIndex  The index of the node in question
     * @param wallThickness  vector of thickness of all nodes in node index order
     * @return Neighbourhood average thickness (will return 0 if the node is not local to this process)
     */
    double GetAveragedThicknessLocalNode(const unsigned nodeIndex, const std::vector<double>& wallThickness) const;

    /**
     * R is the maximum angle between the fibre and the v axis (heart region dependant)
     * @param  nodesRegionsForElement is a small vector containing the region tags of the element's nodes
     * @return  Pi/4 (if the element is in RV), Pi/3 otherwise
     */
   double GetFibreMaxAngle(const c_vector<HeartRegionType, SPACE_DIM+1>& nodesRegionsForElement) const;

public:
    /**
     * Constructor
     *
     * @param  rMesh reference to the tetrahedral mesh of the ventricles
     */
    StreeterFibreGenerator(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh);

    /**
     * Destructor
     */
    ~StreeterFibreGenerator();

    /**
     * Uses the names of files defining the diferent surfaces of the mesh to construct the geometry information class
     * File format: list of triangles
     *
     * @param rEpicardiumFile Epicardium surface
     * @param rRightVentricleFile Right Ventricle surface.
     * @param rLeftVentricleFile Left Ventricle surface
     * @param indexFromZero  Are the nodes in the original mesh file/surfaces indexed from 0?
     *
     * If either rRightVentricleFile or rLeftVentricleFile are the empty string, then it is assumed that this is a
     * wedge preparation for left or right ventricle, respectively.  That is, the ventricle with a non-empty string.
     * If both are empty strings then throws exception in HeartGeometryInformation.
     */
    void SetSurfaceFiles(const std::string &rEpicardiumFile,
                         const std::string &rRightVentricleFile,
                         const std::string &rLeftVentricleFile,
                         bool indexFromZero);

    /**
     * Generates an orthotropic fibre orientation model of the ventricular mesh provided. Assumes that the base-apex axis is x. Based on Streeter 1969 and Potse 2006
     *
     * File format: The first line indicates the number of elements. Each of the following lines contain SPACE_DIM vectors of SPACE_DIM elements for the
     * direction of the myofibre, the transverse to it in the plane of the myocite laminae and the normal to this laminae.
     *
     * @param outputDirectory Output directory
     * @param fibreOrientationFile Output file
     * @param logInfo Tells the method to output extra debug info. To be eliminated once it's fully tested
     *
     */
    void GenerateOrthotropicFibreOrientation(std::string outputDirectory, std::string fibreOrientationFile, bool logInfo=false);

    /**
     * Set the direction from apex to base
     * @param apexToBase  is a non-zero vector.  It will be stored in normalised form
     */
    void SetApexToBase(const c_vector<double, SPACE_DIM>& apexToBase);

    /**
     * Set the direction from apex to base
     * @param axis  is the Cartesian axis from apex to base.
     */
    void SetApexToBase(unsigned axis);
};

#endif /*STREETERFIBREGENERATOR_HPP_*/
