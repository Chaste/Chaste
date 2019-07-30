/*

Copyright (c) 2005-2019, University of Oxford.
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
#ifndef STREETERFIBREGENERATOR_HPP_
#define STREETERFIBREGENERATOR_HPP_

#include <vector>
#include <string>
#include <set>
#include "DistanceMapCalculator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "HeartGeometryInformation.hpp"
#include "AbstractPerElementWriter.hpp"
#include "OutputFileHandler.hpp"

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
 *
 * Output file format: The first line indicates the number of elements.
 * Each of the following lines contain SPACE_DIM vectors of SPACE_DIM elements for the
 * direction of the myofibre, the transverse to it in the plane of
 * the myocyte laminae and the normal to this laminae.
 */
template<unsigned SPACE_DIM>
class StreeterFibreGenerator : public AbstractPerElementWriter<SPACE_DIM, SPACE_DIM, SPACE_DIM*SPACE_DIM>
{
private:
    HeartGeometryInformation<SPACE_DIM>* mpGeometryInfo; /**< Provides a method to calculate the relative position of a node with respect to two (or three) given surfaces*/

    c_vector <double, SPACE_DIM> mApexToBase; /**< Normalised direction from apex to base */

    /**
     * Compute the wall thickness of a given node based on a
     * neighbourhood average of its thickness and of those in the forward star.
     *
     * @param nodeIndex  The index of the node in question
     * @param wallThickness  vector of thickness of all nodes in node index order
     * @return Neighbourhood average thickness (will return 0 if the node is not local to this process)
     */
    double GetAveragedThicknessLocalNode(const unsigned nodeIndex,
                                         const std::vector<double>& wallThickness) const;

    /**
     * R is the maximum angle between the fibre and the v axis (heart region dependant)
     * @param  nodesRegionsForElement is a small vector containing the region tags of the element's nodes
     * @return  Pi/4 (if the element is in RV), Pi/3 otherwise
     */
   double GetFibreMaxAngle(const c_vector<HeartRegionType, SPACE_DIM+1>& nodesRegionsForElement) const;

   /** Wall thickness at each node in the mesh. */
   std::vector<double> mWallThickness;

   /** Wall thickness at each node, smoothed by averaging over local nodes by #GetAveragedThicknessLocalNode().*/
   std::vector<double> mAveragedWallThickness;

   /** Whether to write Streeter generation log files for regions and wall thicknesses */
   bool mLogInfo;

protected:

   /**
    * Associate an element with a fibre direction.
    *
    * This is only called (by abstract class) on processes which own pElement.
    *
    * @param pElement  a locally-owned element for which to calculate or lookup some data
    * @param localElementIndex  the index of pElement in the local vector.
    * @param rData  the double-precision data to write to file (output from the method). Gives Fibre x 3 (space dim), sheet(3), normal(3) in one vector.
    */
   void Visit(Element<SPACE_DIM, SPACE_DIM>* pElement,
                      unsigned localElementIndex,
                      c_vector<double, SPACE_DIM*SPACE_DIM>& rData);

   /**
    * Overridden method to write the header line.
    */
   void WriteHeaderOnMaster();

   /**
    * Does calculations to generate an orthotropic fibre orientation model of the ventricular mesh provided.
    * In particular - populates the member variables #mWallThickness and #mAveragedWallThickness.
    *
    * Also writes these to file if SetLogInfo() has been called.
    *
    * Assumes that the base-apex axis is x. Based on Streeter 1969 and Potse 2006
    *
    * @param rOutputDirectory Handler for output directory
    */
   void PreWriteCalculations(OutputFileHandler& rOutputDirectory);

public:
    /**
     * Constructor
     *
     * @param rMesh  Reference to the tetrahedral mesh of the ventricles.
     */
    StreeterFibreGenerator(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh);

    /**
     * Destructor
     */
    ~StreeterFibreGenerator();

    /**
     * Uses the names of files defining the different surfaces of the mesh to construct the geometry information class
     * File formats: list of nodes, either one per line or multiple (e.g. nodes in each boundary element on surface).
     *
     * @param rEpicardiumFile  Epicardium surface nodes (global indices).
     * @param rRightVentricleFile  Right Ventricle surface nodes (global indices).
     * @param rLeftVentricleFile  Left Ventricle surface nodes (global indices).
     * @param indexFromZero  Are the nodes in the original mesh file/surface files indexed from 0?
     *
     * If either rRightVentricleFile or rLeftVentricleFile are the empty string, then it is assumed that this is a
     * wedge preparation for left or right ventricle, respectively.  That is, the ventricle with a non-empty string.
     * If both are empty strings then throws exception in HeartGeometryInformation.
     */
    void SetSurfaceFiles(const std::string& rEpicardiumFile,
                         const std::string& rRightVentricleFile,
                         const std::string& rLeftVentricleFile,
                         bool indexFromZero);

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

    /**
     * Tells the WriteData method to output extra debug info on nodes in particular regions,
     * and wall thicknesses.
     *
     * @param logInfo  whether or not to write log files.
     */
    void SetLogInfo(bool logInfo = true);
};

#endif /*STREETERFIBREGENERATOR_HPP_*/
