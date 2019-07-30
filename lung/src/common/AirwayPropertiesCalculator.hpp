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


#ifndef AIRWAY_PROPERTIES_CALCULATOR_HPP_
#define AIRWAY_PROPERTIES_CALCULATOR_HPP_

#include "TetrahedralMesh.hpp"
#include "AirwayTreeWalker.hpp"
#include "AirwayBranch.hpp"

/**
 * Calculates a number of morphological properties for an airway tree
 *
 * Properties include adding the airway generation as an attribute to the tree, calculating branching and rotation angles
 * and length/diameter ratios.
 */
class AirwayPropertiesCalculator
{
public:
    /**
     * Constructor
     *
     * @param rAirwaysMesh A mesh containing the airway tree to analyze
     * @param rootIndex The index of the inlet node (trachea entrance)
     * @param radiusOnEdge Whether airway radii are defined element wise or nodally.
     */
    AirwayPropertiesCalculator(TetrahedralMesh<1,3>& rAirwaysMesh,
                               unsigned rootIndex = 0u,
                               bool radiusOnEdge = false);

    /**
     * Destructor
     */
    ~AirwayPropertiesCalculator();


    /**
     * Get branch Generation
     *
     * @param pBranch The branch to get the Generation for
     * @return The branch's Generation order
     */
    unsigned GetBranchGeneration(AirwayBranch* pBranch);

    /**
     * Get branch Horsfield Order
     *
     * @param pBranch The branch to get the Horsfield order for
     * @return The branch's Horsfield order
     */
    unsigned GetBranchHorsfieldOrder(AirwayBranch* pBranch);

    /**
     * Get branch Strahler Order
     *
     * @param pBranch The branch to get the Strahler order for
     * @return The branch's Strahler order
     */
    unsigned GetBranchStrahlerOrder(AirwayBranch* pBranch);

    /** @return */
    double GetLengthOneOverLengthTwoMean() const;

    /** @return */
    double GetLengthOverDiameterMajorChildMean() const;

    /** @return */
    double GetLengthOverDiameterMean() const;

    /** @return */
    double GetLengthOverDiameterMinorChildMean() const;

    /** @return */
    double GetLengthOverLengthParentMean() const;

    /** @return */
    double GetDiameterOverParentDiameterMean() const;

    /** @return */
    double GetMajorDiameterOverParentDiameterMean() const;

    /** @return */
    double GetMinorDiameterOverMajorDiameterMean() const;

    /** @return */
    double GetMinorDiameterOverParentDiameterMean() const;

    /** @return */
    double GetPercentageLengthOverParentLengthLessThanOne() const;

    /** @return */
    double GetPhiMean() const;

    /** @return */
    double GetThetaMajorBranches() const;

    /** @return */
    double GetThetaMean() const;

    /** @return */
    double GetThetaMinorBranches() const;

    /** @return */
    double GetThetaParentDiameter2mmTo1mm() const;

    /** @return */
    double GetThetaParentDiameter3mmTo2mm() const;

    /** @return */
    double GetThetaParentDiameter4mmTo3mm() const;

    /** @return */
    double GetThetaParentDiameterGreaterThan4mm() const;

    /** @return A vector of branches in the airway tree */
    std::vector<AirwayBranch*> GetBranches();

    /** @return A vector of total subtree lengths, defined for each branch of the tree */
    std::vector<double> GetSubtreeBranchLengths();

    /** @return A vector of total subtree volumes, defined for each branch of the tree */
    std::vector<double> GetSubtreeBranchVolumes();

    /** @return A vector of total subtree lateral surface areas, defined for each branch of the tree */
    std::vector<double> GetSubtreeBranchLateralSurfaceAreas();

    /** @return A vector of total subtree Poiseuille resistances, defined for each branch of the tree */
    std::vector<double> GetSubtreePoiseuilleResistances();

    /** @return A vector of total subtree centroids, defined for each branch of the tree */
    std::vector<c_vector<double, 3> > GetSubtreeCentroids();

    /** @return A vector of upstream lengths, defined for each branch of the tree */
    std::vector<double> GetUpstreamBranchLengths();

    /** @return A vector of upstream volumes, defined for each branch of the tree */
    std::vector<double> GetUpstreamBranchVolumes();

    /** @return A vector of upstream lateral surface areas, defined for each branch of the tree */
    std::vector<double> GetUpstreamBranchLateralSurfaceAreas();

    /** @return A vector of upstream Poiseuille resistances, defined for each branch of the tree */
    std::vector<double> GetUpstreamPoiseuilleResistances();

    /** @return The maximum generation number of the terminal branches */
    unsigned GetMaximumTerminalGeneration();

    /** @return The minimum generation number of the terminal branches */
    unsigned GetMinimumTerminalGeneration();

    /** @return The mean generation number of the terminal branches */
    unsigned GetMeanTerminalGeneration();


    /**
     * Calculates the properties of all branches.
     *
     * Must be called prior to calling any of the GetPropertyX methods.
     */
    void CalculateBranchProperties();

    /**
     * Calculates properties of all subtrees.
     *
     * Must be called prior to calling any of the GetSubtreeX methods.
     */
    void CalculateSubtreeProperties();

    /**
     * Calculates properties of upstream paths.
     *
     * Must be called prior to calling any of the GetUpstreamX methods.
     */
    void CalculateUpstreamProperties();

private:

    /** A mesh containing the airways tree.  */
    TetrahedralMesh<1,3>& mAirwaysMesh;

    unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */

    /** Allows easy traversal of the airway tree */
    AirwayTreeWalker mWalker;

    /** Flag indicating whether airway radii are defined nodally or on elements */
    bool mRadiusOnEdge;

    /** An easy access list of airway branches */
    std::vector<AirwayBranch*> mBranches;

    /** The average branch angle of the tree */
    double mThetaMean;

    /** The branch angle spread */
    double mThetaSpread;

    /** The branch angle for branches with a parent diameter greater than 4mm */
    double mThetaParentDiameterGreaterThan4mm;

    /** The branch angle for branches with a parent diameter greater between 4mm and 3mm */
    double mThetaParentDiameter4mmTo3mm;

    /** The branch angle for branches with a parent diameter greater between 3mm and 2mm */
    double mThetaParentDiameter3mmTo2mm;

    /** The branch angle for branches with a parent diameter greater between 2mm and 1mm */
    double mThetaParentDiameter2mmTo1mm;

    /** The branch angle for major branches */
    double mThetaMajorBranches;

    /** The branch angle for major branches spread */
    double mThetaMajorBranchesSpread;

    /** The branch angle for minor branches */
    double mThetaMinorBranches;

    /** The branch angle for minor branches spread */
    double mThetaMinorBranchesSpread;

    /** The average rotation angle */
    double mPhiMean;

    /** The rotation angle spread */
    double mPhiSpread;

    /** The average length over diameter */
    double mLengthOverDiameterMean;

    /** The length over diameter spread */
    double mLengthOverDiameterSpread;

    /** The length over diameter for minor branches mean value */
    double mLengthOverDiameterMinorChildMean;

    /** The length over diameter for minor branches spread */
    double mLengthOverDiameterMinorChildSpread;

    /** The length over diameter for major branches mean value */
    double mLengthOverDiameterMajorChildMean;

    /** The length over diameter for major branches spread */
    double mLengthOverDiameterMajorChildSpread;

    /** The diameter over the parent diameter mean */
    double mDiameterOverParentDiameterMean;

    /** The minor diameter over the major diameter mean */
    double mMinorDiameterOverMajorDiameterMean;

    /** The minor diameter over the major diameter spread */
    double mMinorDiameterOverMajorDiameterSpread;

    /** The minor diameter over the parent diameter mean */
    double mMinorDiameterOverParentDiameterMean;

    /** The minor diameter over the parent diameter spread */
    double mMinorDiameterOverParentDiameterSpread;

    /** The major diameter over the parent diameter mean */
    double mMajorDiameterOverParentDiameterMean;

    /** The major diameter over the parent diameter spread */
    double mMajorDiameterOverParentDiameterSpread;

    /** The length over the parent length mean */
    double mLengthOverLengthParentMean;

    /** The length over the parent length spread */
    double mLengthOverLengthParentSpread;

    /** The percentage of branches that are shorter than their parent */
    double mPercentageLengthOverParentLengthLessThanOne;

    /** The ratio of sibling lengths mean */
    double mLengthOneOverLengthTwoMean;

    /** The ratio of sibling lengths spread */
    double mLengthOneOverLengthTwoSpread;

    /** For each branch, the total length of all branches in the distal direction */
    std::vector<double> mTotalSubtreeBranchLength;

    /** For each branch, the total volume of all branches in the distal direction */
    std::vector<double> mTotalSubtreeBranchVolume;

    /** For each branch, the total lateral surface area of all branches in the distal direction */
    std::vector<double> mTotalSubtreeBranchLateralSurfaceArea;

    /** For each branch, the Poiseuille resistance of all branches in the distal direction */
    std::vector<double> mTotalSubtreePoiseuilleResistance;

    /** For each branch, the centroid (by volume) all branches in the distal direction */
    std::vector<c_vector<double, 3> > mTotalSubtreeCentroid;

    /** For each branch, the total length of all branches in the proximal (tracheal) direction */
    std::vector<double> mUpstreamPathBranchLengths;

    /** For each branch, the total volume of all branches in the proximal (tracheal) direction */
    std::vector<double> mUpstreamPathBranchVolumes;

    /** For each branch, the total lateral surface area of all branches in the proximal (tracheal) direction */
    std::vector<double> mUpstreamPathBranchLateralSurfaceAreas;

    /** For each branch, the Poiseuille resistance of all branches in the proximal (tracheal) direction */
    std::vector<double> mUpstreamPathPoiseuilleResistances;


    /**
     * Recursively setup multiple element branches to allow properties to be calculated easily
     *
     * @param pElement The current element (should not have been added to the branch yet)
     * @param pBranch The current branch
     */
    void SetupBranches(Element<1,3>* pElement, AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate subtree length.
     * Helper method for CalculateSubtreeProperties().
     *
     * @param pBranch The current branch
     * @return The subtree length
     */
    double RecursivelyCalculateSubtreeLength(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate subtree volume.
     * Helper method for CalculateSubtreeProperties().
     *
     * @param pBranch The current branch
     * @return The subtree volume
     */
    double RecursivelyCalculateSubtreeVolume(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate subtree lateral surface area.
     * Helper method for CalculateSubtreeProperties().
     *
     * @param pBranch The current branch
     * @return The subtree lateral surface area
     */
    double RecursivelyCalculateSubtreeLateralSurfaceArea(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate subtree Poiseuille resistance.
     * Helper method for CalculateSubtreeProperties().
     *
     * @param pBranch The current branch
     * @return The subtree Poiseuille resistance
     */
    double RecursivelyCalculateSubtreePoiseuilleResistance(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate subtree centroid.
     * Helper method for CalculateSubtreeProperties().
     *
     * @param pBranch The current branch
     * @return The subtree centroid
     */
    c_vector<double, 3> RecursivelyCalculateSubtreeCentroid(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate upstream lengths.
     * Helper method for CalculateUpstreamProperties().
     *
     * @param pBranch The current branch
     */
    void RecursivelyCalculateUpstreamLengths(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate upstream volumes.
     * Helper method for CalculateUpstreamProperties().
     *
     * @param pBranch The current branch
     */
    void RecursivelyCalculateUpstreamVolumes(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate upstream lateral surface areas.
     * Helper method for CalculateUpstreamProperties().
     *
     * @param pBranch The current branch
     */
    void RecursivelyCalculateUpstreamLateralSurfaceAreas(AirwayBranch* pBranch);

    /**
     * Recursively walk airway tree to calculate upstream Poiseuille resistances.
     * Helper method for CalculateUpstreamProperties().
     *
     * @param pBranch The current branch
     */
    void RecursivelyCalculateUpstreamPoiseuilleResistances(AirwayBranch* pBranch);
};

#endif // AIRWAY_PROPERTIES_CALCULATOR
