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


#ifndef AIRWAY_BRANCH_HPP_
#define AIRWAY_BRANCH_HPP_

#include "TetrahedralMesh.hpp"

#include <list>

/**
 * Represents a single airway made up of multiple <1,3> elements
 */
class AirwayBranch
{
public:

    /**
     * Constructor.
     *
     * @param radiusOnEdge Specifies whether radii are specified on nodes or on elements
     */
    AirwayBranch(bool radiusOnEdge = false);

    /**
     * Adds an element to the branch list
     *
     * @param pElement A pointer to the element to add
     */
    void AddElement(Element<1,3>* pElement);

    /**
     * Returns the list of elements that make up this branch
     *
     * @return A list of pointers to elements that make up the branch
     */
    std::list<Element<1,3>* > GetElements();

    /**
     * @return The true length of the branch (obtained by traversing the constitutive elements)
     */
    double GetLength();

    /**
     * @return the average radius of this branch
     */
    double GetAverageRadius();

    /**
     * @return the Poiseuille resistance, up to a constant, of this branch
     */
    double GetPoiseuilleResistance();

    /**
     * @return The direction of the branch (from the start node to the end node)
     */
    c_vector<double, 3> GetDirection();

    /**
     * @return true if this is a major branch
     */
    bool IsMajor();

    /**
     * @return The angle of this branch with respect to its parent
     */
    double GetBranchAngle();

    /**
     * @return The rotation angle of this branch
     */
    double GetRotationAngle();

    /**
     * @return Vector of all children of this branch
     */
    std::vector<AirwayBranch*> GetAllChildren();

    /**
     * @param pChild A pointer to a child to be added to mAllChildren
     */
    void AddChild(AirwayBranch* pChild);

    /**
     * @return The first child branch of this branch
     */
    AirwayBranch* GetChildOne();

    /**
     * @return The second child branch of this branch
     */
    AirwayBranch* GetChildTwo();

    /**
     * @param pChildOne A pointer to the first child of the branch
     */
    void SetChildOne(AirwayBranch* pChildOne);

    /**
     * @param pChildTwo A pointer to the second child of the branch
     */
    void SetChildTwo(AirwayBranch* pChildTwo);

    /**
     * @return The sibling branch of this branch
     */
    AirwayBranch* GetSibling();

    /**
     * @return The parent branch of this branch
     */
    AirwayBranch* GetParent();

    /**
     * @param pSibling A pointer to the sibling of the branch
     */
    void SetSibling(AirwayBranch* pSibling);

    /**
     * @param pParent A pointer to the parent of the branch
     */
    void SetParent(AirwayBranch* pParent);

    /**
     * @return unsigned index of branch
     */
    unsigned GetIndex();

    /**
     * @param index  New index of branch
     */
    void SetIndex(unsigned index);

    /**
     * Finds proximal node of branch.  Assumes only that elements are ordered
     * Proximal to Distal, which should be assured when branches are set up.
     *
     * @return pointer to proximal node of branch
     */
    Node<3>* GetProximalNode();

    /**
     * Finds distal node of branch.  Assumes only that elements are ordered
     * Proximal to Distal, which should be assured when branches are set up.
     *
     * @return pointer to proximal node of branch
     */
    Node<3>* GetDistalNode();

    /**
     * Calculate branch volume.
     * Note: this is based just on radius at each node; not well defined locally near branch-points.
     *
     * @return total volume of this branch
     */
    double GetBranchVolume();

    /**
     * Calculate lateral surface area of this branch.
     * Note: this is based just on radius at each node; not well defined locally near branch-points.
     *
     * @return lateral surface area this branch
     */
    double GetBranchLateralSurfaceArea();

    /**
     * Calculate centroid of this branch.
     * Note: this is based just on radius at each node; not well defined locally near branch-points.
     *
     * @return centroid of this branch
     */
    c_vector<double, 3> GetBranchCentroid();

    /**
     * @return True if this is a terminal airway branch
     */
    bool IsTerminal();

private:
    /**
     * A list of elements that make up this branch
     * These are assumed to be in order, from start element to end element!
     */
    std::list< Element<1,3>* > mElements;

    /**
     * A vector of all children of current branch
     * There will be more than two if trifurcations are present in the tree.
     */
    std::vector<AirwayBranch*> mAllChildren;

    /** The first child of this branch */
    AirwayBranch* mpChildOne;

    /** The second child of this branch */
    AirwayBranch* mpChildTwo;

    /** The parent of this branch */
    AirwayBranch* mpParent;

    /** The sibling of this branch */
    AirwayBranch* mpSibling;

    /** Branch index.  Currently optional and set in AirwayPropertiesCalculator. \todo: add index into constructors */
    unsigned mIndex;

    /** Flag to indicate whether airway radii are specified on nodes or edges */
    bool mRadiusOnEdge;
};

#endif // AIRWAY_BRANCH
