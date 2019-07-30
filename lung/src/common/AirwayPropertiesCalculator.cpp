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


#include "AirwayPropertiesCalculator.hpp"

AirwayPropertiesCalculator::AirwayPropertiesCalculator(TetrahedralMesh<1,3>& rAirwaysMesh,
                                                       unsigned rootIndex,
                                                       bool radiusOnEdge) :
                                                           mAirwaysMesh(rAirwaysMesh),
                                                           mOutletNodeIndex(rootIndex),
                                                           mWalker(mAirwaysMesh, mOutletNodeIndex),
                                                           mRadiusOnEdge(radiusOnEdge)
{
    // Get the head element & process
    Node<3>* p_node = mAirwaysMesh.GetNode(mOutletNodeIndex);
    Element<1,3>* p_element = mAirwaysMesh.GetElement(*(p_node->ContainingElementsBegin()));

    AirwayBranch* p_head_branch = new AirwayBranch(mRadiusOnEdge);

    // Set index to zero for head branch
    p_head_branch->SetIndex(0);

    mBranches.push_back(p_head_branch);
    SetupBranches(p_element, p_head_branch);
}

AirwayPropertiesCalculator::~AirwayPropertiesCalculator()
{
    for (std::vector<AirwayBranch*>::iterator iter = mBranches.begin();
         iter != mBranches.end();
         ++iter)
    {
        delete (*iter);
    }
}

void AirwayPropertiesCalculator::SetupBranches(Element<1,3>* pElement, AirwayBranch* pBranch)
{
    pBranch->AddElement(pElement);

    if (mWalker.GetNumberOfChildElements(pElement) == 0u)
    {
        return;
    }
    else if (mWalker.GetNumberOfChildElements(pElement) == 1u)
    {
        SetupBranches(mWalker.GetChildElements(pElement)[0], pBranch);
    }
    else // if (mWalker.GetNumberOfChildElements(pElement) >= 2u)
    {
        AirwayBranch* p_child_branch_one = new AirwayBranch(mRadiusOnEdge);
        AirwayBranch* p_child_branch_two = new AirwayBranch(mRadiusOnEdge);
        p_child_branch_one->SetSibling(p_child_branch_two);
        p_child_branch_two->SetSibling(p_child_branch_one);
        p_child_branch_one->SetParent(pBranch);
        p_child_branch_two->SetParent(pBranch);
        pBranch->AddChild(p_child_branch_one);
        pBranch->AddChild(p_child_branch_two);
        pBranch->SetChildOne(p_child_branch_one);
        pBranch->SetChildTwo(p_child_branch_two);

        // Hack to add correct indices to branch.
        ///\todo: add indexing to constructor
        mBranches.push_back(p_child_branch_one);
        p_child_branch_one->SetIndex(mBranches.size()-1);
        mBranches.push_back(p_child_branch_two);
        p_child_branch_two->SetIndex(mBranches.size()-1);

        SetupBranches(mWalker.GetChildElements(pElement)[0], p_child_branch_one);
        SetupBranches(mWalker.GetChildElements(pElement)[1], p_child_branch_two);
    }

    if (mWalker.GetNumberOfChildElements(pElement) >= 3u)
    {
        // std::cout << "Warning: Trifurcation detected. Third (and higher) branch will be ignored when calculating statistics. " << std::endl;

        for (unsigned extra_branch = 2; extra_branch < mWalker.GetNumberOfChildElements(pElement); ++extra_branch)
        {
            AirwayBranch* p_child_branch = new AirwayBranch(mRadiusOnEdge);
            p_child_branch->SetParent(pBranch);
            pBranch->AddChild(p_child_branch);

            // Hack to add unique index to branch
            ///\todo: add indexing to constructor
            mBranches.push_back(p_child_branch);
            p_child_branch->SetIndex(mBranches.size()-1);

            SetupBranches(mWalker.GetChildElements(pElement)[extra_branch], p_child_branch);
        }
    }
}

double AirwayPropertiesCalculator::GetLengthOneOverLengthTwoMean() const
{
    return mLengthOneOverLengthTwoMean;
}

double AirwayPropertiesCalculator::GetLengthOverDiameterMajorChildMean() const
{
    return mLengthOverDiameterMajorChildMean;
}

double AirwayPropertiesCalculator::GetLengthOverDiameterMean() const
{
    return mLengthOverDiameterMean;
}

double AirwayPropertiesCalculator::GetLengthOverDiameterMinorChildMean() const
{
    return mLengthOverDiameterMinorChildMean;
}

double AirwayPropertiesCalculator::GetLengthOverLengthParentMean() const
{
    return mLengthOverLengthParentMean;
}

double AirwayPropertiesCalculator::GetDiameterOverParentDiameterMean() const
{
    return mDiameterOverParentDiameterMean;
}

double AirwayPropertiesCalculator::GetMajorDiameterOverParentDiameterMean() const
{
    return mMajorDiameterOverParentDiameterMean;
}

double AirwayPropertiesCalculator::GetMinorDiameterOverMajorDiameterMean() const
{
    return mMinorDiameterOverMajorDiameterMean;
}

double AirwayPropertiesCalculator::GetMinorDiameterOverParentDiameterMean() const
{
    return mMinorDiameterOverParentDiameterMean;
}

double AirwayPropertiesCalculator::GetPercentageLengthOverParentLengthLessThanOne() const
{
    return mPercentageLengthOverParentLengthLessThanOne;
}

double AirwayPropertiesCalculator::GetPhiMean() const
{
    return mPhiMean;
}

double AirwayPropertiesCalculator::GetThetaMajorBranches() const
{
    return mThetaMajorBranches;
}

double AirwayPropertiesCalculator::GetThetaMean() const
{
    return mThetaMean;
}

double AirwayPropertiesCalculator::GetThetaMinorBranches() const
{
    return mThetaMinorBranches;
}

double AirwayPropertiesCalculator::GetThetaParentDiameter2mmTo1mm() const
{
    return mThetaParentDiameter2mmTo1mm;
}

double AirwayPropertiesCalculator::GetThetaParentDiameter3mmTo2mm() const
{
    return mThetaParentDiameter3mmTo2mm;
}

double AirwayPropertiesCalculator::GetThetaParentDiameter4mmTo3mm() const
{
    return mThetaParentDiameter4mmTo3mm;
}

double AirwayPropertiesCalculator::GetThetaParentDiameterGreaterThan4mm() const
{
    return mThetaParentDiameterGreaterThan4mm;
}

std::vector<AirwayBranch*> AirwayPropertiesCalculator::GetBranches()
{
    return mBranches;
}

unsigned AirwayPropertiesCalculator::GetBranchGeneration(AirwayBranch* pBranch)
{
    return mWalker.GetElementGeneration(pBranch->GetElements().front());
}

unsigned AirwayPropertiesCalculator::GetBranchHorsfieldOrder(AirwayBranch* pBranch)
{
    return mWalker.GetElementHorsfieldOrder(pBranch->GetElements().front());
}

unsigned AirwayPropertiesCalculator::GetBranchStrahlerOrder(AirwayBranch* pBranch)
{
    return mWalker.GetElementStrahlerOrder(pBranch->GetElements().front());
}

std::vector<double> AirwayPropertiesCalculator::GetSubtreeBranchLengths()
{
    return mTotalSubtreeBranchLength;
}

std::vector<double> AirwayPropertiesCalculator::GetSubtreeBranchVolumes()
{
    return mTotalSubtreeBranchVolume;
}

std::vector<double> AirwayPropertiesCalculator::GetSubtreeBranchLateralSurfaceAreas()
{
    return mTotalSubtreeBranchLateralSurfaceArea;
}

std::vector<double> AirwayPropertiesCalculator::GetSubtreePoiseuilleResistances()
{
    return mTotalSubtreePoiseuilleResistance;
}

std::vector<c_vector<double, 3> > AirwayPropertiesCalculator::GetSubtreeCentroids()
{
    return mTotalSubtreeCentroid;
}

std::vector<double> AirwayPropertiesCalculator::GetUpstreamBranchLengths()
{
    return mUpstreamPathBranchLengths;
}

std::vector<double> AirwayPropertiesCalculator::GetUpstreamBranchVolumes()
{
    return mUpstreamPathBranchVolumes;
}

std::vector<double> AirwayPropertiesCalculator::GetUpstreamBranchLateralSurfaceAreas()
{
    return mUpstreamPathBranchLateralSurfaceAreas;
}

std::vector<double> AirwayPropertiesCalculator::GetUpstreamPoiseuilleResistances()
{
    return mUpstreamPathPoiseuilleResistances;
}

unsigned AirwayPropertiesCalculator::GetMaximumTerminalGeneration()
{
    unsigned max_generation = 0;

    for (std::vector<AirwayBranch*>::iterator iter = mBranches.begin();
         iter != mBranches.end();
         ++iter)
    {
        if ((*iter)->IsTerminal())
        {
            unsigned generation = mWalker.GetElementGeneration((*iter)->GetElements().front());
            if (generation > max_generation)
            {
                max_generation = generation;
            }
        }
    }

    return max_generation;
}

unsigned AirwayPropertiesCalculator::GetMinimumTerminalGeneration()
{
    unsigned min_generation = std::numeric_limits<unsigned>::max();

    for (std::vector<AirwayBranch*>::iterator iter = mBranches.begin();
         iter != mBranches.end();
         ++iter)
    {
        if ((*iter)->IsTerminal())
        {
            unsigned generation = mWalker.GetElementGeneration((*iter)->GetElements().front());
            if (generation < min_generation)
            {
                min_generation = generation;
            }
        }
    }

    return min_generation;
}

unsigned AirwayPropertiesCalculator::GetMeanTerminalGeneration()
{
    unsigned mean_generation = 0;
    unsigned terminal_branches_count = 0;

    for (std::vector<AirwayBranch*>::iterator iter = mBranches.begin();
         iter != mBranches.end();
         ++iter)
    {
        if ((*iter)->IsTerminal())
        {
            unsigned generation = mWalker.GetElementGeneration((*iter)->GetElements().front());

            mean_generation += generation;
            terminal_branches_count++;
        }
    }

    mean_generation /= terminal_branches_count;
    return mean_generation;
}

void AirwayPropertiesCalculator::CalculateBranchProperties()
{
    mThetaMean = 0.0;
    mThetaSpread = 0.0;
    mThetaParentDiameterGreaterThan4mm = 0.0;
    unsigned thetaParentDiameterGreaterThan4mmCount = 0u;
    mThetaParentDiameter4mmTo3mm = 0.0;
    unsigned thetaParentDiameter4mmTo3mmCount = 0u;
    mThetaParentDiameter3mmTo2mm = 0.0;
    unsigned thetaParentDiameter3mmTo2mmCount = 0u;
    mThetaParentDiameter2mmTo1mm = 0.0;
    unsigned thetaParentDiameter2mmTo1mmCount = 0u;
    mThetaMajorBranches = 0.0;
    mThetaMajorBranchesSpread = 0.0;
    unsigned majorBranchesCount = 0u;
    mThetaMinorBranches = 0.0;
    mThetaMinorBranchesSpread = 0.0;
    unsigned minorBranchesCount = 0u;
    mPhiMean = 0.0;
    mPhiSpread = 0.0;
    unsigned phiMeanCount = 0u;
    mLengthOverDiameterMean = 0.0;
    mLengthOverDiameterSpread = 0.0;
    mLengthOverDiameterMinorChildMean = 0.0;
    mLengthOverDiameterMinorChildSpread = 0.0;
    mLengthOverDiameterMajorChildMean = 0.0;
    mLengthOverDiameterMajorChildSpread = 0.0;
    mDiameterOverParentDiameterMean = 0.0;
    mMinorDiameterOverMajorDiameterMean = 0.0;
    mMinorDiameterOverMajorDiameterSpread = 0.0;
    mMinorDiameterOverParentDiameterMean = 0.0;
    mMinorDiameterOverParentDiameterSpread = 0.0;
    mMajorDiameterOverParentDiameterMean = 0.0;
    mMajorDiameterOverParentDiameterSpread = 0.0;
    mLengthOverLengthParentMean = 0.0;
    mLengthOverLengthParentSpread = 0.0;
    mPercentageLengthOverParentLengthLessThanOne = 0.0;
    unsigned lengthOverParentLengthLessThanOneCount = 0u;
    mLengthOneOverLengthTwoMean = 0.0;
    mLengthOneOverLengthTwoSpread = 0.0;
    unsigned lengthOneOverLengthTwoCount = 0u;

    for (std::vector<AirwayBranch*>::iterator iter = mBranches.begin();
         iter != mBranches.end();
         ++iter)
    {
        bool is_major = (*iter)->IsMajor();
        double length = (*iter)->GetLength();
        double diameter = 2*(*iter)->GetAverageRadius();

        assert(diameter > 0.0); //The airway tree must have a well defined set of radii to calculate branch properties

        mLengthOverDiameterMean += length/diameter;

        if ((*iter)->GetParent() != nullptr)
        {
            double theta = (*iter)->GetBranchAngle();
            double parent_diameter = 2*(*iter)->GetParent()->GetAverageRadius();

            mLengthOverLengthParentMean += length/(*iter)->GetParent()->GetLength();
            mDiameterOverParentDiameterMean += diameter/parent_diameter;

            mThetaMean += theta;

            if (parent_diameter > 4.0)
            {
                mThetaParentDiameterGreaterThan4mm += theta;
                thetaParentDiameterGreaterThan4mmCount++;
            }
            else if (parent_diameter < 4.0 && parent_diameter > 3.0)
            {
                mThetaParentDiameter4mmTo3mm += theta;
                thetaParentDiameter4mmTo3mmCount++;
            }
            else if (parent_diameter < 3.0 && parent_diameter > 2.0)
            {
                mThetaParentDiameter3mmTo2mm += theta;
                thetaParentDiameter3mmTo2mmCount++;
            }
            else if (parent_diameter < 2.0 && parent_diameter > 1.0)
            {
                mThetaParentDiameter2mmTo1mm += theta;
                thetaParentDiameter2mmTo1mmCount++;
            }

            if ((*iter)->GetSibling() != nullptr)
            {
                if (is_major)
                {
                    mThetaMajorBranches += theta;
                    mLengthOverDiameterMajorChildMean += length/diameter;
                    mMajorDiameterOverParentDiameterMean += diameter/parent_diameter;
                    majorBranchesCount++;
                }
                else
                {
                    mThetaMinorBranches += theta;
                    mLengthOverDiameterMinorChildMean += length/diameter;
                    mMinorDiameterOverMajorDiameterMean += diameter/(2*(*iter)->GetSibling()->GetAverageRadius());
                    mMinorDiameterOverParentDiameterMean += diameter/parent_diameter;
                    minorBranchesCount++;
                }

                if (length < (*iter)->GetSibling()->GetLength())
                {
                    mLengthOneOverLengthTwoMean += length/(*iter)->GetSibling()->GetLength();
                    lengthOneOverLengthTwoCount++;
                }
            }

            if ((*iter)->GetParent() != nullptr && length < (*iter)->GetParent()->GetLength())
            {
                lengthOverParentLengthLessThanOneCount++;
            }
        }

        if ((*iter)->GetParent() != nullptr && (*iter)->GetParent()->GetSibling() != nullptr && (*iter)->GetSibling() != nullptr)
        {
            mPhiMean += (*iter)->GetRotationAngle();
            phiMeanCount++;
        }
    }

    mThetaMean /= (mBranches.size() - 1);

    if (thetaParentDiameterGreaterThan4mmCount > 0u)
    {
        mThetaParentDiameterGreaterThan4mm /= thetaParentDiameterGreaterThan4mmCount;
    }

    if (thetaParentDiameter4mmTo3mmCount > 0u)
    {
        mThetaParentDiameter4mmTo3mm /= thetaParentDiameter4mmTo3mmCount;
    }

    if (thetaParentDiameter3mmTo2mmCount > 0u)
    {
        mThetaParentDiameter3mmTo2mm /= thetaParentDiameter3mmTo2mmCount;
    }

    if (thetaParentDiameter2mmTo1mmCount > 0u)
    {
        mThetaParentDiameter2mmTo1mm /= thetaParentDiameter2mmTo1mmCount;
    }

    if (majorBranchesCount > 0u)
    {
        mThetaMajorBranches /= majorBranchesCount;
        mLengthOverDiameterMajorChildMean /= majorBranchesCount;
        mMajorDiameterOverParentDiameterMean /= majorBranchesCount;
    }

    if (minorBranchesCount > 0u)
    {
        mThetaMinorBranches /= minorBranchesCount;
        mLengthOverDiameterMinorChildMean /= minorBranchesCount;
        mMinorDiameterOverMajorDiameterMean /= minorBranchesCount;
        mMinorDiameterOverParentDiameterMean /= minorBranchesCount;
    }

    if (phiMeanCount > 0u)
    {
        mPhiMean /= phiMeanCount;
    }

    if (mBranches.size() > 0u)
    {
        mLengthOverDiameterMean /= mBranches.size();
    }

    if ((mBranches.size() - 1) > 0u)
    {
        mDiameterOverParentDiameterMean /= (mBranches.size() - 1);
        mLengthOverLengthParentMean /= (mBranches.size() - 1);
        mPercentageLengthOverParentLengthLessThanOne = ((double)lengthOverParentLengthLessThanOneCount)/((double)(mBranches.size() - 1));
    }

    if (lengthOneOverLengthTwoCount > 0u)
    {
        mLengthOneOverLengthTwoMean /= lengthOneOverLengthTwoCount;
    }
}

void AirwayPropertiesCalculator::CalculateSubtreeProperties()
{
    // Lengthen vectors to store a value for each branch
    unsigned num_branches = mBranches.size();
    mTotalSubtreeBranchLength.resize(num_branches);
    mTotalSubtreeBranchVolume.resize(num_branches);
    mTotalSubtreeBranchLateralSurfaceArea.resize(num_branches);
    mTotalSubtreePoiseuilleResistance.resize(num_branches);
    mTotalSubtreeCentroid.resize(num_branches);

    // Get pointer to trachea, and double check it's the correct
    AirwayBranch* trachea = mBranches[0];
    assert(trachea->GetIndex() == 0);

    // Begin the recursions, sanity checking results
    double total_tree_length = RecursivelyCalculateSubtreeLength(trachea);
    assert(fabs(total_tree_length - mTotalSubtreeBranchLength[trachea->GetIndex()]) < 1e-6);
    UNUSED_OPT(total_tree_length);

    double total_tree_volume = RecursivelyCalculateSubtreeVolume(trachea);
    assert(fabs(total_tree_volume - mTotalSubtreeBranchVolume[trachea->GetIndex()]) < 1e-6);
    UNUSED_OPT(total_tree_volume);

    double total_tree_lateral_surface_area = RecursivelyCalculateSubtreeLateralSurfaceArea(trachea);
    assert(fabs(total_tree_lateral_surface_area - mTotalSubtreeBranchLateralSurfaceArea[trachea->GetIndex()]) < 1e-6);
    UNUSED_OPT(total_tree_lateral_surface_area);

    // The variable assigned to is unused; does the method call have desirable side effects?
    double total_tree_poiseuille_resistance = RecursivelyCalculateSubtreePoiseuilleResistance(trachea);
    RecursivelyCalculateSubtreePoiseuilleResistance(trachea);
    assert(fabs(total_tree_poiseuille_resistance - mTotalSubtreePoiseuilleResistance[trachea->GetIndex()]) < 1e-6);
    UNUSED_OPT(total_tree_poiseuille_resistance);

    c_vector<double, 3> total_centroid = RecursivelyCalculateSubtreeCentroid(trachea);
    assert(norm_2(total_centroid - mTotalSubtreeCentroid[trachea->GetIndex()]) < 1e-6);
    UNUSED_OPT(total_centroid);

    // Correct centroid by undoing volume-weighting
    for (unsigned branch_idx = 0 ; branch_idx < num_branches ; branch_idx ++)
    {
        mTotalSubtreeCentroid[branch_idx] /= mTotalSubtreeBranchVolume[branch_idx];
    }
}

void AirwayPropertiesCalculator::CalculateUpstreamProperties()
{
    // Lengthen vectors to store a value for each branch
    unsigned num_branches = mBranches.size();
    mUpstreamPathBranchLengths.resize(num_branches);
    mUpstreamPathBranchVolumes.resize(num_branches);
    mUpstreamPathBranchLateralSurfaceAreas.resize(num_branches);
    mUpstreamPathPoiseuilleResistances.resize(num_branches);

    // Get pointer to trachea, and double check it's the correct
    AirwayBranch* trachea = mBranches[0];
    assert(trachea->GetIndex() == 0);

    // Begin the recursions
    RecursivelyCalculateUpstreamLengths(trachea);
    RecursivelyCalculateUpstreamVolumes(trachea);
    RecursivelyCalculateUpstreamLateralSurfaceAreas(trachea);
    RecursivelyCalculateUpstreamPoiseuilleResistances(trachea);
}

double AirwayPropertiesCalculator::RecursivelyCalculateSubtreeLength(AirwayBranch* pBranch)
{
    mTotalSubtreeBranchLength[pBranch->GetIndex()] = pBranch->GetLength();

    std::vector<AirwayBranch*> all_children = pBranch->GetAllChildren();

    // Loop over all children.  If there are no children, the loop will not be executed.
    for (unsigned child_branch_idx = 0 ; child_branch_idx < all_children.size() ; child_branch_idx++)
    {
        mTotalSubtreeBranchLength[pBranch->GetIndex()] += RecursivelyCalculateSubtreeLength(all_children[child_branch_idx]);
    }

    // Return the total below, including current branch length
    return mTotalSubtreeBranchLength[pBranch->GetIndex()];
}

double AirwayPropertiesCalculator::RecursivelyCalculateSubtreeVolume(AirwayBranch* pBranch)
{
    mTotalSubtreeBranchVolume[pBranch->GetIndex()] = pBranch->GetBranchVolume();

    std::vector<AirwayBranch*> all_children = pBranch->GetAllChildren();

    // Loop over all children.  If there are no children, the loop will not be executed.
    for (unsigned child_branch_idx = 0 ; child_branch_idx < all_children.size() ; child_branch_idx++)
    {
        mTotalSubtreeBranchVolume[pBranch->GetIndex()] += RecursivelyCalculateSubtreeVolume(all_children[child_branch_idx]);
    }

    // Return the total below, including current branch length
    return mTotalSubtreeBranchVolume[pBranch->GetIndex()];
}

double AirwayPropertiesCalculator::RecursivelyCalculateSubtreeLateralSurfaceArea(AirwayBranch* pBranch)
{
    mTotalSubtreeBranchLateralSurfaceArea[pBranch->GetIndex()] = pBranch->GetBranchLateralSurfaceArea();

    std::vector<AirwayBranch*> all_children = pBranch->GetAllChildren();

    // Loop over all children.  If there are no children, the loop will not be executed.
    for (unsigned child_branch_idx = 0 ; child_branch_idx < all_children.size() ; child_branch_idx++)
    {
        mTotalSubtreeBranchLateralSurfaceArea[pBranch->GetIndex()] += RecursivelyCalculateSubtreeLateralSurfaceArea(all_children[child_branch_idx]);
    }

    // Return the total below, including current branch surface area
    return mTotalSubtreeBranchLateralSurfaceArea[pBranch->GetIndex()];
}

double AirwayPropertiesCalculator::RecursivelyCalculateSubtreePoiseuilleResistance(AirwayBranch* pBranch)
{
    mTotalSubtreePoiseuilleResistance[pBranch->GetIndex()] = pBranch->GetPoiseuilleResistance();

    std::vector<AirwayBranch*> all_children = pBranch->GetAllChildren();

    // Only carry on with recursion if there are children
    if (!all_children.empty())
    {
        // If there are children, currently we assume there are exactly two
        assert(all_children.size() == 2);

        double child0_subtree_res = RecursivelyCalculateSubtreePoiseuilleResistance(all_children[0]);
        double child1_subtree_res = RecursivelyCalculateSubtreePoiseuilleResistance(all_children[1]);

        // Update the Poiseuille resistance by reciprocal sum
        mTotalSubtreePoiseuilleResistance[pBranch->GetIndex()] += ( child0_subtree_res * child1_subtree_res / (child0_subtree_res + child1_subtree_res) );
    }

    // Return the total below, including current branch resistance
    return mTotalSubtreePoiseuilleResistance[pBranch->GetIndex()];
}

c_vector<double, 3> AirwayPropertiesCalculator::RecursivelyCalculateSubtreeCentroid(AirwayBranch* pBranch)
{
    // Need to rescale centroid by volume of the element, then divide later by total subtree volume
    mTotalSubtreeCentroid[pBranch->GetIndex()] = pBranch->GetBranchCentroid() * pBranch->GetBranchVolume();

    std::vector<AirwayBranch*> all_children = pBranch->GetAllChildren();

    // Loop over all children.  If there are no children, the loop will not be executed.
    for (unsigned child_branch_idx = 0 ; child_branch_idx < all_children.size() ; child_branch_idx++)
    {
        mTotalSubtreeCentroid[pBranch->GetIndex()] += RecursivelyCalculateSubtreeCentroid(all_children[child_branch_idx]);
    }

    // Return the final value
    return mTotalSubtreeCentroid[pBranch->GetIndex()];
}

void AirwayPropertiesCalculator::RecursivelyCalculateUpstreamLengths(AirwayBranch* pBranch)
{
    unsigned current_index = pBranch->GetIndex();

    // If current branch is trachea, there is nothing upstream, so just need current branch data
    if (pBranch->GetParent() == nullptr)
    {
        mUpstreamPathBranchLengths[current_index] = pBranch->GetLength();
    }
    // Otherwise, simply add current branch data to parent branch data
    else
    {
        mUpstreamPathBranchLengths[current_index] = pBranch->GetLength() + mUpstreamPathBranchLengths[pBranch->GetParent()->GetIndex()];
    }

    // Recursively calculate correct value for each child branch
    if (pBranch->GetChildOne() != nullptr)
    {
        RecursivelyCalculateUpstreamLengths(pBranch->GetChildOne());
    }
    if (pBranch->GetChildTwo() != nullptr)
    {
        RecursivelyCalculateUpstreamLengths(pBranch->GetChildTwo());
    }
}
void AirwayPropertiesCalculator::RecursivelyCalculateUpstreamVolumes(AirwayBranch* pBranch)
{
    unsigned current_index = pBranch->GetIndex();

    // If current branch is trachea, there is nothing upstream, so just need current branch data
    if (pBranch->GetParent() == nullptr)
    {
        mUpstreamPathBranchVolumes[current_index] = pBranch->GetBranchVolume();
    }
    // Otherwise, simply add current branch data to parent branch data
    else
    {
        mUpstreamPathBranchVolumes[current_index] = pBranch->GetBranchVolume() + mUpstreamPathBranchVolumes[pBranch->GetParent()->GetIndex()];
    }

    // Recursively calculate correct value for each child branch
    if (pBranch->GetChildOne() != nullptr)
    {
        RecursivelyCalculateUpstreamVolumes(pBranch->GetChildOne());
    }
    if (pBranch->GetChildTwo() != nullptr)
    {
        RecursivelyCalculateUpstreamVolumes(pBranch->GetChildTwo());
    }
}
void AirwayPropertiesCalculator::RecursivelyCalculateUpstreamLateralSurfaceAreas(AirwayBranch* pBranch)
{
    unsigned current_index = pBranch->GetIndex();

    // If current branch is trachea, there is nothing upstream, so just need current branch data
    if (pBranch->GetParent() == nullptr)
    {
        mUpstreamPathBranchLateralSurfaceAreas[current_index] = pBranch->GetBranchLateralSurfaceArea();
    }
    // Otherwise, simply add current branch data to parent branch data
    else
    {
        mUpstreamPathBranchLateralSurfaceAreas[current_index] = pBranch->GetBranchLateralSurfaceArea() + mUpstreamPathBranchLateralSurfaceAreas[pBranch->GetParent()->GetIndex()];
    }

    // Recursively calculate correct value for each child branch
    if (pBranch->GetChildOne() != nullptr)
    {
        RecursivelyCalculateUpstreamLateralSurfaceAreas(pBranch->GetChildOne());
    }
    if (pBranch->GetChildTwo() != nullptr)
    {
        RecursivelyCalculateUpstreamLateralSurfaceAreas(pBranch->GetChildTwo());
    }
}

void AirwayPropertiesCalculator::RecursivelyCalculateUpstreamPoiseuilleResistances(AirwayBranch* pBranch)
{
    unsigned current_index = pBranch->GetIndex();

    // If current branch is trachea, there is nothing upstream, so just need current branch data
    if (pBranch->GetParent() == nullptr)
    {
        mUpstreamPathPoiseuilleResistances[current_index] = pBranch->GetPoiseuilleResistance();
    }
    // Otherwise, simply add current branch data to parent branch data
    else
    {
        mUpstreamPathPoiseuilleResistances[current_index] = pBranch->GetPoiseuilleResistance() + mUpstreamPathPoiseuilleResistances[pBranch->GetParent()->GetIndex()];
    }

    // Recursively calculate correct value for each child branch
    if (pBranch->GetChildOne() != nullptr)
    {
        RecursivelyCalculateUpstreamPoiseuilleResistances(pBranch->GetChildOne());
    }
    if (pBranch->GetChildTwo() != nullptr)
    {
        RecursivelyCalculateUpstreamPoiseuilleResistances(pBranch->GetChildTwo());
    }
}
