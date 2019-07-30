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

#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"

template<unsigned DIM>
LinearSpringWithVariableSpringConstantsForce<DIM>::LinearSpringWithVariableSpringConstantsForce()
    : GeneralisedLinearSpringForce<DIM>(),
      mUseEdgeBasedSpringConstant(false),
      mUseMutantSprings(false),
      mMutantMutantMultiplier(DOUBLE_UNSET),
      mNormalMutantMultiplier(DOUBLE_UNSET),
      mUseBCatSprings(false),
      mUseApoptoticSprings(false),
      mBetaCatSpringScaler(18.14/6.0), // scale spring constant with beta-catenin level (divided by 6 for heaxagonal cells)
      mApoptoticSpringTensionStiffness(15.0*0.25),
      mApoptoticSpringCompressionStiffness(15.0*0.75)
{
}

template<unsigned DIM>
LinearSpringWithVariableSpringConstantsForce<DIM>::~LinearSpringWithVariableSpringConstantsForce()
{
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2); // LCOV_EXCL_LINE
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier, double normalMutantMultiplier)
{
    mUseMutantSprings = useMutantSprings;
    mMutantMutantMultiplier = mutantMutantMultiplier;
    mNormalMutantMultiplier = normalMutantMultiplier;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetBetaCateninSprings(bool useBCatSprings)
{
    mUseBCatSprings = useBCatSprings;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetApoptoticSprings(bool useApoptoticSprings)
{
    mUseApoptoticSprings = useApoptoticSprings;
}

template<unsigned DIM>
double LinearSpringWithVariableSpringConstantsForce<DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{

    double multiplication_factor = GeneralisedLinearSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex,
                                                                                                            nodeBGlobalIndex,
                                                                                                            rCellPopulation,
                                                                                                            isCloserThanRestLength);

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    /*
     * The next code block computes the edge-dependent spring constant as given by equation
     * (3) in the following reference: van Leeuwen et al. 2009. An integrative computational model
     * for intestinal tissue renewal. Cell Prolif. 42(5):617-636. doi:10.1111/j.1365-2184.2009.00627.x
     */
    if (mUseEdgeBasedSpringConstant)
    {
        assert(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)));
        assert(!mUseBCatSprings);   // don't want to do both (both account for edge length)

        multiplication_factor = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->GetVoronoiEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex)*sqrt(3.0);
    }

    if (mUseMutantSprings)
    {
        unsigned number_of_mutants = 0;

        if (p_cell_A->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || p_cell_A->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
        {
            // If cell A is mutant
            number_of_mutants++;
        }

        if (p_cell_B->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || p_cell_B->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
        {
            // If cell B is mutant
            number_of_mutants++;
        }

        switch (number_of_mutants)
        {
            case 1u:
            {
                multiplication_factor *= mNormalMutantMultiplier;
                break;
            }
            case 2u:
            {
                multiplication_factor *= mMutantMutantMultiplier;
                break;
            }
        }
    }

    /*
     * The next code block computes the beta-catenin dependent spring constant as given by equation
     * (4) in the following reference: van Leeuwen et al. 2009. An integrative computational model
     * for intestinal tissue renewal. Cell Prolif. 42(5):617-636. doi:10.1111/j.1365-2184.2009.00627.x
     */
    if (mUseBCatSprings)
    {
        assert(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)));

        // If using beta-cat dependent springs, both cell-cycle models had better be VanLeeuwen2009WntSwatCellCycleModel
        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model_A = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(p_cell_A->GetCellCycleModel());
        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model_B = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(p_cell_B->GetCellCycleModel());

        assert(!mUseEdgeBasedSpringConstant);   // This already adapts for edge lengths - don't want to do it twice.
        double beta_cat_cell_1 = p_model_A->GetMembraneBoundBetaCateninLevel();
        double beta_cat_cell_2 = p_model_B->GetMembraneBoundBetaCateninLevel();

        MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation));

        double perim_cell_1 = p_static_cast_cell_population->GetSurfaceAreaOfVoronoiElement(nodeAGlobalIndex);
        double perim_cell_2 = p_static_cast_cell_population->GetSurfaceAreaOfVoronoiElement(nodeBGlobalIndex);
        double edge_length_between_1_and_2 = p_static_cast_cell_population->GetVoronoiEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex);

        double beta_cat_on_cell_1_edge = beta_cat_cell_1 *  edge_length_between_1_and_2 / perim_cell_1;
        double beta_cat_on_cell_2_edge = beta_cat_cell_2 *  edge_length_between_1_and_2 / perim_cell_2;

        double min_beta_Cat_of_two_cells = std::min(beta_cat_on_cell_1_edge, beta_cat_on_cell_2_edge);

        multiplication_factor *= min_beta_Cat_of_two_cells / mBetaCatSpringScaler;
    }

    if (mUseApoptoticSprings)
    {
        bool cell_A_is_apoptotic = p_cell_A->HasCellProperty<ApoptoticCellProperty>();
        bool cell_B_is_apoptotic = p_cell_B->HasCellProperty<ApoptoticCellProperty>();

        if (cell_A_is_apoptotic || cell_B_is_apoptotic)
        {
            double spring_a_stiffness = 2.0 * this->GetMeinekeSpringStiffness();
            double spring_b_stiffness = 2.0 * this->GetMeinekeSpringStiffness();

            if (cell_A_is_apoptotic)
            {
                if (!isCloserThanRestLength) // if under tension
                {
                    spring_a_stiffness = mApoptoticSpringTensionStiffness;
                }
                else // if under compression
                {
                    spring_a_stiffness = mApoptoticSpringCompressionStiffness;
                }
            }
            if (cell_B_is_apoptotic)
            {
                if (!isCloserThanRestLength) // if under tension
                {
                    spring_b_stiffness = mApoptoticSpringTensionStiffness;
                }
                else // if under compression
                {
                    spring_b_stiffness = mApoptoticSpringCompressionStiffness;
                }
            }

            multiplication_factor /= (1.0/spring_a_stiffness + 1.0/spring_b_stiffness)*this->GetMeinekeSpringStiffness();
        }
    }

    return multiplication_factor;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a MeshBasedCellPopulation
    if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("LinearSpringWithVariableSpringConstantsForce is to be used with a subclass of MeshBasedCellPopulation only");
    }

    MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

    for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
        spring_iterator != p_static_cast_cell_population->SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = this->CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);
        c_vector<double, DIM> negative_force = -1.0*force;

        spring_iterator.GetNodeB()->AddAppliedForceContribution(negative_force);
        spring_iterator.GetNodeA()->AddAppliedForceContribution(force);
    }
}

template<unsigned DIM>
double LinearSpringWithVariableSpringConstantsForce<DIM>::GetBetaCatSpringScaler()
{
    return mBetaCatSpringScaler;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetBetaCatSpringScaler(double betaCatSpringScaler)
{
    assert(betaCatSpringScaler > 0.0);
    mBetaCatSpringScaler = betaCatSpringScaler;
}

template<unsigned DIM>
double LinearSpringWithVariableSpringConstantsForce<DIM>::GetApoptoticSpringTensionStiffness()
{
    return mApoptoticSpringTensionStiffness;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness)
{
    assert(apoptoticSpringTensionStiffness >= 0.0);
    mApoptoticSpringTensionStiffness = apoptoticSpringTensionStiffness;
}

template<unsigned DIM>
double LinearSpringWithVariableSpringConstantsForce<DIM>::GetApoptoticSpringCompressionStiffness()
{
    return mApoptoticSpringCompressionStiffness;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness)
{
    assert(apoptoticSpringCompressionStiffness >= 0.0);
    mApoptoticSpringCompressionStiffness = apoptoticSpringCompressionStiffness;
}

template<unsigned DIM>
void LinearSpringWithVariableSpringConstantsForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseEdgeBasedSpringConstant>" << mUseEdgeBasedSpringConstant << "</UseEdgeBasedSpringConstant>\n";
    *rParamsFile << "\t\t\t<UseMutantSprings>" << mUseMutantSprings << "</UseMutantSprings>\n";
    *rParamsFile << "\t\t\t<MutantMutantMultiplier>" << mMutantMutantMultiplier << "</MutantMutantMultiplier>\n";
    *rParamsFile << "\t\t\t<NormalMutantMultiplier>" << mNormalMutantMultiplier << "</NormalMutantMultiplier>\n";
    *rParamsFile << "\t\t\t<UseBCatSprings>" << mUseBCatSprings << "</UseBCatSprings>\n";
    *rParamsFile << "\t\t\t<UseApoptoticSprings>" << mUseApoptoticSprings << "</UseApoptoticSprings>\n";
    *rParamsFile << "\t\t\t<BetaCatSpringScaler>" << mBetaCatSpringScaler << "</BetaCatSpringScaler>\n";
    *rParamsFile << "\t\t\t<ApoptoticSpringTensionStiffness>" << mApoptoticSpringTensionStiffness << "</ApoptoticSpringTensionStiffness>\n";
    *rParamsFile << "\t\t\t<ApoptoticSpringCompressionStiffness>" << mApoptoticSpringCompressionStiffness << "</ApoptoticSpringCompressionStiffness>\n";

    // Call method on direct parent class
    GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class LinearSpringWithVariableSpringConstantsForce<1>;
template class LinearSpringWithVariableSpringConstantsForce<2>;
template class LinearSpringWithVariableSpringConstantsForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LinearSpringWithVariableSpringConstantsForce)
