/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_
#define ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

#include "Debug.hpp"

/**
 * An abstract immersed boundary force class, for use in 
 * immersed boundary cell-based simulations.
 * 
 * NOTE: while this class shares similarities with AbstractForce, 
 * it is not passed to a simulation object, but instead to an 
 * ImmersedBoundarySimulationModifier.
 */
template<unsigned DIM>
class AbstractImmersedBoundaryForce : public Identifiable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive& mMultiplicativeNormalNoise;
        archive& mNormalNoiseMean;
        archive& mNormalNoiseStdDev;
    }

protected:

    /**
     * Add the random noise to nodes, if mMultiplicativeNormalNoise is set.  This method handles randomization of the
     * forces in such a way as to ensure no net change to forces across the whole domain.
     *
     * @param rCellPopulation an immersed boundary cell population
     */
    void AddMultiplicativeNormalNoiseToNodes(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned elem_idx = 0; elem_idx < rCellPopulation.GetNumElements(); ++elem_idx)
        {
            unsigned num_nodes_this_elem = rCellPopulation.GetElement(elem_idx)->GetNumNodes();

            // Generate a vector with 0, 1, 2, ..., num_nodes-1
            std::vector<unsigned> random_node_order;
            random_node_order.reserve(num_nodes_this_elem);
            for (unsigned node_idx = 0; node_idx < num_nodes_this_elem; ++node_idx)
            {
                random_node_order.push_back(node_idx);
            }

            // Shuffle the vector
            std::random_shuffle(random_node_order.begin(), random_node_order.end());

            // Only go to half way (rounded down - forget about the very last node for now)
            unsigned half_way = num_nodes_this_elem / 2;
            for (unsigned i = 0; i < half_way; ++i)
            {
                unsigned rand_node_a = random_node_order[2 * i];
                unsigned rand_node_b = random_node_order[2 * i + 1];

                c_vector<double, DIM> random_force;
                for (unsigned dim = 0; dim < DIM; ++dim)
                {
                    random_force[dim] = p_gen->NormalRandomDeviate(mNormalNoiseMean, mNormalNoiseStdDev);
                }

                double avg_magnitude = 0.5 * (norm_2(rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_a)->rGetAppliedForce()) +
                                              norm_2(rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_b)->rGetAppliedForce()));

                random_force *= avg_magnitude;
                rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_a)->AddAppliedForceContribution(random_force);

                random_force *= -1.0;
                rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_b)->AddAppliedForceContribution(random_force);
            }
        }
    }

    /** Whether to apply multiplicative normal noise to the calculated force */
    bool mMultiplicativeNormalNoise;

    /** The mean of the Normal distribution from which random noise variations are drawn */
    double mNormalNoiseMean;

    /** The standard deviation of the Normal distribution from which random noise variations are drawn */
    double mNormalNoiseStdDev;

public:

    /**
     * Default constructor.
     */
    AbstractImmersedBoundaryForce();

    /**
     * Destructor.
     */
    virtual ~AbstractImmersedBoundaryForce();

    /**
     * Calculates the force on each immersed boundary node.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to contribute the force
     * @param rCellPopulation an immersed boundary cell population
     */
    virtual void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
            ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Outputs the name of the immersed boundary force used in 
     * the simulation to file and then calls OutputImmersedBoundaryForceParameters()
     * to output all relevant parameters.
     *
     * \todo At present, this method is not used in simulations. In contrast, the force method
     * OutputForceInfo is called in OffLatticeSimulation::OutputAdditionalSimulationSetup().
     * Call this method in simulations, or else remove it.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceInfo(out_stream& rParamsFile);

    /**
     * Outputs force parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)=0;

    /** @return mMultiplicativeNormalNoise */
    bool GetMultiplicativeNormalNoise() const;

    /** @param multiplicativeNormalNoise whether to include multiplicative normal noise */
    void SetMultiplicativeNormalNoise(bool multiplicativeNormalNoise);

    /** @return mNormalNoiseMean */
    double GetNormalNoiseMean() const;

    /** @param normalNoiseMean the new value of mNormalNoiseMean */
    void SetNormalNoiseMean(double normalNoiseMean);

    /** @return mNormalNoiseStdDev */
    double GetNormalNoiseStdDev() const;

    /** @param normalNoiseStdDev the new value of mNormalNoiseStdDev */
    void SetNormalNoiseStdDev(double normalNoiseStdDev);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractImmersedBoundaryForce)

#endif /*ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_*/
