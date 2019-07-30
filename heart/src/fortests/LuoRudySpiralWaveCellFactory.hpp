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

#ifndef LUORUDYSPIRALWAVECELLFACTORY_HPP_
#define LUORUDYSPIRALWAVECELLFACTORY_HPP_

#include "LuoRudy1991BackwardEuler.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "MultiStimulus.hpp"
#include "SimpleStimulus.hpp"

/**
 * This class sets up Luo Rudy 1991 backward Euler cells with
 * a change to parameters (G_si = 0) and S1-S2 stimulus pattern
 * that together provoke spiral wave initiation.
 */
class LuoRudySpiralWaveCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    /** The stimulus to be applied at the S1 site (simple one-off stimulus). */
    boost::shared_ptr<SimpleStimulus> mpS1Stimulus;
    /** The stimulus to be applied at the S2 site (simple one-off stimulus). */
    boost::shared_ptr<SimpleStimulus> mpS2Stimulus;
    /** Both of the S1 and S2 stimuli. */
    boost::shared_ptr<MultiStimulus> mpBothStimulus;
    /** Width of the domain. */
    double mXExtent;
    /** Height of the domain. */
    double mYExtent;
public:
    /**
     * Constructor for spiral wave cell factory.
     *
     * @param xExtent  The size of the domain (width)
     * @param yExtent  The size of the domain (height)
     */
    LuoRudySpiralWaveCellFactory(double xExtent, double yExtent)
        : AbstractCardiacCellFactory<2>(),
          mpS1Stimulus(new SimpleStimulus(-50000.0, 2, 0)),
          mpS2Stimulus(new SimpleStimulus(-70000.0, 2, 45)),
          mpBothStimulus(new MultiStimulus()),
          mXExtent(xExtent),
          mYExtent(yExtent)
    {
        mpBothStimulus->AddStimulus(mpS1Stimulus);
        mpBothStimulus->AddStimulus(mpS2Stimulus);
    }

    /**
     * This method is provided by all factories, and returns an action potential model
     * with the correct stimulus for that location. Here we set up some cells with S1 and S2,
     * some with S1 and some with S2, to initiate a spiral wave. We also alter one of the
     * conductances in the model to make the wave easier to initiate on this small domain.
     *
     * @param pNode  the pointer to Node object we would like an action potential model for
     * @return a cardiac cell for this node
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];

        double x_threshold_for_S1 = 0.1 + 1e-6;
        double x_threshold_for_S2 = mXExtent*0.6;
        double y_threshold_for_S2 = mYExtent*0.5;

        AbstractCardiacCell* p_cell;
        if (x < x_threshold_for_S1)
        {
            if (y < y_threshold_for_S2)
            {
                p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpBothStimulus);
            }
            else
            {
                p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpS1Stimulus);
            }
        }
        else if ((x < x_threshold_for_S2) && (y < y_threshold_for_S2))
        {
            p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpS2Stimulus);
        }
        else
        {
            p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, this->mpZeroStimulus);
        }

        // Alter parameters to match those used in
        // Qu et al. Origins of Spiral Wave Meander and Breakup... Annals of Biomedical Eng. 28:755-771 (2000).
        p_cell->SetParameter("membrane_L_type_calcium_current_conductance",0); // slow inward current in original model.

        unsigned node_index = pNode->GetIndex();
        if (node_index==0)
        {
//            std::cout << "[K_o] = " << p_cell->GetAnyVariable("extracellular_potassium_concentration") << "mM \n";
//            std::cout << "[K_i] = " << p_cell->GetAnyVariable("cytosolic_potassium_concentration") << "mM \n";
//            std::cout << "[Na_o] = " << p_cell->GetAnyVariable("extracellular_sodium_concentration") << "mM \n";
//            std::cout << "[Na_i] = " << p_cell->GetAnyVariable("cytosolic_sodium_concentration") << "mM \n";
            std::cout << "G_si = " << p_cell->GetAnyVariable("membrane_L_type_calcium_current_conductance") << "\n";
//            std::cout << "G_Na = " << p_cell->GetAnyVariable("membrane_fast_sodium_current_conductance") << "mM \n";
//            std::cout << "G_K = " << p_cell->GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current_conductance") << "mM \n";
//            std::cout << "G_K1 = " << p_cell->GetAnyVariable("membrane_inward_rectifier_potassium_current_conductance") << "mM \n";
        }
        return p_cell;
    }
};

#endif // LUORUDYSPIRALWAVECELLFACTORY_HPP_
