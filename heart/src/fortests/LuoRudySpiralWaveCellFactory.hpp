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

#ifndef LUORUDYSPIRALWAVECELLFACTORY_HPP_
#define LUORUDYSPIRALWAVECELLFACTORY_HPP_

#include "LuoRudy1991BackwardEuler.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "MultiStimulus.hpp"

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
     * @param nodeIndex  the global node index we would like an action potential model for
     * @return a cardiac cell for this node
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];

        double x_threshold_for_S1 = 0.1 + 1e-6;
        double x_threshold_for_S2 = mXExtent*0.6;
        double y_threshold_for_S2 = mYExtent*0.5;

        AbstractCardiacCell* p_cell;
        if ( x < x_threshold_for_S1 )
        {
            if (y<y_threshold_for_S2)
            {
                p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpBothStimulus);
            }
            else
            {
                p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpS1Stimulus);
            }
        }
        else if ( (x < x_threshold_for_S2) && (y < y_threshold_for_S2) )
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

        if (nodeIndex==0)
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
