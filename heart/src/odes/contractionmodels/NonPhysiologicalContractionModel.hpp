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

#ifndef NONPHYSIOLOGICALCONTRACTIONMODEL_HPP_
#define NONPHYSIOLOGICALCONTRACTIONMODEL_HPP_

#include "AbstractAlgebraicContractionModel.hpp"


/**
 *  3 algebraic (non-ODE based), non-physiological contraction models. The returned active
 *  tension, sigma_a, is:
 *  (1) sigma_a = |5 sin(t/4)|
 *  (2) sigma_a = |5 lam sin(t/4)|
 *  (3) sigma_a = |5 exp(1-lam) sin(t/4)|
 *  where lam is the fibre stretch.
 */
class NonPhysiologicalContractionModel : public AbstractAlgebraicContractionModel
{
private:
    /** Number between 1 and 3 indicating which model is to be used */
    unsigned mOption;
    /** Fibre stretch */
    double mStretch;

public:
    /** Constructor.
     *  @param option Number between 1 and 3 which determines which model is used. See class documentations
     */
    NonPhysiologicalContractionModel(unsigned option);

    /**
     *  Set the input parameters (none of which are used in this model
     *  @param rInputParameters Structure containing voltage, [Ca]_i.
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters);

    /** Set the fibre stretch and stretch-rate (only stretch is used)
     *  @param stretch stretch in the fibre direction
     *  @param stretchRate rate of change of stretch in the fibre direction (not used in this model).
     */
    void SetStretchAndStretchRate(double stretch, double stretchRate);

    /**
     *  Get the active tension given the current stretch and time
     */
    double GetActiveTension();

    /**
     *  This model is stretch-rate-independent
     */
    bool IsStretchDependent()
    {
        return true;
    }

    /**
     *  This model is stretch-rate-independent
     */
    bool IsStretchRateDependent()
    {
        return false;
    }
};


#endif /*NONPHYSIOLOGICALCONTRACTIONMODEL_HPP_*/

