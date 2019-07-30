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


#ifndef _PURKINJEVENTRICULARJUNCTIONSTIMULUS_HPP_
#define _PURKINJEVENTRICULARJUNCTIONSTIMULUS_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * Provides a stimulus between two cells dependent on the difference between their
 * transmembrane potential. This stimulus represents current flow across a Purkine-
 * Ventricular junction.
 *
 * I_pvj = (Vm - Vp)/R
 *
 * where Vm is the transmembrane potential in the ventricular myocyte, Vp is the transmembrane potential
 * in the Purkinje myocyte and R is the resistance across the junction.
 */
class PurkinjeVentricularJunctionStimulus : public AbstractStimulusFunction
{


private:
    /**
     * The resistance across the Purkinje-ventricular junction
     */
    double mJunctionResistance;

    /**
     * Flag to set if this Purkinje-ventricular junction is stimulating a ventricular cell model
     * or a Purkinje cell model (defaults to true)
     */
    bool mAppliedToVentricularCellModel;

    /**
     * Pointer to the cell model on the ventricular side of the junction
     */
    AbstractCardiacCellInterface* mpVentricularCellModel;

    /**
     * Pointer to the cell model on the Purkinje side of the junction
     */
    AbstractCardiacCellInterface* mpPurkinjeCellModel;

public:

    /**
     * Constructor.
     *
     * Note that Purkinje-ventricular junctions default to generating stimuli
     * for ventricular models, not Purkinje models.
     *
     * @param rJunctionResistance The resistance across the junction
     */
    PurkinjeVentricularJunctionStimulus(const double& rJunctionResistance);

    /**
     * Sets the pointer to the cell model on the ventricular side of the junction
     *
     * @param pVentricularModel Pointer to the ventricular cell model
     */
    void SetVentricularCellModel(AbstractCardiacCellInterface* pVentricularModel);

    /**
     * Sets the pointer to the cell model on the Purkinje side of the junction
     *
     * @param pPurkinjeModel Pointer to the Purkinje cell model
     */
    void SetPurkinjeCellModel(AbstractCardiacCellInterface* pPurkinjeModel);

    /**
     * Sets the Purkinje ventricular junction to generate stimuli for Purkinje cell models instead of
     * ventricular cell models.
     */
    void SetAppliedToPurkinjeCellModel();

    /**
     * @return the stimulus at a given time.
     *
     * @param time  time at which to return the stimulus (note that this is ignored in PurkinjeVentricularStimuli)
     */
    double GetStimulus(double time);
};

#endif //_PURKINJEVENTRICULARJUNCTIONSTIMULUS_HPP_

