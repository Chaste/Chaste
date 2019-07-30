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


#ifndef PLANESTIMULUSCELLFACTORY_HPP_
#define PLANESTIMULUSCELLFACTORY_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellFactory.hpp"
#include "LogFile.hpp"
#include "SimpleStimulus.hpp"

/**
 * PlaneStimulusCellFactory provides cells within 1e-5 of x=0 with a SimpleStimulus.
 */
template<class CELL, unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class PlaneStimulusCellFactory : public AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>
{
protected:
    /** The stimulus to apply at stimulated nodes */
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    /**
     * Constructor
     * @param stimulusMagnitude  The magnitude of the simple stimulus to be applied (defaults to -600).
     * @param stimulusDuration  The duration of the simple stimulus to be applied (defaults to 0.5ms).
     */
    PlaneStimulusCellFactory(double stimulusMagnitude=-600, double stimulusDuration=0.5)
        : AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>()
    {
        mpStimulus.reset(new SimpleStimulus(stimulusMagnitude, stimulusDuration));
        LOG(1, "Defined a PlaneStimulusCellFactory<"<<SPACE_DIM<<"> with SimpleStimulus("<<stimulusMagnitude<<","<< stimulusDuration<< ")\n");
    }

    /**
     * @param pNode  Pointer to the node.
     * @return  A cardiac cell which corresponds to this node.
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<SPACE_DIM>* pNode)
    {
        double x = pNode->GetPoint()[0];

        ///\todo remove magic number? (#1884)
        if (x*x<=1e-10)
        {
            return new CELL(this->mpSolver, mpStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


#endif /*PLANESTIMULUSCELLFACTORY_HPP_*/
