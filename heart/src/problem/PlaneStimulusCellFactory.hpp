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
     * @param node  The global index of a node
     * @return  A cardiac cell which corresponds to this node.
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        double x = this->GetMesh()->GetNode(node)->GetPoint()[0];

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
