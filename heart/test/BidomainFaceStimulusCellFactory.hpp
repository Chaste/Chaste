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


#ifndef BIDOMAINFACESTIMULUSCELLFACTORY_HPP_
#define BIDOMAINFACESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"

class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    boost::shared_ptr<RegularStimulus> mpRegStimulus;

public:
    //Pdetime step is (by default) 0.01
    //Odetime step set below to 0.001 (10:1)
    BidomainFaceStimulusCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpRegStimulus(new RegularStimulus(-900.0*1000, 0.5, 100.0, 0.0))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpRegStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

#endif /*BIDOMAINFACESTIMULUSCELLFACTORY_HPP_*/
