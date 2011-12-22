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

#ifndef ELECTRODESSTIMULSUFACTORY_HPP_
#define ELECTRODESSTIMULSUFACTORY_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedVector.hpp"
#include "AbstractChasteRegion.hpp"
#include "AbstractStimulusFactory.hpp"
#include "LinearBasisFunction.hpp"

/**
 * This class implements the specification of tow electrodes with a RegularStimulus applied to them.
 * It makes sure the compatibility conditions is verified by scaling the magnitude of the second electrode
 * according to the flux flowing thorugh each of them.
 *
 * User needs to pass in a vector of pairs of electrodes and, for each pair, the stimulation parameters.
 * This class will then implement the CreateStimulusForNode accordingly by creating a RegularStimulus object
 * for each node. Each node belonging to any "first electrode" will have an extracellular volume stimulus equals to the value that is passed in,
 * while each node belonging to any "second elecrode" will have an extracellular volume stimulus corrected by the flux (proportional to electrode volume)
 * in such a way that compatibility conditions of the extended bidomain equations are followed.
 *
 * It is possible to ground all second electrodes by calling GroundSecondElectrode. In this case, the member variable
 * mGroundSecondElectrode will have size > 0 and other classes can behave accordingly (e.g., problem class will set Dirichlet
 * boundary conditions accordingly).
 *
 */
template<unsigned DIM>
class ElectrodesStimulusFactory : public AbstractStimulusFactory<DIM>
{
protected:

    friend class TestStimulusFactory;//for testing

    /**
     * Vector of pairs, each pair representing a pair of electrodes.
     */
    std::vector<std::pair<AbstractChasteRegion<DIM>*, AbstractChasteRegion<DIM>*> >& mrElectrodePairs;

    /**
     * Vector of stimuli magnitudes (microA / cm^3). Must be the same size as mrElectrodePairs.
     */
    std::vector<double>& mrMagnitudes;
    /**
     * Vector of stimuli durations (ms). Must be the same size as mrElectrodePairs.
     */
    std::vector<double>& mrDurations;
    /**
     * Vector of stimuli periods (ms). Must be the same size as mrElectrodePairs.
     */
    std::vector<double>& mrPeriods;
    /**
     * Vector of stimuli start times (ms). Must be the same size as mrElectrodePairs.
     */
    std::vector<double>& mrStarts;
    /**
     * Vector of stimuli end times (ms). Must be the same size as mrElectrodePairs.
     */
    std::vector<double>& mrEnds;

    /**
     * Used to store temporarily magnitudes of electrode 1 for correction
     */
    std::vector<double> mMagnitudesElectrode1;

    /**
     * Used to store temporarily magnitudes of electrode 2 for correction
     */
    std::vector<double> mMagnitudesElectrode2;


    /** Whether the second electrode is grounded */
    bool mGroundSecondElectrode;

    /**
     * Helper method to check if any of the electrodes intersects with the other.
     * it needs to check whether any element that conatins any node of electrode 1
     * does not conatin also nodes of eledtrode 2
     */
    void CheckForElectrodesIntersection();

    /** typedef for basis function*/
    typedef LinearBasisFunction<DIM> BasisFunction;


    /**
     * Helper method to compute electrode flux. This method is the main functionality of this class.
     *
     * It essentially computes the contribution of each electrode to the RHS vector of the linear system.
     * The loops over elements and gauss points mimicks the ones in AssembleOnElement.
     * Interpolation is done for the value of the magnitude of the electrode.
     *
     * NOTE that this method assumes the use of a GaussianQuadratureRule with 2 Gauss points (the default value of Chaste FE assemblers).
     *
     * @param pRegion  the electrode
     * @param stimulusMagnitude the stimulus magnitude
     */
    double ComputeElectrodeTotalFlux(AbstractChasteRegion<DIM>* pRegion, double stimulusMagnitude);


public:


    /**
     * Constructor. Electrodes and stimulation parameters need to be passed in.
     *
     * @param rElectrodePairs the pairs of electrodes
     * @param rStimulusMagnitudes the magnitudes of the stimuli (microA / cm^3). First electrode will have magnitude value and second electrode will have -magnitude (before being corrected to ensure equal flux).
     * @param rDurations the duration of each stimulus (ms)
     * @param rPeriods the period of each stimulus (ms)
     * @param rStarts the start time of each stimulus (ms).
     * @param rEnds the end of each stimulation (ms)
     */
    ElectrodesStimulusFactory(std::vector<std::pair<AbstractChasteRegion<DIM>*, AbstractChasteRegion<DIM>*> >& rElectrodePairs,
                              std::vector<double>& rStimulusMagnitudes,
                              std::vector<double>& rDurations,
                              std::vector<double>& rPeriods,
                              std::vector<double>& rStarts,
                              std::vector<double>& rEnds);

    /**
     * Creates an appropriate stimuus for each node as to abide compatibility conditions.
     *
     * @param nodeIndex the node index for which to create the stimulus
     */
    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(unsigned nodeIndex);

    /**
     * Destructor
     */
    ~ElectrodesStimulusFactory();

    /**
     * This method checks whether the users wanted a grounded electrode.
     * If so, it fills in mGroundedRegions accordingly.
     *
     * If not, it calls ComputeElectrodeTotalFlux and scales the magnitudes of the electrodes stimuli accordingly
     */
    void SetCompatibleExtracellularStimulus();

    /**
     * Allow the user to ground the second electrode
     *
     * @param grounded : true if second electrode is grounded
     */
    void GroundSecondElectrode(bool grounded);

};

#endif /*ELECTRODESSTIMULSUFACTORY_HPP_*/
