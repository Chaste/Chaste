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
     * NOTE that this method assumes the use of a 2nd order Gauss quadrature rule (the default value of Chaste FE assemblers).
     *
     * @param pRegion  the electrode
     * @param stimulusMagnitude the stimulus magnitude
     * @return electrode flux
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
     * @return Creates an appropriate stimulus for each node as to abide compatibility conditions.
     *
     * @param pNode the node for which to create the stimulus
     */
    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<DIM>* pNode);

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
