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

#ifndef CRYPTPROJECTIONFORCE_HPP_
#define CRYPTPROJECTIONFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law for use in crypt projection simulations.
 */
class CryptProjectionForce : public GeneralisedLinearSpringForce<2>
{
    friend class TestCryptProjectionForce;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<2> >(*this);
        archive & mA;
        archive & mB;
        archive & mIncludeWntChemotaxis;
        archive & mWntChemotaxisStrength;
    }

    /**
     * The value of the constant a in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mA;

    /**
     * The value of the constant b in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mB;

    /**
     * Whether to include Wnt-dependent chemotaxis for stem cells.
     */
    bool mIncludeWntChemotaxis;

    /**
     * Strength of Wnt-based chemotactic force.
     */
    double mWntChemotaxisStrength;

    /**
     * Map node indices to 3D locations on the crypt surface.
     */
    std::map<unsigned, c_vector<double, 3> > mNode3dLocationMap;

    /**
     * Fix up the mappings between node indices and 3D locations.
     *
     * @param rCellPopulation the cell population
     */
    void UpdateNode3dLocationMap(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex the index of the first node
     * @param nodeBGlobalIndex the index of the second node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double,2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<2>& rCellPopulation);

public:

    /**
     * Constructor.
     */
    CryptProjectionForce();

    /**
     * Destructor.
     */
    ~CryptProjectionForce();

    /**
     * @return mA.
     */
    double GetA() const;

    /**
     * @return mB.
     */
    double GetB() const;

    /**
     * @return mWntChemotaxisStrength
     */
    double GetWntChemotaxisStrength();

    /**
     * Set mWntChemotaxisStrength.
     *
     * @param wntChemotaxisStrength the new value of mWntChemotaxisStrength
     */
    void SetWntChemotaxisStrength(double wntChemotaxisStrength);

    /**
     * Set mIncludeWntChemotaxis.
     *
     * @param includeWntChemotaxis whether to include Wnt-dependent chemotaxis
     */
    void SetWntChemotaxis(bool includeWntChemotaxis);

    /**
     * Calculates the height of the crypt surface given by
     *     z = f(r) = a*r^b
     * at a point whose 2D position is a distance r from the centre of the cell population.
     * This assumes that the cell population is centred at the origin.
     *
     * @param rNodeLocation
     * @return the z component corresponding to rNodeLocation
     */
    double CalculateCryptSurfaceHeightAtPoint(const c_vector<double,2>& rNodeLocation);

    /**
     * Calculates the derivative df/dr of the crypt surface function z=f(r) at a point
     * whose 2D position is a distance r from the centre of the cell_population, which we assume
     * to be at (0,0).
     *
     * @param rNodeLocation the 2D location of a node
     * @return the gradient
     */
    double CalculateCryptSurfaceDerivativeAtPoint(const c_vector<double,2>& rNodeLocation);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptProjectionForce)

#endif /*CRYPTPROJECTIONFORCE_HPP_*/
