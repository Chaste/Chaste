/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef PATHMANATHANINTERACTIONFORCE_HPP_
#define PATHMANATHANINTERACTIONFORCE_HPP_

#include "AbstractVariableSizeTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A two-body force law implementing the model from Pathmanathan et al (2009)
 * (doi:10.1088/1478-3975/6/3/036001), also described in Osborne et al (2017)
 * (doi:10.1242/dev.126359).
 *
 * For two cells whose separation is less than the sum of their radii
 * (i.e. overlap < 0, cells are compressed), a logarithmic repulsion force is used:
 * \f[
 * \mathbf{F} = \mu \hat{\mathbf{r}} s \ln\!\left(1 + \frac{d - s}{s}\right)
 * \f]
 *
 * For two cells whose separation is greater than the sum of their radii
 * (i.e. overlap < 0, cells are stretched), an exponential attraction force is used:
 * \f[
 * \mathbf{F} = \mu \hat{\mathbf{r}} (d - s) e^{-\alpha(d-s)/s}
 * \f]
 *
 * Here \f$\mu\f$ is the spring stiffness, \f$s\f$ is the rest length (sum of cell radii),
 * \f$d\f$ is the distance between cell centres, and \f$\alpha\f$ controls the
 * range of attraction.
 *
 * This force is designed for use with NodeBasedCellPopulation.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class PathmanathanInteractionForce : public AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mAlpha;
    }

protected:

    /**
     * Parameter controlling the range of attraction between cells.
     * A larger value makes the attractive force decay more rapidly with distance.
     * Defaults to 5.0.
     */
    double mAlpha;

    /**
     * Overridden CalculateLinkInteraction() method.
     *
     * Calculates the Pathmanathan force law expression after shared
     * rest-length mechanics have been computed by AbstractVariableSizeTwoBodyInteractionForce.
     *
     *  @param overlap the amount by which the distance between nodes is less than the rest length
     *  @param restLength the rest length of the spring between the nodes
     *  @param rUnitDifference the unit vector pointing from one node to the other
     *  @param multiplicationFactor a multiplication factor for the spring constant
     *  @return the force vector between the two nodes
     */
    c_vector<double, SPACE_DIM> CalculateLinkInteraction(double overlap,
                             double restLength,
                             const c_vector<double, SPACE_DIM>& rUnitDifference,
                             double multiplicationFactor);

public:

    /**
     * Constructor.
     */
    PathmanathanInteractionForce();

    /**
     * Destructor.
     */
    virtual ~PathmanathanInteractionForce();

    /**
     * @return mAlpha
     */
    double GetAlpha();

    /**
     * Set mAlpha.
     *
     * @param alpha the new value of mAlpha
     */
    void SetAlpha(double alpha);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PathmanathanInteractionForce)

#endif /*PATHMANATHANINTERACTIONFORCE_HPP_*/
