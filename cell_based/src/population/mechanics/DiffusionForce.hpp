/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef DIFFUSIONFORCE_HPP_
#define DIFFUSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A force class to model diffusion.
 */
template<unsigned DIM>
class DiffusionForce : public AbstractForce<DIM>
{
private :

    /**
     * Diffusion constant.
     */
    double mDiffusionConstant;

    /**
     * Absolute temperature (in Kelvin).
     */
    double mAbsoluteTemperature;

    /**
     * Viscosity of media. We assume that this is measured in units of
     * kg microns^(-1) h^(-1), and that cell diameters are scaled with
     * a characteristic length of 1 micron.
     */
    double mViscosity;

    /**
     * Mechanics cut off length.
     */
    double mMechanicsCutOffLength;

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDiffusionConstant;
        archive & mAbsoluteTemperature;
        archive & mViscosity;
        archive & mMechanicsCutOffLength;
    }

public :

    /**
     * Constructor.
     */
    DiffusionForce();

    /**
     * Destructor.
     */
    ~DiffusionForce();

    /**
     * Set the diffusion constant for the cells.
     *
     * @param diffusionConstant the diffusion constant to use
     */
    void SetDiffusionConstant(double diffusionConstant);

    /**
     * Set the absolute temperature, which affects the
     * diffusion constant.
     *
     * @param absoluteTemperature the temperature in Kelvin
     */
    void SetAbsoluteTemperature(double absoluteTemperature);

    /**
      * Set the media viscosity (dynamic), which affects the
      * diffusion constant.
      *
      * @param viscosity the viscosity
      */
    void SetViscosity(double viscosity);

    /**
     * Use a cutoff point, i.e. specify zero force if two cells are greater
     * than the cutoff distance apart.
     *
     * @param cutOffLength the cutoff to use
     */
    void SetCutOffLength(double cutOffLength);

    /**
     * Get the diffusion coefficient.
     *
     * @return mDiffusionConstant
     */
    double GetDiffusionConstant();

    /**
     * Get the absolute temperature.
     *
     * @return mAbsoluteTemperature
     */
    double GetAbsoluteTemperature();

    /**
      * Get the viscosity.
      *
      * @return mViscosity
      */
    double GetViscosity();

    /**
     * Get the cutoff length.
     *
     * @return mCutOffLenght
     */
    double GetCutOffLength();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
            AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionForce)

#endif /*DIFFUSIONFORCE_HPP_*/
