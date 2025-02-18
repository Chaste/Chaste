/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYSVGWRITER_HPP_
#define IMMERSEDBOUNDARYSVGWRITER_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellPopulation.hpp"

/**
 * A modifier class which after each simulation time step exports the simulation geometry to a scalable vector graphic
 * (SVG) file.
 */
template <unsigned DIM>
class ImmersedBoundarySvgWriter : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    /** The sampling frequency for exporting svg frames */
    unsigned mSamplingMultiple;

    /** The width and height in pixels of the svg file */
    double mSvgSize;

    /** Set to the output directory during SetupSolve() */
    std::string mOutputDirectory;

    /** The svg file header, which will be constant for a given simulation */
    std::string mSvgHeader;

    /** The svg file footer, which will be constant for a given simulation */
    std::string mSvgFooter;

    friend class TestImmersedBoundarySvgWriter;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& mSamplingMultiple;
        archive& mSvgSize;
        archive& mOutputDirectory;
        archive& mSvgHeader;
        archive& mSvgFooter;
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
    }

public:
    /**
     * Default constructor.
     */
    ImmersedBoundarySvgWriter();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundarySvgWriter();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method for UpdateAtEndOfTimeStep().
     *
     * Add a circle representing a point location to the svg file.
     *
     * @param rSvgFile the svg file stream to add to
     * @param location the location of the point
     * @param region metadata to allow different colouring for points
     * @param rad the radius of circle representing the point
     */
    void AddPointToSvgFile(out_stream& rSvgFile, c_vector<double, DIM> location, unsigned region, double rad);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /** @return mSamplingMultiple **/
    unsigned GetSamplingMultiple() const;

    /** @param samplingMultiple the new value of mSamplingMultiple */
    void SetSamplingMultiple(unsigned samplingMultiple);

    /** @return mSvgSize **/
    double GetSvgSize() const;

    /** @param svgSize the new value of mSvgSize */
    void SetSvgSize(double svgSize);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySvgWriter)

#endif /*IMMERSEDBOUNDARYSVGWRITER_HPP_*/
