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

#ifndef VTKSCENEMODIFIER_HPP_
#define VTKSCENEMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "VtkScene.hpp"

/**
 * A cell-based simulation modifier that renders the simulation with a VtkScene.
 *
 * Holds a VtkScene and, at each time step that is a multiple of the update
 * frequency, has the scene render the current cell population, writing a
 * per-step image/animation frame and/or updating an interactive window.
 * 1D simulations are not rendered.
 */
template <unsigned DIM>
class VtkSceneModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
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
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
    }

    /**
     * The scene that renders the cell population
     */
    std::shared_ptr<VtkScene<DIM> > mpScene;

    /**
     * The scene is rendered every mUpdateFrequency time steps
     */
    unsigned mUpdateFrequency;

public:
    /**
     * Default constructor.
     */
    VtkSceneModifier();

    /**
     * Destructor.
     */
    virtual ~VtkSceneModifier() = default;

    /**
     * @return the number of time steps between scene renders
     */
    unsigned GetUpdateFrequency() const;

    /**
     * @return the scene used to render the cell population
     */
    std::shared_ptr<VtkScene<DIM> > GetVtkScene();

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /**
     * Set how often the scene is rendered.
     * @param frequency the number of time steps between renders
     */
    void SetUpdateFrequency(unsigned frequency);

    /**
     * Set the scene used to render the cell population.
     * @param pScene the scene
     */
    void SetVtkScene(std::shared_ptr<VtkScene<DIM> > pScene);

    /**
     * Overridden SetupSolve() method.
     *
     * Called before the time loop: updates the cell population and renders the
     * initial scene.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Updates the cell population and, if the step is due, renders the scene.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Ensure the cell population is up to date before the scene is rendered.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

private:
    /**
     * Render the scene for the current time step, if a scene has been set and
     * the step is a multiple of the update frequency.
     */
    void RenderSceneIfDue();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VtkSceneModifier)

#endif // VTKSCENEMODIFIER_HPP_
