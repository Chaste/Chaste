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

#include "VtkScene.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr/make_shared.hpp>

#include <vtkActor.h>
#include <vtkActorCollection.h>
#include <vtkAutoInit.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPNGWriter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVersion.h>
#include <vtkWindowToImageFilter.h>

// Initialise the VTK rendering backend's object factories. VTK_MODULE_INIT
// defines a static initialiser, so it must live in exactly one translation
// unit (here, not in a header) to register the modules once for the library.
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);

namespace
{
/** Default render window dimensions, in pixels. */
constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 600;

/** Camera azimuth applied to 3D scenes, in degrees. */
constexpr double DEFAULT_3D_AZIMUTH_DEGREES = 45.0;

/** Animation frame rate, in frames per second. */
constexpr double ANIMATION_FRAME_RATE = 1.0;
} // namespace

template <unsigned DIM>
VtkScene<DIM>::VtkScene()
        : mpRenderer(vtkSmartPointer<vtkRenderer>::New()),
          mpRenderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
          mpRenderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
          mAnimationWriter(vtkSmartPointer<vtkOggTheoraWriter>::New()),
          mWindowToImageFilter(vtkSmartPointer<vtkWindowToImageFilter>::New()),
          mpCellPopulationGenerator(boost::make_shared<CellPopulationPyChasteActorGenerator<DIM> >())
{
    mpRenderer->SetBackground(1.0, 1.0, 1.0); // white
    mpRenderWindow->AddRenderer(mpRenderer);
    mpRenderWindow->SetSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    mpRenderWindowInteractor->SetRenderWindow(mpRenderWindow);

    // An explicit interactor style is needed to capture mouse events
    auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    mpRenderWindowInteractor->SetInteractorStyle(style);
}

template <unsigned DIM>
VtkScene<DIM>::~VtkScene()
{
    End();
}

template <unsigned DIM>
void VtkScene<DIM>::End()
{
    if (mAnimationWriterStarted)
    {
        mAnimationWriter->End();
        mAnimationWriterStarted = false;
    }
}

template <unsigned DIM>
boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > VtkScene<DIM>::GetCellPopulationActorGenerator() const
{
    return mpCellPopulationGenerator;
}

template <unsigned DIM>
vtkSmartPointer<vtkRenderWindow> VtkScene<DIM>::GetRenderWindow() const
{
    return mpRenderWindow;
}

template <unsigned DIM>
vtkSmartPointer<vtkRenderer> VtkScene<DIM>::GetRenderer() const
{
    return mpRenderer;
}

template <unsigned DIM>
vtkSmartPointer<vtkUnsignedCharArray> VtkScene<DIM>::GetSceneAsCharBuffer()
{
    // Rebuild and render the scene to an in-memory PNG, without writing any
    // output files (unlike RenderFrame).
    if (!mHasStarted)
    {
        Start();
    }
    else
    {
        UpdateScene();
    }

    RenderForCapture();

    auto p_writer = MakePngWriter();
    p_writer->Write();

    return p_writer->GetResult();
}

template <unsigned DIM>
void VtkScene<DIM>::RenderFrame(unsigned timeStep)
{
    // The first call lazily runs Start() (which builds the scene); later calls
    // just rebuild the actors to reflect the current cell population.
    if (!mHasStarted)
    {
        Start();
    }
    else
    {
        UpdateScene();
    }

    // Write this step's output: a numbered PNG (images) and/or an animation frame.
    if (mSaveAsImages || mSaveAsAnimation)
    {
        RenderForCapture();

        if (mSaveAsImages)
        {
            // Redirect the writer from its in-memory buffer to a per-step file.
            auto p_writer = MakePngWriter();
            if (!mOutputFilePath.empty())
            {
                p_writer->SetWriteToMemory(0);
                p_writer->SetFileName((mOutputFilePath + "_" + boost::lexical_cast<std::string>(timeStep) + ".png").c_str());
                p_writer->Write();
            }
        }

        if (mSaveAsAnimation)
        {
            mAnimationWriter->Write();
        }
    }

    // Interactive mode: show the freshly-rebuilt scene in the on-screen window.
    if (mIsInteractive)
    {
        mpRenderWindow->SetOffScreenRendering(0);
        mpRenderWindow->Render();
    }
}

template <unsigned DIM>
void VtkScene<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation)
{
    mpCellPopulationGenerator->SetCellPopulation(pCellPopulation);
}

template <unsigned DIM>
void VtkScene<DIM>::SetIsInteractive(bool isInteractive)
{
    mIsInteractive = isInteractive;
}

template <unsigned DIM>
void VtkScene<DIM>::SetOutputFilePath(const std::string& rPath)
{
    mOutputFilePath = rPath;
}

template <unsigned DIM>
void VtkScene<DIM>::SetSaveAsAnimation(bool saveAsAnimation)
{
    mSaveAsAnimation = saveAsAnimation;
}

template <unsigned DIM>
void VtkScene<DIM>::SetSaveAsImages(bool saveAsImages)
{
    mSaveAsImages = saveAsImages;
}

template <unsigned DIM>
void VtkScene<DIM>::Start()
{
    // Build the scene before the setup render below, so the window-to-image
    // filter and animation writer are initialised from a populated frame.
    UpdateScene();
    if constexpr (DIM == 3)
    {
        mpRenderer->GetActiveCamera()->Azimuth(DEFAULT_3D_AZIMUTH_DEGREES);
    }

    // The window-to-image filter feeds every capture path (the in-memory PNG
    // buffer as well as saved images/animation), so wire it to the render
    // window unconditionally, not only when saving output. Read the back buffer,
    // which is where off-screen rendering draws.
    mWindowToImageFilter->SetInput(mpRenderWindow);
    mWindowToImageFilter->SetReadFrontBuffer(0);

    // Render off-screen for non-interactive display and whenever writing output.
    if (!mIsInteractive || mSaveAsImages || mSaveAsAnimation)
    {
        mpRenderWindow->SetOffScreenRendering(1);

        // Prime the window-to-image filter from a populated, realized frame. Every
        // off-screen capture path needs this, including GetSceneAsCharBuffer()
        // (vtk_show): on VTK 7 capturing an unprimed off-screen window segfaults.
        mpRenderWindow->Render();
        mWindowToImageFilter->Update();

        if (mSaveAsAnimation)
        {
            mAnimationWriter->SetInputConnection(mWindowToImageFilter->GetOutputPort());
            mAnimationWriter->SetFileName((mOutputFilePath + ".ogg").c_str());
            mAnimationWriter->SetRate(ANIMATION_FRAME_RATE);
            mAnimationWriter->Start();
            mAnimationWriterStarted = true;
        }
    }

    mHasStarted = true;

    if (mIsInteractive)
    {
        mpRenderWindowInteractor->Initialize();
    }
}

template <unsigned DIM>
void VtkScene<DIM>::StartInteractiveEventHandler()
{
    mpRenderWindowInteractor->Start();
}

template <unsigned DIM>
void VtkScene<DIM>::UpdateScene()
{
    // Rebuild the scene's actors to reflect the current cell population, and
    // reset the camera. Unlike RenderFrame(), this writes no output.
    vtkSmartPointer<vtkActor> p_actor;
    vtkSmartPointer<vtkActorCollection> p_actors = mpRenderer->GetActors();

    for (p_actors->InitTraversal(); (p_actor = p_actors->GetNextItem()) != nullptr;)
    {
        mpRenderer->RemoveActor(p_actor);
    }

    if (mpCellPopulationGenerator)
    {
        mpCellPopulationGenerator->AddActor(mpRenderer);
    }
    mpRenderer->ResetCamera();
}

template <unsigned DIM>
vtkSmartPointer<vtkPNGWriter> VtkScene<DIM>::MakePngWriter()
{
    auto p_writer = vtkSmartPointer<vtkPNGWriter>::New();
    p_writer->SetWriteToMemory(1);
    p_writer->SetInputConnection(mWindowToImageFilter->GetOutputPort());
    return p_writer;
}

template <unsigned DIM>
void VtkScene<DIM>::RenderForCapture()
{
    mpRenderWindow->SetOffScreenRendering(1);
    mpRenderWindow->Render();
    mWindowToImageFilter->Modified();
}

template class VtkScene<2>;
template class VtkScene<3>;
