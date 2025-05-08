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

#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr/make_shared.hpp>

#include <vtkActor.h>
#include <vtkActorCollection.h>
#include <vtkCamera.h>
#include <vtkCell.h>
#include <vtkConvexPointSet.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkExtractEdges.h>
#include <vtkFeatureEdges.h>
#include <vtkGeometryFilter.h>
#include <vtkGlyph2D.h>
#include <vtkGlyph3D.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLine.h>
#include <vtkNamedColors.h>
#include <vtkObjectFactory.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkWindowToImageFilter.h>

#include "Exception.hpp"
#include "UblasIncludes.hpp"
#include "UblasVectorInclude.hpp"

#include "VtkScene.hpp"

// For some reason an explicit interactor style is needed capture mouse events
class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static customMouseInteractorStyle* New();
    vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    virtual void OnLeftButtonDown()
    {
        // Forward events
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    virtual void OnMiddleButtonDown()
    {
        // Forward events
        vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }

    virtual void OnRightButtonDown()
    {
        // Forward events
        vtkInteractorStyleTrackballCamera::OnRightButtonDown();
    }
};

vtkStandardNewMacro(customMouseInteractorStyle);

template <unsigned DIM>
VtkScene<DIM>::VtkScene()
        : mpRenderer(vtkSmartPointer<vtkRenderer>::New()),
          mpRenderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
          mpRenderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
          mOutputFilePath(),
          mAnimationWriter(vtkSmartPointer<vtkOggTheoraWriter>::New()),
          mWindowToImageFilter(vtkSmartPointer<vtkWindowToImageFilter>::New()),
          mIsInteractive(false),
          mSaveAsAnimation(false),
          mSaveAsImages(false),
          mHasStarted(false),
          mAddAnnotations(false),
          mOutputFrequency(1),
          mpCellPopulationGenerator(boost::make_shared<CellPopulationPyChasteActorGenerator<DIM> >())
{
    mpRenderer->SetBackground(1.0, 1.0, 1.0);
    mpRenderWindow->AddRenderer(mpRenderer);
    mpRenderWindow->SetSize(800.0, 600.0);
    mpRenderWindowInteractor->SetRenderWindow(mpRenderWindow);

    auto style = vtkSmartPointer<customMouseInteractorStyle>::New();
    mpRenderWindowInteractor->SetInteractorStyle(style);
}

template <unsigned DIM>
VtkScene<DIM>::~VtkScene()
{
}

template <unsigned DIM>
boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > VtkScene<DIM>::GetCellPopulationActorGenerator()
{
    return mpCellPopulationGenerator;
}

template <unsigned DIM>
void VtkScene<DIM>::SetIsInteractive(bool isInteractive)
{
    mIsInteractive = isInteractive;
}

template <unsigned DIM>
vtkSmartPointer<vtkRenderer> VtkScene<DIM>::GetRenderer()
{
    return mpRenderer;
}

template <unsigned DIM>
void VtkScene<DIM>::SetSaveAsImages(bool saveAsImages)
{
    mSaveAsImages = saveAsImages;
}

template <unsigned DIM>
vtkSmartPointer<vtkUnsignedCharArray> VtkScene<DIM>::GetSceneAsCharBuffer()
{
    ResetRenderer(0);

    mpRenderWindow->SetOffScreenRendering(1);
    mpRenderWindow->Render();
    mWindowToImageFilter->Modified();

    auto p_writer = vtkSmartPointer<vtkPNGWriter>::New();
    p_writer->SetWriteToMemory(1);
    p_writer->SetInputConnection(mWindowToImageFilter->GetOutputPort());
    p_writer->Write();

    return p_writer->GetResult();
}

template <unsigned DIM>
void VtkScene<DIM>::ResetRenderer(unsigned time_step)
{
    if (!mHasStarted)
    {
        Start();
    }

    vtkSmartPointer<vtkActor> p_actor;
    vtkSmartPointer<vtkActorCollection> p_actors = mpRenderer->GetActors();

    for (p_actors->InitTraversal(); (p_actor = p_actors->GetNextItem()) != NULL;)
    {
        mpRenderer->RemoveActor(p_actor);
    }

    if (mpCellPopulationGenerator)
    {
        mpCellPopulationGenerator->AddActor(mpRenderer);
    }
    mpRenderer->ResetCamera();

    if (mSaveAsImages)
    {
        mpRenderWindow->SetOffScreenRendering(1);
        mpRenderWindow->Render();
        mWindowToImageFilter->Modified();

        auto p_writer = vtkSmartPointer<vtkPNGWriter>::New();
        p_writer->SetWriteToMemory(1);
        p_writer->SetInputConnection(mWindowToImageFilter->GetOutputPort());

        if (!mOutputFilePath.empty())
        {
            p_writer->SetWriteToMemory(0);
            p_writer->SetFileName((mOutputFilePath + "_" + boost::lexical_cast<std::string>(time_step) + ".png").c_str());
            p_writer->Write();
        }
    }

    if (mSaveAsAnimation)
    {
        if (!mSaveAsImages)
        {
            mpRenderWindow->SetOffScreenRendering(1);
            mpRenderWindow->Render();
            mWindowToImageFilter->Modified();
        }
        mAnimationWriter->Write();
    }

    if (mIsInteractive)
    {
        mpRenderWindow->SetOffScreenRendering(0);
        mpRenderWindow->Render();
    }
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
void VtkScene<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation)
{
    mpCellPopulationGenerator->SetCellPopulation(pCellPopulation);
}

template <unsigned DIM>
void VtkScene<DIM>::End()
{
    if (mSaveAsAnimation and mHasStarted)
    {
        mAnimationWriter->End();
    }
}

template <unsigned DIM>
void VtkScene<DIM>::Start()
{
    mpRenderer->ResetCamera();
    if (DIM == 3)
    {
        mpRenderer->GetActiveCamera()->Azimuth(45.0);
    }

    if (!mIsInteractive)
    {
        mpRenderWindow->SetOffScreenRendering(1);
    }

    if (mSaveAsImages || mSaveAsAnimation)
    {
        mpRenderWindow->SetOffScreenRendering(1);
        mpRenderWindow->Render();
        mWindowToImageFilter->SetInput(mpRenderWindow);
        mWindowToImageFilter->Update();
    }

    if (mSaveAsAnimation)
    {
        mAnimationWriter->SetInputConnection(mWindowToImageFilter->GetOutputPort());
        mAnimationWriter->SetFileName((mOutputFilePath + ".ogg").c_str());
        mAnimationWriter->SetRate(1.0);
        mAnimationWriter->Start();
    }

    mHasStarted = true;
    ResetRenderer();

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

template class VtkScene<2>;
template class VtkScene<3>;
