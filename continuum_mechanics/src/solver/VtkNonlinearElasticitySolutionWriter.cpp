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

#include "VtkNonlinearElasticitySolutionWriter.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
void VtkNonlinearElasticitySolutionWriter<DIM>::Write()
{
    if (mpSolver->mOutputDirectory=="")
    {
        EXCEPTION("No output directory was given to the mechanics solver");
    }

#ifdef CHASTE_VTK
    VtkMeshWriter<DIM, DIM> mesh_writer(mpSolver->mOutputDirectory + "/vtk", "solution", true);

    // write the displacement
    std::vector<c_vector<double,DIM> > displacement(mpSolver->mrQuadMesh.GetNumNodes());
    std::vector<c_vector<double,DIM> >& r_spatial_solution = mpSolver->rGetSpatialSolution();
    for (unsigned i=0; i<mpSolver->mrQuadMesh.GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            displacement[i](j) = r_spatial_solution[i](j)- mpSolver->mrQuadMesh.GetNode(i)->rGetLocation()[j];
        }
    }
    mesh_writer.AddPointData("Displacement", displacement);

    // write pressures
    if (mpSolver->mCompressibilityType==INCOMPRESSIBLE)
    {
        mesh_writer.AddPointData("Pressure", mpSolver->rGetPressures());
    }

    // write the element attribute as cell data.
    std::vector<double> element_attribute;
    for (typename QuadraticMesh<DIM>::ElementIterator iter = mpSolver->mrQuadMesh.GetElementIteratorBegin();
        iter != mpSolver->mrQuadMesh.GetElementIteratorEnd();
        ++iter)
    {
        element_attribute.push_back(iter->GetAttribute());
    }
    mesh_writer.AddCellData("Attribute", element_attribute);

    // write strains if requested
    if (mWriteElementWiseStrains)
    {
        mTensorData.clear();
        mTensorData.resize(mpSolver->mrQuadMesh.GetNumElements());

        std::string name;
        switch(mElementWiseStrainType)
        {
            case DEFORMATION_GRADIENT_F:
            {
                name = "deformation_gradient_F";
                break;
            }
            case DEFORMATION_TENSOR_C:
            {
                name = "deformation_tensor_C";
                break;
            }
            case LAGRANGE_STRAIN_E:
            {
                name = "Lagrange_strain_E";
                break;
            }
            default:
            {
                NEVER_REACHED;
            }
        }

        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = mpSolver->mrQuadMesh.GetElementIteratorBegin();
             iter != mpSolver->mrQuadMesh.GetElementIteratorEnd();
             ++iter)
        {
            mpSolver->GetElementCentroidStrain(mElementWiseStrainType, *iter, mTensorData[iter->GetIndex()]);
        }

        mesh_writer.AddTensorCellData(name, mTensorData);
    }
//// Future..
//        if (mWriteNodeWiseStresses)
//        {
//            std::vector<c_matrix<double,DIM,DIM> > tensor_data;
//            // use recoverer
//            mesh_writer.AddTensorCellData("Stress_NAME_ME", tensor_data);
//        }

    // final write
    mesh_writer.WriteFilesUsingMesh(mpSolver->mrQuadMesh);
#endif // CHASTE_VTK
}

// Explicit instantiation
template class VtkNonlinearElasticitySolutionWriter<2>;
template class VtkNonlinearElasticitySolutionWriter<3>;
