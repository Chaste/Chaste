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
#include "BielmeierForce.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include "UblasCustomFunctions.hpp"
#include "VertexBasedCellPopulation.hpp"

BielmeierForce::BielmeierForce(const double springConstant, const double ExternalSurfaceTensionParameter)
        : GeneralMonolayerVertexMeshForce(),
          mEcmSpringConstant(springConstant),
          mExternalSurfaceTensionParameter(ExternalSurfaceTensionParameter)
{
    this->SetApicalParameters(0.18, 3.1);
    this->SetBasalParameters(0.18, 6.95);
    this->SetLateralParameter(0, 1);
    this->SetVolumeParameters(1000);
}

void BielmeierForce::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
    // Adding internal forces
    this->GeneralMonolayerVertexMeshForce::AddForceContribution(rCellPopulation);

    if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("BielmeierForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);
    MutableVertexMesh<3, 3>& rMesh = p_cell_population->rGetMesh();

    // Adding external forces
    for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
    {
        Node<3>* p_node = rMesh.GetNode(i);

        if (IsBasalNode(p_node))
        {
            c_vector<double, 3> result = Create_c_vector(0, 0, p_node->rGetLocation()[2]);
            result *= -1 * mEcmSpringConstant;
            p_node->AddAppliedForceContribution(result);
        }
    }

    for (unsigned i = 0; i < rMesh.GetNumFaces(); ++i)
    {
        VertexElement<2, 3>* p_face = rMesh.GetFace(i);

        if (IsApicalFace(p_face))
        {
            for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
            {
                c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, j);
                result *= mExternalSurfaceTensionParameter;
                p_face->GetNode(j)->AddAppliedForceContribution(result);
            }
        }
    }
}

void BielmeierForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EcmSpringConstant>" << mEcmSpringConstant << "</EcmSpringConstant>\n";
    *rParamsFile << "\t\t\t<ExternalSurfaceTensionParameter>" << mExternalSurfaceTensionParameter << "</ExternalSurfaceTensionParameter>\n";
    GeneralMonolayerVertexMeshForce::OutputForceParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BielmeierForce)
