/*
 * BielmeierForce.cpp
 *
 *  Created on: 23 Mar 2017
 *      Author: Weijie
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

    if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("GeneralMonolayerVertexMeshForce is to be used with a VertexBasedCellPopulation only");
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
