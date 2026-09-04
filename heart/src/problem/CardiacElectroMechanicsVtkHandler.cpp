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


#include "CardiacElectroMechanicsVtkHandler.hpp"

template<unsigned DIM, unsigned ELEC_PROB_DIM>
CardiacElectroMechanicsVtkHandler<DIM,ELEC_PROB_DIM>::CardiacElectroMechanicsVtkHandler(AbstractNonlinearElasticitySolver<DIM>& rMechanicsSolver,
    QuadraticMesh<DIM>& rQuadMesh,
    TetrahedralMesh<DIM,DIM>& rElectricsMesh,
    ReplicatableVector& rElectricsInitialCondition,
    const std::string& rOutputDir)
   : mrMechanicsSolver(rMechanicsSolver),
     mpVtkElastictyWriter(nullptr),
     mpInterpolater(nullptr)
{

#ifdef CHASTE_VTK // Requires "sudo aptitude install libvtk5-dev" or similar
    //create an internal copy of the quadratic mesh (we will modify node locations for oputput)
    mpVtkOutputMesh = new QuadraticMesh<DIM>();
    mpVtkOutputMesh->ConstructFromMesh(rQuadMesh);
    assert(mpVtkOutputMesh->GetNumNodes() == rQuadMesh.GetNumNodes());

    //cretae the writer, using the output mesh voltage and other quantities
    mpVtkWriter =  new VtkDeformedMeshWriter<DIM>(mpVtkOutputMesh, rOutputDir + "/vtk","deformed_mechanics_mesh_",true);
    //set up voltage interpolating object
    mpInterpolater = new VoltageInterpolaterOntoMechanicsMesh<DIM>(rElectricsMesh,*mpVtkOutputMesh);
    mpVtkElastictyWriter = new VtkNonlinearElasticitySolutionWriter<DIM>(rMechanicsSolver);

    //allocate memory for the interpolated voltages (once on initialize)
    mInterpolatedVoltagesNodeWise.assign(mpVtkOutputMesh->GetNumNodes(),0.0);
    //allocate memory for displacements
    c_vector<double,DIM> no_displ;
    c_matrix<double,DIM,DIM> zero_strains;
    for (unsigned i=0u ; i < DIM; ++i)
    {
        no_displ(i) = 0.0;
        for (unsigned j=0u; j < DIM; ++j)
        {
            zero_strains(i,j) = 0.0;
        }
    }
    mDisplacements.assign(mpVtkOutputMesh->GetNumNodes(),no_displ);
    mStrains.assign(mpVtkOutputMesh->GetNumElements(),zero_strains);

    //Create initial condition vector for the voltages (passed in)
    ReplicatableVector init_cond(rElectricsMesh.GetNumNodes());
    for (unsigned i = 0u; i < init_cond.GetSize(); ++i)
    {
        //the following line assumes interleaved solution for ELEC_PROB_DIM>1 (e.g, [Vm_0, phi_e_0, Vm1, phi_e_1...])
        init_cond[i] = rElectricsInitialCondition[ELEC_PROB_DIM*i];
    }
    mpInterpolater->InterpolateOnCoarseMesh(mInterpolatedVoltagesNodeWise,init_cond);

    mpVtkWriter->SetOutputBaseFileName("deformed_mechanics_mesh_" + std::to_string(0));
    mpVtkWriter->AddPointData("V", mInterpolatedVoltagesNodeWise);
    mpVtkWriter->AddPointData("displacements", mDisplacements);
    mpVtkWriter->AddTensorCellData("deformation_gradient_F", mStrains);
    assert(mInterpolatedVoltagesNodeWise.size()==rQuadMesh.GetNumNodes());
    mpVtkWriter->WriteDeformedFiles();

#endif //CHASTE_VTK
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
CardiacElectroMechanicsVtkHandler<DIM,ELEC_PROB_DIM>::~CardiacElectroMechanicsVtkHandler()
{
#ifdef CHASTE_VTK // Requires "sudo aptitude install libvtk5-dev" or similar
    delete mpVtkWriter;
    delete mpVtkOutputMesh;
#endif
    delete mpInterpolater;
    delete mpVtkElastictyWriter;
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsVtkHandler<DIM,ELEC_PROB_DIM>::WriteSolution(unsigned counter, ReplicatableVector& rElectricsSolution)
{
#ifdef CHASTE_VTK // Requires "sudo aptitude install libvtk5-dev" or similar
    //Apply deformation solution to mechanics mesh (the one in the writer object, not the one used by the solver!)
    mpVtkWriter->ApplyDeformation(mrMechanicsSolver.rGetDeformedPosition());
    mpVtkWriter->SetOutputBaseFileName("deformed_mechanics_mesh_" + std::to_string(counter));
    mpVtkWriter->SetWriteMeshCells(false);

    //voltage
    mpInterpolater->InterpolateOnCoarseMesh(mInterpolatedVoltagesNodeWise,rElectricsSolution);//VTK output. Determine voltage on coarse nodes
    mpVtkWriter->AddPointData("V", mInterpolatedVoltagesNodeWise);

    //displacement
    mpVtkElastictyWriter->CalculateDisplacements(mDisplacements);
    mpVtkWriter->AddPointData("displacements", mDisplacements);

    //Strains
    mpVtkElastictyWriter->SetWriteElementWiseStrains(DEFORMATION_GRADIENT_F);
    std::string name;
    mpVtkElastictyWriter->CalculateStrains(name, mStrains);
    mpVtkWriter->AddTensorCellData(name,mStrains);

    //write to file
    mpVtkWriter->WriteDeformedFiles();
#endif
}


// Explicit instantiation
template class CardiacElectroMechanicsVtkHandler<2,1>;
template class CardiacElectroMechanicsVtkHandler<3,1>;
template class CardiacElectroMechanicsVtkHandler<2,2>;
template class CardiacElectroMechanicsVtkHandler<3,2>;
