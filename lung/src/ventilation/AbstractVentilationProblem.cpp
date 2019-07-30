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

#include "AbstractVentilationProblem.hpp"
#include "MathsCustomFunctions.hpp"
#include "Warnings.hpp"
#include "AirwayTreeWalker.hpp"
#include "AirwayPropertiesCalculator.hpp"

AbstractVentilationProblem::AbstractVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex)
    : mOutletNodeIndex(rootIndex),
      mViscosity(1.92e-5),
      mDensity(1.15),
      mLengthScaling(1.0),
      mDynamicResistance(false),
      mPerGenerationDynamicResistance(false),
      mRadiusOnEdge(false),
      mNodesInGraphOrder(true)
{
    TrianglesMeshReader<1,3> mesh_reader(rMeshDirFilePath);
    mMesh.ConstructFromMeshReader(mesh_reader);

    if (mMesh.GetNode(mOutletNodeIndex)->IsBoundaryNode() == false)
    {
        EXCEPTION("Outlet node is not a boundary node");
    }
    Initialise();
}
void
AbstractVentilationProblem::Initialise()
{
    /// We might want to make the tree walker a member of this class
    AirwayTreeWalker walker(mMesh, mOutletNodeIndex);
    mNodesInGraphOrder = walker.GetNodesAreGraphOrdered();

    // Reset edge attributes
    bool intermediate_nodes = false;
    for (AbstractTetrahedralMesh<1,3>::ElementIterator iter = mMesh.GetElementIteratorBegin();
         iter != mMesh.GetElementIteratorEnd();
         ++iter)
    {
        // Add a place holder for radius if necessary
        if (iter->GetNumElementAttributes() == 0u)
        {
            // Radius on edge is not going to be used but we'll put in a placeholder
            iter->AddElementAttribute(DOUBLE_UNSET);
            assert( iter->rGetElementAttributes()[RADIUS] == DOUBLE_UNSET);
        }

        // Add a length attribute
        assert (iter->GetNumElementAttributes() == 1u); // Only radius used so far
        // Set the segment length to be the actual length
        c_vector<double, 3> dummy;
        double length;
        unsigned element_index = iter->GetIndex();
        mMesh.GetWeightedDirectionForElement(element_index, dummy, length);
        iter->AddElementAttribute(length);
        assert( iter->rGetElementAttributes()[SEGMENT_LENGTH] == length);

        // Check for intermediate nodes
        if (walker.GetChildElementIndices(&*iter).size() == 1u)
        {
            intermediate_nodes = true;
        }
    }
    if (intermediate_nodes)
    {
        ///\todo #2300 There's redundancy here because there are now two walkers
        AirwayPropertiesCalculator properties_calculator(mMesh, mOutletNodeIndex);
        std::vector<AirwayBranch*> branches = properties_calculator.GetBranches();
        for (std::vector<AirwayBranch*>::iterator branch_it=branches.begin(); branch_it != branches.end(); branch_it++)
        {
            std::list<Element<1,3>* > branch_elements = (*branch_it)->GetElements();
            if (branch_elements.size() > 1u)
            {
                double branch_length = (*branch_it)->GetLength();

                for (std::list<Element<1,3>* >::const_iterator element_iterator=branch_elements.begin();
                        element_iterator != branch_elements.end(); ++element_iterator)
                {
                    (*element_iterator)->rGetElementAttributes()[SEGMENT_LENGTH] = branch_length;
                }
            }
//            else
//            {
//                assert( (*branch_it)->GetLength() == branch_elements.front()->rGetElementAttributes()[SEGMENT_LENGTH]);
//            }
        }
    }
}

void
AbstractVentilationProblem::SetMeshInMilliMetres()
{
    mLengthScaling = 1e-3;
}

void AbstractVentilationProblem::SetRadiusOnEdge(bool isOnEdges)
{
    mRadiusOnEdge = isOnEdges;
}

TetrahedralMesh<1, 3>&
AbstractVentilationProblem::rGetMesh()
{
    return mMesh;
}

double AbstractVentilationProblem::CalculateResistance(Element<1,3>& rElement, bool usePedley, double flux)
{
    //Resistance is based on radius, length and viscosity
    double radius = 0.0;
    if (mRadiusOnEdge)
    {
        radius = rElement.GetAttribute();
    }
    else
    {
        radius = ( rElement.GetNode(0)->rGetNodeAttributes()[0] + rElement.GetNode(1)->rGetNodeAttributes()[0]) / 2.0;
    }

    c_vector<double, 3> dummy;
    double length;
    mMesh.GetWeightedDirectionForElement(rElement.GetIndex(), dummy, length);

    radius *= mLengthScaling;
    length *= mLengthScaling;

    double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
    if (usePedley)
    {
        /* Pedley et al. 1970
         * http://dx.doi.org/10.1016/0034-5687(70)90094-0
         * also Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
         * Standard Poiseuille equation is similar to Pedley's modified Eq1. and matches Swan Eq5.
         * R_p = 128*mu*L/(pi*d^4) = 8*mu*L/(pi*r^4)
         *
         * Pedley Eq 2 and Swan Eq4:
         *  Z = C/(4*sqrt(2)) * sqrt(Re*d/l) = (C/4)*sqrt(Re*r/l)
         * Pedley suggests that C = 1.85
         * R_r = Z*R_p
         *
         * Reynold's number in a pipe is
         * Re = rho*v*d/mu (where d is a characteristic length scale - diameter of pipe)
         * since flux = v*area
         * Re = Q * d/(mu*area) = 2*rho*Q/(mu*pi*r) ... (see Swan p 224)
         *
         *
         * The upshot of this calculation is that the resistance is scaled with sqrt(Q)
         */
        std::vector<double>& r_attributes = rElement.rGetElementAttributes();
        double segment_length = r_attributes [SEGMENT_LENGTH] * mLengthScaling;
        double reynolds_number = fabs( 2.0 * mDensity * flux / (mViscosity * M_PI * radius) );
        double pedley_c = 1.85;
        if (mPerGenerationDynamicResistance)
        {
            pedley_c = r_attributes[PEDLEY_CORRECTION];
        }
        // assert(mDynamicResistance);

        double z = (pedley_c/4.0) * sqrt(reynolds_number * radius / segment_length);

        // Pedley's method will only increase the resistance
        if (z > 1.0)
        {
            resistance *= z;
        }
        //if (rElement.GetIndex() == 0) PRINT_3_VARIABLES(reynolds_number, z, resistance);
    }
    return resistance;
}
void AbstractVentilationProblem::SetConstantInflowPressures(double pressure)
{
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter =mMesh.GetBoundaryNodeIteratorBegin();
          iter != mMesh.GetBoundaryNodeIteratorEnd();
          ++iter)
     {
         if ((*iter)->GetIndex() != mOutletNodeIndex)
         {
             //Boundary conditions at each boundary/leaf node
             SetPressureAtBoundaryNode(*(*iter), pressure);
         }
     }
}

void AbstractVentilationProblem::SetConstantInflowFluxes(double flux)
{
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter =mMesh.GetBoundaryNodeIteratorBegin();
          iter != mMesh.GetBoundaryNodeIteratorEnd();
          ++iter)
     {
         if ((*iter)->GetIndex() != mOutletNodeIndex)
         {
             SetFluxAtBoundaryNode(*(*iter), flux);
         }
     }
}
void AbstractVentilationProblem::SetOutflowPressure(double pressure)
{
    SetPressureAtBoundaryNode(*(mMesh.GetNode(mOutletNodeIndex)), pressure);
}

void AbstractVentilationProblem::SetOutflowFlux(double flux)
{
    SetFluxAtBoundaryNode(*(mMesh.GetNode(mOutletNodeIndex)), flux);
}


void AbstractVentilationProblem::SolveOverTime(TimeStepper& rTimeStepper,
        void (*pBoundaryConditionFunction)(AbstractVentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&),
        const std::string& rDirName, const std::string& rFileBaseName)
{
#ifdef CHASTE_VTK
    VtkMeshWriter<1, 3> vtk_writer(rDirName, rFileBaseName, false);
#endif

    bool first_step=true;
    while (!rTimeStepper.IsTimeAtEnd())
    {
        if (first_step)
        {
            // Do a solve at t=0 before advancing time.
            first_step = false;
        }
        else
        {
            rTimeStepper.AdvanceOneTimeStep();
        }
        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
                 iter != mMesh.GetBoundaryNodeIteratorEnd();
                 ++iter )
        {
            if ((*iter)->GetIndex() != mOutletNodeIndex)
            {
                //Boundary conditions at each boundary/leaf node
                pBoundaryConditionFunction(this, rTimeStepper, *(*iter));
            }
        }

        // Regular solve
        Solve();

        std::ostringstream suffix_name;
        suffix_name <<  "_" << std::setw(6) << std::setfill('0') << rTimeStepper.GetTotalTimeStepsTaken();
#ifdef CHASTE_VTK
        AddDataToVtk(vtk_writer, suffix_name.str());
#endif
    }

#ifdef CHASTE_VTK
    vtk_writer.WriteFilesUsingMesh(mMesh);
#endif
}

void AbstractVentilationProblem::SolveProblemFromFile(const std::string& rInFilePath, const std::string& rOutFileDir, const std::string& rOutFileName)
{
    std::ifstream file(FileFinder(rInFilePath).GetAbsolutePath().c_str(), std::ios::binary);
    if (!file.is_open())
    {
        EXCEPTION("Could not open file "+rInFilePath);
    }
    std::string key, unit;
    double value;
    while (!file.eof())
    {
        file >> key >> value >> unit;
        if (file.fail())
        {
            break;
        }
        if (key == "RHO_AIR")
        {
            SetDensity(value);
        }
        else if (key == "MU_AIR")
        {
            SetViscosity(value);
        }
        else if (key == "PRESSURE_OUT")
        {
            SetOutflowPressure(value);
        }
        else if (key == "PRESSURE_IN")
        {
            SetConstantInflowPressures(value);
        }
        else
        {
            WARNING("The key "+ key+ " is not recognised yet");
        }
    }
    Solve();
#ifdef CHASTE_VTK
    WriteVtk(rOutFileDir, rOutFileName);
#endif
}

#ifdef CHASTE_VTK
void AbstractVentilationProblem::WriteVtk(const std::string& rDirName, const std::string& rFileBaseName)
{
    VtkMeshWriter<1, 3> vtk_writer(rDirName, rFileBaseName, false);
    AddDataToVtk(vtk_writer, "");
    vtk_writer.WriteFilesUsingMesh(mMesh);

}

void AbstractVentilationProblem::AddDataToVtk(VtkMeshWriter<1, 3>& rVtkWriter,
        const std::string& rSuffix)
{
    std::vector<double> pressures;
    std::vector<double> fluxes;
    GetSolutionAsFluxesAndPressures(fluxes, pressures);
    rVtkWriter.AddCellData("Flux"+rSuffix, fluxes);
    rVtkWriter.AddPointData("Pressure"+rSuffix, pressures);
}
#endif // CHASTE_VTK


void AbstractVentilationProblem::SetPerGenerationDynamicResistance()
{
    mDynamicResistance = true;
    mPerGenerationDynamicResistance = true;

    AirwayTreeWalker walker(mMesh, mOutletNodeIndex);
    /*
     * According to van Ertbruggen 2005 DOI: 10.1152/japplphysiol.00795.2004 equation 3,
     * Pedley's original correction used a value
     * gamma =  0.327
     * which is equivalent to C/(4*sqrt(2)) = 1.85/(4*sqrt(2)) in Pedley's notation.
     * van Ertbruggen's generational-based gamma's vary between 0.162 and 0.566 (Fig. 8B).
     * Ismail et al. 2013 DOI: 10.1002/cnm.2577 table 2 tabulates:
     * Generation: 0     1     2     3     4     5     6     7     >7
     * Gamma:      0.162 0.239 0.244 0.295 0.175 0.303 0.356 0.566 0.327
     */
    double per_generation_pedley[9] =
            {0.162, 0.239, 0.244, 0.295, 0.175, 0.303, 0.356, 0.566, 0.327};

    // Convert from gamma to C
    for (unsigned i=0; i<9; i++)
    {
        per_generation_pedley[i] *= 4.0*sqrt(2.0);
    }
    for (AbstractTetrahedralMesh<1,3>::ElementIterator iter = mMesh.GetElementIteratorBegin();
         iter != mMesh.GetElementIteratorEnd();
         ++iter)
    {
        unsigned gen = walker.GetElementGeneration((*iter).GetIndex());
        // Add an attribute for Pedley and check it's in the correct place
        double pedley_c = per_generation_pedley[8]; // Lowest possible by default for deep branches
        if (gen < 8)
        {
            pedley_c = per_generation_pedley[gen];
        }
        iter->AddElementAttribute(pedley_c);
        assert( iter->rGetElementAttributes()[PEDLEY_CORRECTION] == pedley_c);
    }
}
