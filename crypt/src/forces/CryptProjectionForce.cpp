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

#include "CryptProjectionForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "WntConcentration.hpp"
#include "IsNan.hpp"
#include "StemCellProliferativeType.hpp"

CryptProjectionForce::CryptProjectionForce()
    : GeneralisedLinearSpringForce<2>(),
      mIncludeWntChemotaxis(false),
      mWntChemotaxisStrength(100.0)
{
    mA = WntConcentration<2>::Instance()->GetCryptProjectionParameterA();
    mB = WntConcentration<2>::Instance()->GetCryptProjectionParameterB();
}

CryptProjectionForce::~CryptProjectionForce()
{
}

void CryptProjectionForce::UpdateNode3dLocationMap(AbstractCellPopulation<2>& rCellPopulation)
{
    mNode3dLocationMap.clear();

    c_vector<double, 2> node_location_2d;
    c_vector<double, 3> node_location_3d;

    // Only consider nodes corresponding to real cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get node index
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Get 3D location
        node_location_2d = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        node_location_3d[0] = node_location_2d[0];
        node_location_3d[1] = node_location_2d[1];
        node_location_3d[2] = CalculateCryptSurfaceHeightAtPoint(node_location_2d);

        // Add to map
        mNode3dLocationMap[node_index] = node_location_3d;
    }
}

double CryptProjectionForce::GetA() const
{
    return mA;
}

double CryptProjectionForce::GetB() const
{
    return mB;
}

void CryptProjectionForce::SetWntChemotaxisStrength(double wntChemotaxisStrength)
{
    assert(wntChemotaxisStrength >= 0.0);
    mWntChemotaxisStrength = wntChemotaxisStrength;
}
double CryptProjectionForce::GetWntChemotaxisStrength()
{
    return mWntChemotaxisStrength;
}

void CryptProjectionForce::SetWntChemotaxis(bool includeWntChemotaxis)
{
    mIncludeWntChemotaxis = includeWntChemotaxis;
}

double CryptProjectionForce::CalculateCryptSurfaceHeightAtPoint(const c_vector<double,2>& rNodeLocation)
{
    return mA*pow(norm_2(rNodeLocation), mB); // =z_coord;
}

double CryptProjectionForce::CalculateCryptSurfaceDerivativeAtPoint(const c_vector<double,2>& rNodeLocation)
{
    return mA*mB*pow(norm_2(rNodeLocation), (mB - 1.0));
}

c_vector<double,2> CryptProjectionForce::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<2>& rCellPopulation)
{
    MeshBasedCellPopulation<2>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations in 2D
    c_vector<double,2> node_a_location_2d = rCellPopulation.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double,2> node_b_location_2d = rCellPopulation.GetNode(nodeBGlobalIndex)->rGetLocation();

    // "Get the unit vector parallel to the line joining the two nodes" [GeneralisedLinearSpringForce]

    // Create a unit vector in the direction of the 3D spring
    c_vector<double,3> unit_difference = mNode3dLocationMap[nodeBGlobalIndex] - mNode3dLocationMap[nodeAGlobalIndex];

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    // If mUseCutOffLength has been set, then there is zero force between
    // two nodes located a distance apart greater than mUseCutOffLength
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= mMechanicsCutOffLength)
        {
            // Return zero (2D projected) force
            return zero_vector<double>(2);
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until mMeinekeSpringGrowthDuration hour after division.
     */
    if (ageA < mMeinekeSpringGrowthDuration && ageB < mMeinekeSpringGrowthDuration)
    {
        /*
         * The spring rest length increases from a predefined small parameter
         * to a normal rest length of 1.0, over a period of one hour.
         */
        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);
        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            double lambda = mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (1.0 - lambda) * ageA/mMeinekeSpringGrowthDuration;
        }
        if (ageA+SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration)
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;

    // Assert that the rest length does not exceed 1
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = true;

    if (distance_between_nodes - rest_length > 0)
    {
        is_closer_than_rest_length = false;
    }

    /*
     * Although in this class the 'spring constant' is a constant parameter, in
     * subclasses it can depend on properties of each of the cells.
     */
    double multiplication_factor = 1.0;
    multiplication_factor *= VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);

    // Calculate the 3D force between the two points
    c_vector<double,3> force_between_nodes = multiplication_factor * this->GetMeinekeSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);

    // Calculate an outward normal unit vector to the tangent plane of the crypt surface at the 3D point corresponding to node B
    c_vector<double,3> outward_normal_unit_vector;

    double dfdr = CalculateCryptSurfaceDerivativeAtPoint(node_b_location_2d);
    double theta_B = atan2(node_b_location_2d[1], node_b_location_2d[0]); // use atan2 to determine the quadrant
    double normalization_factor = sqrt(1 + dfdr*dfdr);

    outward_normal_unit_vector[0] = dfdr*cos(theta_B)/normalization_factor;
    outward_normal_unit_vector[1] = dfdr*sin(theta_B)/normalization_factor;
    outward_normal_unit_vector[2] = -1.0/normalization_factor;

    // Calculate the projection of the force onto the plane z=0
    c_vector<double,2> projected_force_between_nodes_2d;
    double force_dot_normal = inner_prod(force_between_nodes, outward_normal_unit_vector);

    for (unsigned i=0; i<2; i++)
    {
        projected_force_between_nodes_2d[i] = force_between_nodes[i]
                                              - force_dot_normal*outward_normal_unit_vector[i]
                                              + force_dot_normal*outward_normal_unit_vector[2];
    }

    return projected_force_between_nodes_2d;
}

void CryptProjectionForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
    // First work out the 3D location of each cell
    UpdateNode3dLocationMap(rCellPopulation);

    // Throw an exception message if not using a MeshBasedCellPopulation
    if (dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("CryptProjectionForce is to be used with a subclass of MeshBasedCellPopulation only");
    }

    MeshBasedCellPopulation<2>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
         spring_iterator != p_static_cast_cell_population->SpringsEnd();
         ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, 2> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);
        c_vector<double, 2> negative_force = -1.0 * force;
        spring_iterator.GetNodeB()->AddAppliedForceContribution(negative_force);
        spring_iterator.GetNodeA()->AddAppliedForceContribution(force);
    }

    if (mIncludeWntChemotaxis)
    {
        assert(WntConcentration<2>::Instance()->IsWntSetUp());

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                c_vector<double, 2> wnt_chemotactic_force = mWntChemotaxisStrength*WntConcentration<2>::Instance()->GetWntGradient(*cell_iter);
                unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                rCellPopulation.GetNode(index)->AddAppliedForceContribution(wnt_chemotactic_force);
            }
        }
    }
}

void CryptProjectionForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<A>" << mA << "</A>\n";
    *rParamsFile << "\t\t\t<B>" << mB << "</B>\n";
    *rParamsFile << "\t\t\t<IncludeWntChemotaxis>" << mIncludeWntChemotaxis << "</IncludeWntChemotaxis>\n";
    *rParamsFile << "\t\t\t<WntChemotaxisStrength>" << mWntChemotaxisStrength << "</WntChemotaxisStrength>\n";

    // Call method on direct parent class
    GeneralisedLinearSpringForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptProjectionForce)
