/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "DiscreteSystemForceCalculator.hpp"

DiscreteSystemForceCalculator::DiscreteSystemForceCalculator(MeshBasedTissue<2>& rTissue,
                                                             std::vector<AbstractTwoBodyInteractionForce<2>*> forceCollection)
    : mrTissue(rTissue),
      mForceCollection(forceCollection),
      mEpsilon(0.01)
{
}


std::vector< std::vector<double> > DiscreteSystemForceCalculator::CalculateExtremalNormalForces()
{
    unsigned num_nodes = mrTissue.GetNumNodes();

    std::vector< std::vector<double> > extremal_normal_forces;
    std::vector<double> minimum_normal_forces(num_nodes);
    std::vector<double> maximum_normal_forces(num_nodes);

    for (unsigned i=0; i<num_nodes; i++)
    {
        std::vector<double> sampling_angles = GetSamplingAngles(i);
        std::vector<double> extremal_angles = GetExtremalAngles(i, sampling_angles);

        double minimum_normal_force_for_node_i = DBL_MAX;
        double maximum_normal_force_for_node_i = -DBL_MAX;

        for (unsigned j=0; j<extremal_angles.size(); j++)
        {
            double current_normal_force = CalculateFtAndFn(i, extremal_angles[j])[1];

            if (current_normal_force > maximum_normal_force_for_node_i)
            {
                maximum_normal_force_for_node_i = current_normal_force;
            }
            if (current_normal_force < minimum_normal_force_for_node_i)
            {
                minimum_normal_force_for_node_i = current_normal_force;
            }
        }

        assert( minimum_normal_force_for_node_i <= maximum_normal_force_for_node_i);

        minimum_normal_forces[i] = minimum_normal_force_for_node_i;
        maximum_normal_forces[i] = maximum_normal_force_for_node_i;
    }

    extremal_normal_forces.push_back(minimum_normal_forces);
    extremal_normal_forces.push_back(maximum_normal_forces);

    return extremal_normal_forces;
}


void DiscreteSystemForceCalculator::WriteResultsToFile(std::string simulationOutputDirectory)
{
    double time = SimulationTime::Instance()->GetTime();
    std::ostringstream time_string;
    time_string << time;
    std::string results_directory = simulationOutputDirectory + "/results_from_time_" + time_string.str();

    OutputFileHandler output_file_handler2(results_directory+"/", false);
    mpVizStressResultsFile = output_file_handler2.OpenOutputFile("results.vizstress");

    (*mpVizStressResultsFile) <<  time << "\t";

    double global_index;
    double x;
    double y;
    double minimum;
    double maximum;

    TetrahedralMesh<2,2>& r_mesh = mrTissue.rGetMesh();

    std::vector< std::vector<double> > extremal_normal_forces = CalculateExtremalNormalForces();

    for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
    {
        global_index = (double) i;

        x = r_mesh.GetNode(i)->rGetLocation()[0];
        y = r_mesh.GetNode(i)->rGetLocation()[1];

        minimum = extremal_normal_forces[0][i];
        maximum = extremal_normal_forces[1][i];

        (*mpVizStressResultsFile) << global_index << " " << x << " " << y << " " << minimum << " " << maximum << " ";
    }

    (*mpVizStressResultsFile) << "\n";
    mpVizStressResultsFile->close();
}


std::set<unsigned> DiscreteSystemForceCalculator::GetNeighbouringNodeIndices(unsigned index)
{
    TetrahedralMesh<2,2>& r_mesh = mrTissue.rGetMesh();

    Node<2>* p_node = r_mesh.GetNode(index);

    std::set<unsigned> neighbouring_node_indices;

    for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
         it != p_node->ContainingElementsEnd();
         ++it)
    {
        Element<2,2>* p_element = r_mesh.GetElement(*it);
        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(i);
            if (node_index!=index)
            {
                neighbouring_node_indices.insert(node_index);
            }
        }
    }
    return neighbouring_node_indices;
}


std::vector<double> DiscreteSystemForceCalculator::CalculateFtAndFn(unsigned index, double theta)
{
    TetrahedralMesh<2,2>& r_mesh = mrTissue.rGetMesh();

    std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);

    double tangential_force = 0.0;
    double normal_force = 0.0;
    double alpha;

    c_vector<double,2> unit_vec_between_nodes(2);

    for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
         iter != neighbouring_node_indices.end();
         ++iter)
    {
        // The method GetAngleBetweenNodes() returns an angle in the range (-pi,pi]
        alpha = r_mesh.GetAngleBetweenNodes(index, *iter);

        assert(alpha <= M_PI);
        assert(alpha > -M_PI);

        if ( sin(alpha-theta) > DBL_EPSILON )
        {
            // Initialise a zero force vector between neighbouring nodes
            c_vector<double,2> force_between_nodes = zero_vector<double>(2);

            // Iterate over vector of forces present and add up forces between nodes
            for (std::vector<AbstractTwoBodyInteractionForce<2>*>::iterator force_iter = mForceCollection.begin();
                 force_iter != mForceCollection.end();
                 ++force_iter)
            {
               force_between_nodes += (*force_iter)->CalculateForceBetweenNodes(index, *iter, mrTissue);
            }

            unit_vec_between_nodes[0] = cos(alpha);
            unit_vec_between_nodes[1] = sin(alpha);

            double plusminus_norm_force = inner_prod(force_between_nodes,unit_vec_between_nodes);
            tangential_force += plusminus_norm_force * cos(alpha-theta);
            normal_force += plusminus_norm_force * sin(alpha-theta);
        }
    }

    std::vector<double> ret(2);
    ret[0] = tangential_force;
    ret[1] = normal_force;

    return ret;
}


std::vector<double> DiscreteSystemForceCalculator::GetSamplingAngles(unsigned index)
{
    TetrahedralMesh<2,2>& r_mesh = mrTissue.rGetMesh();

    std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);

    std::vector<double> sampling_angles(4*neighbouring_node_indices.size());

    unsigned i=0;

    for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
         iter != neighbouring_node_indices.end();
         ++iter)
    {
        // The method GetAngleBetweenNodes() returns an angle in the range (-pi,pi]
        double alpha = r_mesh.GetAngleBetweenNodes(index, *iter);

        double alpha_minus_epsilon = alpha - mEpsilon;
        double alpha_plus_epsilon = alpha + mEpsilon;
        double alpha_plus_pi_minus_epsilon = alpha + M_PI - mEpsilon;
        double alpha_plus_pi_plus_epsilon = alpha + M_PI + mEpsilon;

        // Calculate sampling angles in the range (-pi,pi]

        #define COVERAGE_IGNORE
        if (alpha_minus_epsilon <= -M_PI)
        {
            alpha_minus_epsilon += 2*M_PI;
        }
        #undef COVERAGE_IGNORE
        sampling_angles[i] = alpha_minus_epsilon;

        assert(sampling_angles[i] <= M_PI);
        assert(sampling_angles[i] > -M_PI);
        i++;

        if (alpha_plus_epsilon > M_PI)
        {
            alpha_plus_epsilon -= 2*M_PI;
        }
        sampling_angles[i] = alpha_plus_epsilon;

        assert(sampling_angles[i] <= M_PI);
        assert(sampling_angles[i] > -M_PI);
        i++;

        if (alpha_plus_pi_minus_epsilon > M_PI)
        {
            alpha_plus_pi_minus_epsilon -= 2*M_PI;
        }
        sampling_angles[i] = alpha_plus_pi_minus_epsilon;

        assert(sampling_angles[i] <= M_PI);
        assert(sampling_angles[i] > -M_PI);
        i++;

        if (alpha_plus_pi_plus_epsilon > M_PI)
        {
            alpha_plus_pi_plus_epsilon -= 2*M_PI;
        }
        sampling_angles[i] = alpha_plus_pi_plus_epsilon;

        assert(sampling_angles[i] <= M_PI);
        assert(sampling_angles[i] > -M_PI);
        i++;
    }

    sort(sampling_angles.begin(), sampling_angles.end());
    return sampling_angles;
}


double DiscreteSystemForceCalculator::GetLocalExtremum(unsigned index, double angle1, double angle2)
{
    // We always pass in angle1 and angle2 such that angle1<angle2,
    // but note that angle1 may be <M_PI
    assert(angle1 < angle2);

    double tolerance = 1e-5;
    unsigned counter = 0;

    double previous_angle;
    double current_error;
    double current_angle = angle1;

    current_error = angle2 - angle1;
    std::vector<double> current_ft_and_fn(2);

    while (current_error > tolerance)
    {
        previous_angle = current_angle;
        current_ft_and_fn = CalculateFtAndFn(index, current_angle);
        current_angle -= current_ft_and_fn[0]/current_ft_and_fn[1];
        current_error = fabs(current_angle - previous_angle);
        counter++;
    }

    assert(current_angle>angle1 && current_angle<angle2);
    assert(current_error < tolerance);

    return current_angle;
}


std::vector<double> DiscreteSystemForceCalculator::GetExtremalAngles(unsigned index, std::vector<double> samplingAngles)
{
    std::vector<double> extremal_angles;
    std::vector<double> ft_and_fn(2);
    std::vector<double> tangential_force(samplingAngles.size());

    for (unsigned i=0; i<samplingAngles.size(); i++)
    {
        ft_and_fn = CalculateFtAndFn(index, samplingAngles[i]);
        tangential_force[i] = ft_and_fn[0];
    }

    unsigned n = samplingAngles.size()-1;

    for (unsigned i=0; i<n; i++)
    {
        if ( ( tangential_force[i%n]>0 && tangential_force[(i+1)%n]<0 ) ||
             ( tangential_force[i%n]<0 && tangential_force[(i+1)%n]>0 ) )
        {
            double next_extremal_angle;

            // If we are in the interval that crosses the branch line at pi,
            // then subtract 2*pi from the positive angle
            if (i==n-1)
            {
                samplingAngles[i%n] -= 2*M_PI;
            }

            if (samplingAngles[(i+1)%n] - samplingAngles[i%n] < 2*mEpsilon + 1e-6 )
            {
                // If we find a jump through zero, then the local extremum is
                // simply at the mid-point of the interval
                next_extremal_angle = 0.5*(samplingAngles[(i+1)%n] + samplingAngles[i%n]);
            }
            else
            {
                // Otherwise we need to find it using Newton's method
                next_extremal_angle = GetLocalExtremum(index, samplingAngles[i%n], samplingAngles[(i+1)%n]);
            }

            if (next_extremal_angle <= -M_PI)
            {
                next_extremal_angle += 2*M_PI;
            }
            assert(next_extremal_angle>-M_PI && next_extremal_angle<=M_PI);
            extremal_angles.push_back(next_extremal_angle);
        }
    }

    return extremal_angles;
}
