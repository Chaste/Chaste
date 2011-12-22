/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef HEARTCONFIGDEFAULTS_HPP_
#define HEARTCONFIGDEFAULTS_HPP_

/**
 * @file
 * This file is designed to be included within HeartConfig.cpp, and defines a single function,
 * CreateDefaultParameters.  The only reason it is a separate file is for easy identification
 * of what the default parameters are.
 */


/**
 * Create the default Chaste parameters.
 *
 * It sets up an object equivalent to the following XML file:
 * \verbatim
<?xml version="1.0" encoding="UTF-8"?>
<ChasteParameters>
    <Simulation>
        <SpaceDimension>3</SpaceDimension>
        <Domain>Mono</Domain>
        <IonicModels>
            <Default><Hardcoded>LuoRudyI</Hardcoded></Default>
        </IonicModels>
        <OutputDirectory>ChasteResults</OutputDirectory>
        <OutputFilenamePrefix>SimulationResults</OutputFilenamePrefix>
    </Simulation>

    <Physiological>
        <IntracellularConductivities longi="1.75" trans="1.75" normal="1.75" unit="mS/cm" />
        <ExtracellularConductivities longi="7.0"  trans="7.0"  normal="7.0" unit="mS/cm" />
        <BathConductivity unit="mS/cm"> 7.0 </BathConductivity>
        <SurfaceAreaToVolumeRatio unit="1/cm"> 1400 </SurfaceAreaToVolumeRatio>
        <Capacitance unit="uF/cm^2"> 1.0 </Capacitance>
    </Physiological>

    <Numerical>
        <TimeSteps ode="0.01" pde="0.01" printing="0.01" unit="ms" />
        <KSPTolerances>
            <KSPAbsolute>2e-4</KSPAbsolute>
        </KSPTolerances>
        <KSPSolver>cg</KSPSolver>
        <KSPPreconditioner>bjacobi</KSPPreconditioner>
        <MeshPartitioning>metis</MeshPartitioning>
        <UseStateVariableInterpolation>no</UseStateVariableInterpolation>
    </Numerical>

    <PostProcessing>
    </PostProcessing>
</ChasteParameters>
\endverbatim
 */
boost::shared_ptr<cp::chaste_parameters_type> CreateDefaultParameters()
{
    // Simulation parameters
    cp::simulation_type simulation_params;
    simulation_params.SpaceDimension().set(3);
    cp::domain_type domain("Mono");
    simulation_params.Domain().set(domain);
    cp::ionic_model_selection_type default_ionic_model;
    cp::ionic_models_available_type ionic_model("LuoRudyI");
    default_ionic_model.Hardcoded().set(ionic_model);
    cp::ionic_models_type ionic_models(default_ionic_model);
    simulation_params.IonicModels().set(ionic_models);
    simulation_params.OutputDirectory().set("ChasteResults");
    simulation_params.OutputFilenamePrefix().set("SimulationResults");

    // Physiological parameters
    cp::physiological_type phys_params;
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra_conductivities,
                                1.75, 1.75, 1.75, "mS/cm");
    phys_params.IntracellularConductivities().set(intra_conductivities);
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra_conductivities,
                                7.0, 7.0, 7.0, "mS/cm");
    phys_params.ExtracellularConductivities().set(extra_conductivities);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::conductivity_type, bath_conductivity, 7.0, "mS/cm");
    phys_params.BathConductivity().set(bath_conductivity);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::inverse_length_type, surface_area_to_volume_ratio, 1400, "1/cm");
    phys_params.SurfaceAreaToVolumeRatio().set(surface_area_to_volume_ratio);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::capacitance_type, capacitance, 1.0, "uF/cm^2");
    phys_params.Capacitance().set(capacitance);

    // Numerical parameters
    cp::numerical_type numerical_params;
    XSD_CREATE_WITH_FIXED_ATTR3(cp::time_steps_type, timesteps, 0.01, 0.01, 0.01, "ms");
    cp::ksp_tolerances_type tolerances;
    tolerances.KSPAbsolute().set(2e-4);
    cp::ksp_solver_type ksp_solver("cg");
    cp::ksp_preconditioner_type ksp_precond("bjacobi");
    cp::mesh_partitioning_type mesh_partitioning("metis");
    numerical_params.TimeSteps().set(timesteps);
    numerical_params.KSPTolerances().set(tolerances);
    numerical_params.KSPSolver().set(ksp_solver);
    numerical_params.KSPPreconditioner().set(ksp_precond);
    numerical_params.MeshPartitioning().set(mesh_partitioning);
    numerical_params.UseStateVariableInterpolation().set(cp::yesno_type::no);

    // Postprocessing
    cp::postprocessing_type postproc;

    // Full default parameters
    boost::shared_ptr<cp::chaste_parameters_type> p_defaults(new cp::chaste_parameters_type(phys_params, numerical_params));
    p_defaults->Simulation().set(simulation_params);
    p_defaults->PostProcessing().set(postproc);
    return p_defaults;
}

#endif /*HEARTCONFIGDEFAULTS_HPP_*/
