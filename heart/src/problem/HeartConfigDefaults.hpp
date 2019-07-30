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

#ifndef HEARTCONFIGDEFAULTS_HPP_
#define HEARTCONFIGDEFAULTS_HPP_

/**
 * @file
 * This file is designed to be included within HeartConfig.cpp, and defines a single function,
 * CreateDefaultParameters.  The only reason it is a separate file is for easy identification
 * of what the default parameters are.
 */


/**
 * @return the default Chaste parameters.
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
        <Purkinje>
            <SurfaceAreaToVolumeRatio unit="1/cm"> 2800 </SurfaceAreaToVolumeRatio>
            <Capacitance unit="uF/cm^2"> 1.0 </Capacitance>
            <Conductivity unit="mS/cm"> 1.75 </Conductivity>
        </Purkinje>
    </Physiological>

    <Numerical>
        <TimeSteps ode="0.01" pde="0.01" printing="0.01" unit="ms" />
        <KSPTolerances>
            <KSPAbsolute>2e-4</KSPAbsolute>
        </KSPTolerances>
        <KSPSolver>cg</KSPSolver>
        <KSPPreconditioner>bjacobi</KSPPreconditioner>
        <MeshPartitioning>parmetis</MeshPartitioning>
        <UseStateVariableInterpolation>no</UseStateVariableInterpolation>
    </Numerical>
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

    cp::purkinje_physiological_type purkinje_phys_params;
    XSD_CREATE_WITH_FIXED_ATTR1(cp::inverse_length_type, purkinje_Am, 2800, "1/cm");
    purkinje_phys_params.SurfaceAreaToVolumeRatio().set(purkinje_Am);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::capacitance_type, purkinje_Cm, 1.0, "uF/cm^2");
    purkinje_phys_params.Capacitance().set(purkinje_Cm);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::conductivity_type, purkinje_conductivity, 1.75, "mS/cm");
    purkinje_phys_params.Conductivity().set(purkinje_conductivity);
    phys_params.Purkinje().set(purkinje_phys_params);

    // Numerical parameters
    cp::numerical_type numerical_params;
    XSD_CREATE_WITH_FIXED_ATTR3(cp::time_steps_type, timesteps, 0.01, 0.01, 0.01, "ms");
    cp::ksp_tolerances_type tolerances;
    tolerances.KSPAbsolute().set(2e-4);
    cp::ksp_solver_type ksp_solver("cg");
    cp::ksp_preconditioner_type ksp_precond("bjacobi");
    cp::mesh_partitioning_type mesh_partitioning("parmetis");
    numerical_params.TimeSteps().set(timesteps);
    numerical_params.KSPTolerances().set(tolerances);
    numerical_params.KSPSolver().set(ksp_solver);
    numerical_params.KSPPreconditioner().set(ksp_precond);
    numerical_params.MeshPartitioning().set(mesh_partitioning);
    numerical_params.UseStateVariableInterpolation().set(cp::yesno_type::no);

    // Postprocessing - empty is equivalent to missing, so don't include it
    //cp::postprocessing_type postproc;

    // Full default parameters
    boost::shared_ptr<cp::chaste_parameters_type> p_defaults(new cp::chaste_parameters_type(phys_params, numerical_params));
    p_defaults->Simulation().set(simulation_params);
    //p_defaults->PostProcessing().set(postproc);
    return p_defaults;
}


/**
 * If the given parameter is not present in the user parameters, but is in the defaults,
 * then copy its value from defaults to user parameters.
 *
 * @param path  the XSD data model path to the given parameter
 */
#define MERGE_PARAM(path)                                 \
    if (!pParams->path().present()) {                     \
        if (pDefaults->path().present()) {                \
            pParams->path().set(pDefaults->path().get()); \
        }                                                 \
    }
/**
 * An "else if" clause that tests if the given parameter is present in the defaults.
 *
 * @param path  the XSD data model path to the given parameter
 */
#define ELSE_IF_DEFAULT(path)                             \
    else if (pDefaults->path().present())

/**
 * Merge the default parameters (as defined by CreateDefaultParameters above) into given
 * user parameters.  Any parameter that is in the defaults but not the user parameters
 * will have its value copied over.
 *
 * @param pParams  the user parameters
 * @param pDefaults  the default parameters, which must have been created by CreateDefaultParameters
 *    (or be a subset thereof) for this method to work as intended
 */
void MergeDefaults(boost::shared_ptr<cp::chaste_parameters_type> pParams,
                   boost::shared_ptr<cp::chaste_parameters_type> pDefaults)
{
    // Simulation and ResumeSimulation are mutually exclusive
    if (!pParams->ResumeSimulation().present())
    {
        MERGE_PARAM(Simulation)
        ELSE_IF_DEFAULT(Simulation) // Simulation() exists in both
        {
            MERGE_PARAM(Simulation()->SpaceDimension)
            MERGE_PARAM(Simulation()->Domain)
            MERGE_PARAM(Simulation()->IonicModels)
            ELSE_IF_DEFAULT(Simulation()->IonicModels) // Simulation()->IonicModels() exists in both
            {
                //MERGE_PARAM(Simulation()->IonicModels()->Default) // This must be present if IonicModels is
            }
            MERGE_PARAM(Simulation()->OutputDirectory)
            MERGE_PARAM(Simulation()->OutputFilenamePrefix)
        }

        // Physiological() is mandatory
        //MERGE_PARAM(Physiological)
        //ELSE_IF_DEFAULT(Physiological) // Physiological() exists in both
        {
            MERGE_PARAM(Physiological().IntracellularConductivities)
            MERGE_PARAM(Physiological().ExtracellularConductivities)
            MERGE_PARAM(Physiological().BathConductivity)
            MERGE_PARAM(Physiological().SurfaceAreaToVolumeRatio)
            MERGE_PARAM(Physiological().Capacitance)
            MERGE_PARAM(Physiological().Purkinje)
            ELSE_IF_DEFAULT(Physiological().Purkinje) // Physiological()->Purkinje() exists in both
            {
                MERGE_PARAM(Physiological().Purkinje()->SurfaceAreaToVolumeRatio)
                MERGE_PARAM(Physiological().Purkinje()->Capacitance)
                MERGE_PARAM(Physiological().Purkinje()->Conductivity)
            }
        }

        // Numerical() is mandatory
        //MERGE_PARAM(Numerical)
        //ELSE_IF_DEFAULT(Numerical) // Numerical() exists in both
        {
            MERGE_PARAM(Numerical().TimeSteps)
            MERGE_PARAM(Numerical().KSPTolerances)
            ELSE_IF_DEFAULT(Numerical().KSPTolerances) // Numerical()->KSPTolerances() exists in both
            {
                // Note that we aren't allowed both absolute and relative tolerances
                if (!pParams->Numerical().KSPTolerances()->KSPRelative().present())
                {
                    MERGE_PARAM(Numerical().KSPTolerances()->KSPAbsolute)
                }
            }
            MERGE_PARAM(Numerical().KSPSolver)
            MERGE_PARAM(Numerical().KSPPreconditioner)
            MERGE_PARAM(Numerical().MeshPartitioning)
            MERGE_PARAM(Numerical().UseStateVariableInterpolation)
        }
    }
}

#endif /*HEARTCONFIGDEFAULTS_HPP_*/
