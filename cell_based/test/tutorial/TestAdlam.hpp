#ifndef TESTADLAM_HPP_
#define TESTADLAM_HPP_


#include "AbstractCellBasedTestSuite.hpp"
#include "CellDataItemWriter.hpp"

#include "BloatModifier.hpp"
#include "CellLabelWriter.hpp"
#include "CellReaper.hpp"
#include "SusceptibleCellCycleModel.hpp"

#include "IntercellularPressureModifier.hpp"
#include "CellularOxygenationModifier.hpp"
#include "CellularInfectionModifier.hpp"
#include "CustomStoppingSimulation.hpp"


// Structurally it's more 'chaste-appropriate' to store thigns like latency
// or infected as cell properties.


constexpr unsigned int sample_step = 25;
constexpr double simulation_time = 144.0;
constexpr double initial_age_range = 7.0;
constexpr double radius = 6;
constexpr unsigned length = 20;

constexpr unsigned mercy = 20;

constexpr bool attenuate_inwards = false;

constexpr double preponderance = 0.35;
constexpr double mutation_risk = 0.15;
constexpr double latency_radial_factor = 0.0;

constexpr double bloat = 0.3;
const double contagion = 0.8 * std::exp(-10.0);
const double apotheosis = 0.15 * std::exp(-10.0);
const double quiescence = 1.5 * std::exp(-10.0);
const double consumption_rate = 0.5 * std::exp(-10.0);

constexpr double oxygen_attenuation = 0.4;
constexpr double oxygen_baseline = 0.6;
constexpr double oxygen_range = 1.0;
const double membrane_permeability = 7 * std::exp(-10.0);

inline std::string output_directory = "TestLargerInfectionWithCellDeath";


class TestAdlam : public AbstractCellBasedTestSuite
{
public:
    static void TestInfectionWithMechanics()
    {
        RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

        HoneycombMeshGenerator generator(length, length);
        const boost::shared_ptr<MutableMesh<2,2>> p_mesh = generator.GetCircularMesh(radius + 3);

        std::vector<CellPtr> cells;
        std::vector<unsigned> reals;

        cells.reserve(p_mesh->GetNumNodes());
        reals.reserve(p_mesh->GetNumNodes());

        for (MutableMesh<2,2>::NodeIterator p_node = p_mesh->GetNodeIteratorBegin();
            p_node != p_mesh->GetNodeIteratorEnd();
            ++p_node)
        {
            if (p_node->IsDeleted())
            {
                continue;
            }

            const c_vector<double, 2> location = p_node->rGetLocation();

            if (const double r_squared = location[0] * location[0] + location[1] * location[1];
                r_squared <= radius * radius + DBL_EPSILON)
            {
                reals.push_back(p_node->GetIndex());

                const auto p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
                auto* p_model = new SusceptibleCellCycleModel(mutation_risk);

                CellPtr p_cell(new Cell(p_state, p_model));

                p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
                p_cell->SetBirthTime(-(p_model->GetMDuration() + initial_age_range * p_rng->ranf()));

                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes populace(*p_mesh, cells, reals, false, 15.0, 15.0, bloat * 2);

        if constexpr (attenuate_inwards)
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                const boost::shared_ptr<CellData> p_data = (*pp_cell)->GetCellData();

                if (p_rng->ranf() < preponderance * std::exp(- r * latency_radial_factor))
                {
                    p_data->SetItem("latent", 1.0);
                }

                p_data->SetItem("oxygen", oxygen_baseline + oxygen_range * (1.0 - std::exp(- oxygen_attenuation * r)));
            }
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                const boost::shared_ptr<CellData> p_data = (*pp_cell)->GetCellData();

                if (p_rng->ranf() < preponderance * std::exp(- r * latency_radial_factor))
                {
                    p_data->SetItem("latent", 1.0);
                }

                p_data->SetItem("oxygen", oxygen_baseline + oxygen_range * std::exp(- oxygen_attenuation * r));
            }
        }

        populace.AddPopulationWriter<VoronoiDataWriter>();
        populace.SetBoundVoronoiTessellation(true);

        const auto p_blitter = new CellDataItemWriter<2, 2>("oxygen");
        populace.AddCellWriter(static_cast<boost::shared_ptr<CellDataItemWriter<2, 2>>>(p_blitter));

        const auto p_marker = new CellLabelWriter<2, 2>;
        populace.AddCellWriter(static_cast<boost::shared_ptr<CellLabelWriter<2, 2>>>(p_marker));

        CustomStoppingSimulation simulation(populace);

        simulation.SetMercyThreshold(mercy);
        simulation.SetUpdateCellPopulationRule(true);

        const MAKE_PTR(IntercellularPressureModifier, p_pressure_modifier);
        simulation.AddTopologyUpdateSimulationModifier(p_pressure_modifier);

        const MAKE_PTR_ARGS(BloatModifier, p_area_modifier, (bloat));
        simulation.AddSimulationModifier(p_area_modifier);

        const MAKE_PTR_ARGS(CellularOxygenationModifier, p_oxygenation_modifier, (membrane_permeability));
        simulation.AddSimulationModifier(p_oxygenation_modifier);

        const MAKE_PTR_ARGS(CellularInfectionModifier, p_infection_modifier,
            (contagion, apotheosis, quiescence, consumption_rate));
        simulation.AddSimulationModifier(p_infection_modifier);

        const MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulation.AddForce(p_force);

        const MAKE_PTR_ARGS(CellReaper, p_reaper,
            (static_cast<AbstractCellPopulation<2>*>(&populace), mercy));
        simulation.AddCellKiller(p_reaper);

        simulation.SetSamplingTimestepMultiple(sample_step);
        simulation.SetEndTime(simulation_time);
        simulation.SetOutputDirectory(output_directory);

        simulation.Solve();
    }

    static void nTestVoronoiMechanics()
    {
        HoneycombMeshGenerator generator(length, length);
        const boost::shared_ptr<MutableMesh<2,2>> p_mesh = generator.GetCircularMesh(radius);

        const unsigned n_cells = p_mesh->GetNumNodes();

        std::vector<CellPtr> cells;
        cells.reserve(n_cells);

        const auto p_generator = new CellsGenerator<UniformG1GenerationalCellCycleModel, 2>();
        p_generator->GenerateBasicRandom(cells, n_cells, CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        MeshBasedCellPopulation populace(*p_mesh, cells);
        populace.AddPopulationWriter<VoronoiDataWriter>();
        populace.SetBoundVoronoiTessellation(true);

        OffLatticeSimulation simulation(populace);
        simulation.SetUpdateCellPopulationRule(true);

        const MAKE_PTR(IntercellularPressureModifier, p_pressure_modifier);
        simulation.AddTopologyUpdateSimulationModifier(p_pressure_modifier);

        const MAKE_PTR(SimpleTargetAreaModifier<2>, p_area_modifier);
        simulation.AddSimulationModifier(p_area_modifier);

        const MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulation.AddForce(p_force);

        simulation.SetSamplingTimestepMultiple(sample_step);
        simulation.SetEndTime(simulation_time);
        simulation.SetOutputDirectory(output_directory);

        simulation.Solve();
    }

    static void nTestInfectionModifier()
    {
        RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

        HoneycombMeshGenerator generator(length, length);
        const boost::shared_ptr<MutableMesh<2,2>> p_mesh = generator.GetCircularMesh(radius);

        const unsigned n_cells = p_mesh->GetNumNodes();

        std::vector<CellPtr> cells;
        cells.reserve(n_cells);

        for (unsigned int i = 0; i < n_cells; i++)
        {
            const auto p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            auto* p_model = new SusceptibleCellCycleModel(mutation_risk);

            CellPtr p_cell(new Cell(p_state, p_model));

            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            p_cell->SetBirthTime(-(p_model->GetMDuration() + initial_age_range * p_rng->ranf()));

            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation populace(*p_mesh, cells);

        if constexpr (attenuate_inwards)
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                const boost::shared_ptr<CellData> p_data = (*pp_cell)->GetCellData();

                if (p_rng->ranf() < preponderance * std::exp(- r * latency_radial_factor))
                {
                    p_data->SetItem("latent", 1.0);
                }

                p_data->SetItem("oxygen", oxygen_baseline + oxygen_range * (1.0 - std::exp(- oxygen_attenuation * r)));
            }
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                const boost::shared_ptr<CellData> p_data = (*pp_cell)->GetCellData();

                if (p_rng->ranf() < preponderance * std::exp(- r * latency_radial_factor))
                {
                    p_data->SetItem("latent", 1.0);
                }

                p_data->SetItem("oxygen", oxygen_baseline + oxygen_range * std::exp(- oxygen_attenuation * r));
            }
        }

        populace.AddPopulationWriter<VoronoiDataWriter>();
        populace.SetBoundVoronoiTessellation(true);

        OffLatticeSimulation simulation(populace);

        const MAKE_PTR_ARGS(CellularOxygenationModifier, p_oxygenation_modifier, (membrane_permeability));
        simulation.AddSimulationModifier(p_oxygenation_modifier);

        const MAKE_PTR_ARGS(CellularInfectionModifier, p_infection_modifier,
            (contagion, apotheosis, quiescence, consumption_rate));
        simulation.AddSimulationModifier(p_infection_modifier);

        simulation.SetSamplingTimestepMultiple(sample_step);
        simulation.SetEndTime(simulation_time);
        simulation.SetOutputDirectory(output_directory);

        simulation.Solve();
    }

    static void nTestOxygenationModifier()
    {
        RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

        HoneycombMeshGenerator generator(length, length);
        const boost::shared_ptr<MutableMesh<2,2>> p_mesh = generator.GetCircularMesh(radius);

        const unsigned n_cells = p_mesh->GetNumNodes();

        std::vector<CellPtr> cells;
        cells.reserve(n_cells);

        for (unsigned int i = 0; i < n_cells; i++)
        {
            const auto p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            auto* p_model = new SusceptibleCellCycleModel();

            CellPtr p_cell(new Cell(p_state, p_model));

            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            p_cell->SetBirthTime(-(p_model->GetMDuration() + initial_age_range * p_rng->ranf()));

            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation populace(*p_mesh, cells);

        if constexpr (attenuate_inwards)
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                (*pp_cell)->GetCellData()->SetItem("oxygen", oxygen_baseline + oxygen_range * (1.0 - std::exp(- oxygen_attenuation * r)));
            }
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator pp_cell = populace.Begin();
                pp_cell != populace.End();
                ++pp_cell)
            {
                const double x = populace.GetLocationOfCellCentre(*pp_cell)[0];
                const double y = populace.GetLocationOfCellCentre(*pp_cell)[1];
                const double r = std::sqrt(x*x + y*y);

                (*pp_cell)->GetCellData()->SetItem("oxygen", oxygen_baseline + oxygen_range * std::exp(- oxygen_attenuation * r));
            }
        }

        populace.AddPopulationWriter<VoronoiDataWriter>();
        populace.SetBoundVoronoiTessellation(true);

        OffLatticeSimulation simulation(populace);

        const MAKE_PTR_ARGS(CellularOxygenationModifier, p_oxygenation_modifier, (membrane_permeability));
        simulation.AddSimulationModifier(p_oxygenation_modifier);

        simulation.SetSamplingTimestepMultiple(sample_step);
        simulation.SetEndTime(simulation_time);
        simulation.SetOutputDirectory(output_directory);

        simulation.Solve();
    }
};


#endif // TESTADLAM_HPP_
