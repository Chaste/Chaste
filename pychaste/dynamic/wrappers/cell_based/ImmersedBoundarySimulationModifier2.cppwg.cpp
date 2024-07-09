#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"

#include "ImmersedBoundarySimulationModifier2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundarySimulationModifier<2 > ImmersedBoundarySimulationModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundarySimulationModifier2_Overrides : public ImmersedBoundarySimulationModifier2{
    public:
    using ImmersedBoundarySimulationModifier2::ImmersedBoundarySimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundarySimulationModifier2_class(py::module &m){
py::class_<ImmersedBoundarySimulationModifier2 , ImmersedBoundarySimulationModifier2_Overrides , boost::shared_ptr<ImmersedBoundarySimulationModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "ImmersedBoundarySimulationModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ImmersedBoundarySimulationModifier2::*)(::AbstractCellPopulation<2> &)) &ImmersedBoundarySimulationModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ImmersedBoundarySimulationModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ImmersedBoundarySimulationModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ImmersedBoundarySimulationModifier2::*)(::out_stream &)) &ImmersedBoundarySimulationModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetNodeNeighbourUpdateFrequency",
            (void(ImmersedBoundarySimulationModifier2::*)(unsigned int)) &ImmersedBoundarySimulationModifier2::SetNodeNeighbourUpdateFrequency,
            " " , py::arg("newFrequency") )
        .def(
            "GetNodeNeighbourUpdateFrequency",
            (unsigned int(ImmersedBoundarySimulationModifier2::*)()) &ImmersedBoundarySimulationModifier2::GetNodeNeighbourUpdateFrequency,
            " "  )
        .def(
            "AddImmersedBoundaryForce",
            (void(ImmersedBoundarySimulationModifier2::*)(::boost::shared_ptr<AbstractImmersedBoundaryForce<2>>)) &ImmersedBoundarySimulationModifier2::AddImmersedBoundaryForce,
            " " , py::arg("pForce") )
        .def(
            "AddNormalNoise",
            (void(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::AddNormalNoise,
            " "  )
        .def(
            "GetZeroFieldSums",
            (bool(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetZeroFieldSums,
            " "  )
        .def(
            "SetZeroFieldSums",
            (void(ImmersedBoundarySimulationModifier2::*)(bool)) &ImmersedBoundarySimulationModifier2::SetZeroFieldSums,
            " " , py::arg("zeroFieldSums") )
        .def(
            "SetReynoldsNumber",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetReynoldsNumber,
            " " , py::arg("reynoldsNumber") )
        .def(
            "GetReynoldsNumber",
            (double(ImmersedBoundarySimulationModifier2::*)()) &ImmersedBoundarySimulationModifier2::GetReynoldsNumber,
            " "  )
        .def(
            "GetAdditiveNormalNoise",
            (bool(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetAdditiveNormalNoise,
            " "  )
        .def(
            "SetAdditiveNormalNoise",
            (void(ImmersedBoundarySimulationModifier2::*)(bool)) &ImmersedBoundarySimulationModifier2::SetAdditiveNormalNoise,
            " " , py::arg("additiveNormalNoise") )
        .def(
            "GetNoiseStrength",
            (double(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseStrength,
            " "  )
        .def(
            "SetNoiseStrength",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetNoiseStrength,
            " " , py::arg("noiseStrength") )
        .def(
            "GetNoiseSkip",
            (unsigned int(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseSkip,
            " "  )
        .def(
            "SetNoiseSkip",
            (void(ImmersedBoundarySimulationModifier2::*)(unsigned int)) &ImmersedBoundarySimulationModifier2::SetNoiseSkip,
            " " , py::arg("noiseSkip") )
        .def(
            "GetNoiseLengthScale",
            (double(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseLengthScale,
            " "  )
        .def(
            "SetNoiseLengthScale",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetNoiseLengthScale,
            " " , py::arg("noiseLengthScale") )
    ;
}
