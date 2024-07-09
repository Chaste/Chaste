#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"

#include "ImmersedBoundaryMorseMembraneForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseMembraneForce<3 > ImmersedBoundaryMorseMembraneForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseMembraneForce3_Overrides : public ImmersedBoundaryMorseMembraneForce3{
    public:
    using ImmersedBoundaryMorseMembraneForce3::ImmersedBoundaryMorseMembraneForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseMembraneForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseMembraneForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryMorseMembraneForce3_class(py::module &m){
py::class_<ImmersedBoundaryMorseMembraneForce3 , ImmersedBoundaryMorseMembraneForce3_Overrides , boost::shared_ptr<ImmersedBoundaryMorseMembraneForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryMorseMembraneForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryMorseMembraneForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(::out_stream &)) &ImmersedBoundaryMorseMembraneForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetElementWellDepth",
            (double(ImmersedBoundaryMorseMembraneForce3::*)() const ) &ImmersedBoundaryMorseMembraneForce3::GetElementWellDepth,
            " "  )
        .def(
            "SetElementWellDepth",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(double)) &ImmersedBoundaryMorseMembraneForce3::SetElementWellDepth,
            " " , py::arg("elementWellDepth") )
        .def(
            "GetElementRestLength",
            (double(ImmersedBoundaryMorseMembraneForce3::*)() const ) &ImmersedBoundaryMorseMembraneForce3::GetElementRestLength,
            " "  )
        .def(
            "SetElementRestLength",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(double)) &ImmersedBoundaryMorseMembraneForce3::SetElementRestLength,
            " " , py::arg("elementRestLength") )
        .def(
            "GetLaminaWellDepth",
            (double(ImmersedBoundaryMorseMembraneForce3::*)() const ) &ImmersedBoundaryMorseMembraneForce3::GetLaminaWellDepth,
            " "  )
        .def(
            "SetLaminaWellDepth",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(double)) &ImmersedBoundaryMorseMembraneForce3::SetLaminaWellDepth,
            " " , py::arg("laminaWellDepth") )
        .def(
            "GetLaminaRestLength",
            (double(ImmersedBoundaryMorseMembraneForce3::*)() const ) &ImmersedBoundaryMorseMembraneForce3::GetLaminaRestLength,
            " "  )
        .def(
            "SetLaminaRestLength",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(double)) &ImmersedBoundaryMorseMembraneForce3::SetLaminaRestLength,
            " " , py::arg("laminaRestLength") )
        .def(
            "GetWellWidth",
            (double(ImmersedBoundaryMorseMembraneForce3::*)() const ) &ImmersedBoundaryMorseMembraneForce3::GetWellWidth,
            " "  )
        .def(
            "SetWellWidth",
            (void(ImmersedBoundaryMorseMembraneForce3::*)(double)) &ImmersedBoundaryMorseMembraneForce3::SetWellWidth,
            " " , py::arg("wellWidth") )
    ;
}
