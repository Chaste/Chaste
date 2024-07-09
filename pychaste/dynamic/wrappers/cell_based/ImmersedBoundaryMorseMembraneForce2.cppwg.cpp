#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"

#include "ImmersedBoundaryMorseMembraneForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseMembraneForce<2 > ImmersedBoundaryMorseMembraneForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseMembraneForce2_Overrides : public ImmersedBoundaryMorseMembraneForce2{
    public:
    using ImmersedBoundaryMorseMembraneForce2::ImmersedBoundaryMorseMembraneForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseMembraneForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseMembraneForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryMorseMembraneForce2_class(py::module &m){
py::class_<ImmersedBoundaryMorseMembraneForce2 , ImmersedBoundaryMorseMembraneForce2_Overrides , boost::shared_ptr<ImmersedBoundaryMorseMembraneForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryMorseMembraneForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryMorseMembraneForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(::out_stream &)) &ImmersedBoundaryMorseMembraneForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetElementWellDepth",
            (double(ImmersedBoundaryMorseMembraneForce2::*)() const ) &ImmersedBoundaryMorseMembraneForce2::GetElementWellDepth,
            " "  )
        .def(
            "SetElementWellDepth",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(double)) &ImmersedBoundaryMorseMembraneForce2::SetElementWellDepth,
            " " , py::arg("elementWellDepth") )
        .def(
            "GetElementRestLength",
            (double(ImmersedBoundaryMorseMembraneForce2::*)() const ) &ImmersedBoundaryMorseMembraneForce2::GetElementRestLength,
            " "  )
        .def(
            "SetElementRestLength",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(double)) &ImmersedBoundaryMorseMembraneForce2::SetElementRestLength,
            " " , py::arg("elementRestLength") )
        .def(
            "GetLaminaWellDepth",
            (double(ImmersedBoundaryMorseMembraneForce2::*)() const ) &ImmersedBoundaryMorseMembraneForce2::GetLaminaWellDepth,
            " "  )
        .def(
            "SetLaminaWellDepth",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(double)) &ImmersedBoundaryMorseMembraneForce2::SetLaminaWellDepth,
            " " , py::arg("laminaWellDepth") )
        .def(
            "GetLaminaRestLength",
            (double(ImmersedBoundaryMorseMembraneForce2::*)() const ) &ImmersedBoundaryMorseMembraneForce2::GetLaminaRestLength,
            " "  )
        .def(
            "SetLaminaRestLength",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(double)) &ImmersedBoundaryMorseMembraneForce2::SetLaminaRestLength,
            " " , py::arg("laminaRestLength") )
        .def(
            "GetWellWidth",
            (double(ImmersedBoundaryMorseMembraneForce2::*)() const ) &ImmersedBoundaryMorseMembraneForce2::GetWellWidth,
            " "  )
        .def(
            "SetWellWidth",
            (void(ImmersedBoundaryMorseMembraneForce2::*)(double)) &ImmersedBoundaryMorseMembraneForce2::SetWellWidth,
            " " , py::arg("wellWidth") )
    ;
}
