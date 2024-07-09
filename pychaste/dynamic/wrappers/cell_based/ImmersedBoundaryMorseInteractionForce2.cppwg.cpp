#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"

#include "ImmersedBoundaryMorseInteractionForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseInteractionForce<2 > ImmersedBoundaryMorseInteractionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseInteractionForce2_Overrides : public ImmersedBoundaryMorseInteractionForce2{
    public:
    using ImmersedBoundaryMorseInteractionForce2::ImmersedBoundaryMorseInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryMorseInteractionForce2_class(py::module &m){
py::class_<ImmersedBoundaryMorseInteractionForce2 , ImmersedBoundaryMorseInteractionForce2_Overrides , boost::shared_ptr<ImmersedBoundaryMorseInteractionForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryMorseInteractionForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryMorseInteractionForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(::out_stream &)) &ImmersedBoundaryMorseInteractionForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetWellDepth",
            (double(ImmersedBoundaryMorseInteractionForce2::*)() const ) &ImmersedBoundaryMorseInteractionForce2::GetWellDepth,
            " "  )
        .def(
            "SetWellDepth",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(double)) &ImmersedBoundaryMorseInteractionForce2::SetWellDepth,
            " " , py::arg("wellDepth") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryMorseInteractionForce2::*)() const ) &ImmersedBoundaryMorseInteractionForce2::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(double)) &ImmersedBoundaryMorseInteractionForce2::SetRestLength,
            " " , py::arg("restLength") )
        .def(
            "GetLaminaWellDepthMult",
            (double(ImmersedBoundaryMorseInteractionForce2::*)() const ) &ImmersedBoundaryMorseInteractionForce2::GetLaminaWellDepthMult,
            " "  )
        .def(
            "SetLaminaWellDepthMult",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(double)) &ImmersedBoundaryMorseInteractionForce2::SetLaminaWellDepthMult,
            " " , py::arg("laminaWellDepthMult") )
        .def(
            "GetLaminaRestLengthMult",
            (double(ImmersedBoundaryMorseInteractionForce2::*)() const ) &ImmersedBoundaryMorseInteractionForce2::GetLaminaRestLengthMult,
            " "  )
        .def(
            "SetLaminaRestLengthMult",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(double)) &ImmersedBoundaryMorseInteractionForce2::SetLaminaRestLengthMult,
            " " , py::arg("laminaRestLengthMult") )
        .def(
            "GetWellWidth",
            (double(ImmersedBoundaryMorseInteractionForce2::*)() const ) &ImmersedBoundaryMorseInteractionForce2::GetWellWidth,
            " "  )
        .def(
            "SetWellWidth",
            (void(ImmersedBoundaryMorseInteractionForce2::*)(double)) &ImmersedBoundaryMorseInteractionForce2::SetWellWidth,
            " " , py::arg("wellWidth") )
    ;
}
