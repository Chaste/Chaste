#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"

#include "ImmersedBoundaryMorseInteractionForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseInteractionForce<3 > ImmersedBoundaryMorseInteractionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseInteractionForce3_Overrides : public ImmersedBoundaryMorseInteractionForce3{
    public:
    using ImmersedBoundaryMorseInteractionForce3::ImmersedBoundaryMorseInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryMorseInteractionForce3_class(py::module &m){
py::class_<ImmersedBoundaryMorseInteractionForce3 , ImmersedBoundaryMorseInteractionForce3_Overrides , boost::shared_ptr<ImmersedBoundaryMorseInteractionForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryMorseInteractionForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryMorseInteractionForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(::out_stream &)) &ImmersedBoundaryMorseInteractionForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetWellDepth",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetWellDepth,
            " "  )
        .def(
            "SetWellDepth",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetWellDepth,
            " " , py::arg("wellDepth") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetRestLength,
            " " , py::arg("restLength") )
        .def(
            "GetLaminaWellDepthMult",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetLaminaWellDepthMult,
            " "  )
        .def(
            "SetLaminaWellDepthMult",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetLaminaWellDepthMult,
            " " , py::arg("laminaWellDepthMult") )
        .def(
            "GetLaminaRestLengthMult",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetLaminaRestLengthMult,
            " "  )
        .def(
            "SetLaminaRestLengthMult",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetLaminaRestLengthMult,
            " " , py::arg("laminaRestLengthMult") )
        .def(
            "GetWellWidth",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetWellWidth,
            " "  )
        .def(
            "SetWellWidth",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetWellWidth,
            " " , py::arg("wellWidth") )
    ;
}
