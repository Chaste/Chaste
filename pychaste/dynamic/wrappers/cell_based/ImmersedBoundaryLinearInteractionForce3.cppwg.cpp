#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"

#include "ImmersedBoundaryLinearInteractionForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearInteractionForce<3 > ImmersedBoundaryLinearInteractionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearInteractionForce3_Overrides : public ImmersedBoundaryLinearInteractionForce3{
    public:
    using ImmersedBoundaryLinearInteractionForce3::ImmersedBoundaryLinearInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearInteractionForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearInteractionForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearInteractionForce3_class(py::module &m){
py::class_<ImmersedBoundaryLinearInteractionForce3 , ImmersedBoundaryLinearInteractionForce3_Overrides , boost::shared_ptr<ImmersedBoundaryLinearInteractionForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryLinearInteractionForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryLinearInteractionForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(::out_stream &)) &ImmersedBoundaryLinearInteractionForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetSpringConst",
            (double(ImmersedBoundaryLinearInteractionForce3::*)() const ) &ImmersedBoundaryLinearInteractionForce3::GetSpringConst,
            " "  )
        .def(
            "SetSpringConst",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(double)) &ImmersedBoundaryLinearInteractionForce3::SetSpringConst,
            " " , py::arg("springConst") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryLinearInteractionForce3::*)() const ) &ImmersedBoundaryLinearInteractionForce3::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(double)) &ImmersedBoundaryLinearInteractionForce3::SetRestLength,
            " " , py::arg("restLength") )
        .def(
            "GetLaminaSpringConstMult",
            (double(ImmersedBoundaryLinearInteractionForce3::*)() const ) &ImmersedBoundaryLinearInteractionForce3::GetLaminaSpringConstMult,
            " "  )
        .def(
            "SetLaminaSpringConstMult",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(double)) &ImmersedBoundaryLinearInteractionForce3::SetLaminaSpringConstMult,
            " " , py::arg("laminaSpringConstMult") )
        .def(
            "GetLaminaRestLengthMult",
            (double(ImmersedBoundaryLinearInteractionForce3::*)() const ) &ImmersedBoundaryLinearInteractionForce3::GetLaminaRestLengthMult,
            " "  )
        .def(
            "SetLaminaRestLengthMult",
            (void(ImmersedBoundaryLinearInteractionForce3::*)(double)) &ImmersedBoundaryLinearInteractionForce3::SetLaminaRestLengthMult,
            " " , py::arg("laminaRestLengthMult") )
    ;
}
