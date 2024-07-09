#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"

#include "ImmersedBoundaryLinearInteractionForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearInteractionForce<2 > ImmersedBoundaryLinearInteractionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearInteractionForce2_Overrides : public ImmersedBoundaryLinearInteractionForce2{
    public:
    using ImmersedBoundaryLinearInteractionForce2::ImmersedBoundaryLinearInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearInteractionForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearInteractionForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearInteractionForce2_class(py::module &m){
py::class_<ImmersedBoundaryLinearInteractionForce2 , ImmersedBoundaryLinearInteractionForce2_Overrides , boost::shared_ptr<ImmersedBoundaryLinearInteractionForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryLinearInteractionForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryLinearInteractionForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(::out_stream &)) &ImmersedBoundaryLinearInteractionForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetSpringConst",
            (double(ImmersedBoundaryLinearInteractionForce2::*)() const ) &ImmersedBoundaryLinearInteractionForce2::GetSpringConst,
            " "  )
        .def(
            "SetSpringConst",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(double)) &ImmersedBoundaryLinearInteractionForce2::SetSpringConst,
            " " , py::arg("springConst") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryLinearInteractionForce2::*)() const ) &ImmersedBoundaryLinearInteractionForce2::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(double)) &ImmersedBoundaryLinearInteractionForce2::SetRestLength,
            " " , py::arg("restLength") )
        .def(
            "GetLaminaSpringConstMult",
            (double(ImmersedBoundaryLinearInteractionForce2::*)() const ) &ImmersedBoundaryLinearInteractionForce2::GetLaminaSpringConstMult,
            " "  )
        .def(
            "SetLaminaSpringConstMult",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(double)) &ImmersedBoundaryLinearInteractionForce2::SetLaminaSpringConstMult,
            " " , py::arg("laminaSpringConstMult") )
        .def(
            "GetLaminaRestLengthMult",
            (double(ImmersedBoundaryLinearInteractionForce2::*)() const ) &ImmersedBoundaryLinearInteractionForce2::GetLaminaRestLengthMult,
            " "  )
        .def(
            "SetLaminaRestLengthMult",
            (void(ImmersedBoundaryLinearInteractionForce2::*)(double)) &ImmersedBoundaryLinearInteractionForce2::SetLaminaRestLengthMult,
            " " , py::arg("laminaRestLengthMult") )
    ;
}
