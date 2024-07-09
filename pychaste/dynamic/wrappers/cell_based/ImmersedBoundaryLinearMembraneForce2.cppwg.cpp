#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"

#include "ImmersedBoundaryLinearMembraneForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearMembraneForce<2 > ImmersedBoundaryLinearMembraneForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearMembraneForce2_Overrides : public ImmersedBoundaryLinearMembraneForce2{
    public:
    using ImmersedBoundaryLinearMembraneForce2::ImmersedBoundaryLinearMembraneForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearMembraneForce2_class(py::module &m){
py::class_<ImmersedBoundaryLinearMembraneForce2 , ImmersedBoundaryLinearMembraneForce2_Overrides , boost::shared_ptr<ImmersedBoundaryLinearMembraneForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryLinearMembraneForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryLinearMembraneForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(::out_stream &)) &ImmersedBoundaryLinearMembraneForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetElementSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce2::*)() const ) &ImmersedBoundaryLinearMembraneForce2::GetElementSpringConst,
            " "  )
        .def(
            "SetElementSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(double)) &ImmersedBoundaryLinearMembraneForce2::SetElementSpringConst,
            " " , py::arg("elementSpringConst") )
        .def(
            "GetElementRestLength",
            (double(ImmersedBoundaryLinearMembraneForce2::*)() const ) &ImmersedBoundaryLinearMembraneForce2::GetElementRestLength,
            " "  )
        .def(
            "SetElementRestLength",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(double)) &ImmersedBoundaryLinearMembraneForce2::SetElementRestLength,
            " " , py::arg("elementRestLength") )
        .def(
            "GetLaminaSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce2::*)() const ) &ImmersedBoundaryLinearMembraneForce2::GetLaminaSpringConst,
            " "  )
        .def(
            "SetLaminaSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(double)) &ImmersedBoundaryLinearMembraneForce2::SetLaminaSpringConst,
            " " , py::arg("laminaSpringConst") )
        .def(
            "GetLaminaRestLength",
            (double(ImmersedBoundaryLinearMembraneForce2::*)() const ) &ImmersedBoundaryLinearMembraneForce2::GetLaminaRestLength,
            " "  )
        .def(
            "SetLaminaRestLength",
            (void(ImmersedBoundaryLinearMembraneForce2::*)(double)) &ImmersedBoundaryLinearMembraneForce2::SetLaminaRestLength,
            " " , py::arg("laminaRestLength") )
    ;
}
