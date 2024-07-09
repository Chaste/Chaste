#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryLinearDifferentialAdhesionForce.hpp"

#include "ImmersedBoundaryLinearDifferentialAdhesionForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearDifferentialAdhesionForce<2 > ImmersedBoundaryLinearDifferentialAdhesionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearDifferentialAdhesionForce2_Overrides : public ImmersedBoundaryLinearDifferentialAdhesionForce2{
    public:
    using ImmersedBoundaryLinearDifferentialAdhesionForce2::ImmersedBoundaryLinearDifferentialAdhesionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearDifferentialAdhesionForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearDifferentialAdhesionForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearDifferentialAdhesionForce2_class(py::module &m){
py::class_<ImmersedBoundaryLinearDifferentialAdhesionForce2 , ImmersedBoundaryLinearDifferentialAdhesionForce2_Overrides , boost::shared_ptr<ImmersedBoundaryLinearDifferentialAdhesionForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryLinearDifferentialAdhesionForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(::out_stream &)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetLabelledCellToLabelledCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetLabelledCellToLabelledCellSpringConst,
            " "  )
        .def(
            "SetLabelledCellToLabelledCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetLabelledCellToLabelledCellSpringConst,
            " " , py::arg("labelledCellToLabelledCellSpringConst") )
        .def(
            "GetLabelledCellToCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetLabelledCellToCellSpringConst,
            " "  )
        .def(
            "SetLabelledCellToCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetLabelledCellToCellSpringConst,
            " " , py::arg("labelledCellToCellSpringConst") )
        .def(
            "GetCellToCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetCellToCellSpringConst,
            " "  )
        .def(
            "SetCellToCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetCellToCellSpringConst,
            " " , py::arg("cellToCellSpringConst") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetRestLength,
            " " , py::arg("restLength") )
    ;
}
