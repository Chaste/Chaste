#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PythonUblasObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PlaneBoundaryCondition.hpp"

#include "PlaneBoundaryCondition3_3.cppwg.hpp"

namespace py = pybind11;
typedef PlaneBoundaryCondition<3,3 > PlaneBoundaryCondition3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class PlaneBoundaryCondition3_3_Overrides : public PlaneBoundaryCondition3_3{
    public:
    using PlaneBoundaryCondition3_3::PlaneBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            PlaneBoundaryCondition3_3,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            PlaneBoundaryCondition3_3,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PlaneBoundaryCondition3_3,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_PlaneBoundaryCondition3_3_class(py::module &m){
py::class_<PlaneBoundaryCondition3_3 , PlaneBoundaryCondition3_3_Overrides , boost::shared_ptr<PlaneBoundaryCondition3_3 >  , AbstractCellPopulationBoundaryCondition<3>  >(m, "PlaneBoundaryCondition3_3")
        .def(py::init<::AbstractCellPopulation<3> *, ::boost::numeric::ublas::c_vector<double, 3>, ::boost::numeric::ublas::c_vector<double, 3> >(), py::arg("pCellPopulation"), py::arg("point"), py::arg("normal"))
        .def(
            "rGetPointOnPlane",
            (::boost::numeric::ublas::c_vector<double, 3> const &(PlaneBoundaryCondition3_3::*)() const ) &PlaneBoundaryCondition3_3::rGetPointOnPlane,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetNormalToPlane",
            (::boost::numeric::ublas::c_vector<double, 3> const &(PlaneBoundaryCondition3_3::*)() const ) &PlaneBoundaryCondition3_3::rGetNormalToPlane,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetUseJiggledNodesOnPlane",
            (void(PlaneBoundaryCondition3_3::*)(bool)) &PlaneBoundaryCondition3_3::SetUseJiggledNodesOnPlane,
            " " , py::arg("useJiggledNodesOnPlane") )
        .def(
            "GetUseJiggledNodesOnPlane",
            (bool(PlaneBoundaryCondition3_3::*)()) &PlaneBoundaryCondition3_3::GetUseJiggledNodesOnPlane,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(PlaneBoundaryCondition3_3::*)(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const &)) &PlaneBoundaryCondition3_3::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(PlaneBoundaryCondition3_3::*)()) &PlaneBoundaryCondition3_3::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(PlaneBoundaryCondition3_3::*)(::out_stream &)) &PlaneBoundaryCondition3_3::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
