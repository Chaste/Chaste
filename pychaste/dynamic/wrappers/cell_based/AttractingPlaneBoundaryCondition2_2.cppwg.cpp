#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PythonUblasObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AttractingPlaneBoundaryCondition.hpp"

#include "AttractingPlaneBoundaryCondition2_2.cppwg.hpp"

namespace py = pybind11;
typedef AttractingPlaneBoundaryCondition<2,2 > AttractingPlaneBoundaryCondition2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AttractingPlaneBoundaryCondition2_2_Overrides : public AttractingPlaneBoundaryCondition2_2{
    public:
    using AttractingPlaneBoundaryCondition2_2::AttractingPlaneBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            AttractingPlaneBoundaryCondition2_2,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            AttractingPlaneBoundaryCondition2_2,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AttractingPlaneBoundaryCondition2_2,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_AttractingPlaneBoundaryCondition2_2_class(py::module &m){
py::class_<AttractingPlaneBoundaryCondition2_2 , AttractingPlaneBoundaryCondition2_2_Overrides , boost::shared_ptr<AttractingPlaneBoundaryCondition2_2 >  , AbstractCellPopulationBoundaryCondition<2>  >(m, "AttractingPlaneBoundaryCondition2_2")
        .def(py::init<::AbstractCellPopulation<2> *, ::boost::numeric::ublas::c_vector<double, 2>, ::boost::numeric::ublas::c_vector<double, 2> >(), py::arg("pCellPopulation"), py::arg("point"), py::arg("normal"))
        .def(
            "rGetPointOnPlane",
            (::boost::numeric::ublas::c_vector<double, 2> const &(AttractingPlaneBoundaryCondition2_2::*)() const ) &AttractingPlaneBoundaryCondition2_2::rGetPointOnPlane,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetNormalToPlane",
            (::boost::numeric::ublas::c_vector<double, 2> const &(AttractingPlaneBoundaryCondition2_2::*)() const ) &AttractingPlaneBoundaryCondition2_2::rGetNormalToPlane,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetUseJiggledNodesOnPlane",
            (void(AttractingPlaneBoundaryCondition2_2::*)(bool)) &AttractingPlaneBoundaryCondition2_2::SetUseJiggledNodesOnPlane,
            " " , py::arg("useJiggledNodesOnPlane") )
        .def(
            "GetUseJiggledNodesOnPlane",
            (bool(AttractingPlaneBoundaryCondition2_2::*)()) &AttractingPlaneBoundaryCondition2_2::GetUseJiggledNodesOnPlane,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(AttractingPlaneBoundaryCondition2_2::*)(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const &)) &AttractingPlaneBoundaryCondition2_2::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(AttractingPlaneBoundaryCondition2_2::*)()) &AttractingPlaneBoundaryCondition2_2::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(AttractingPlaneBoundaryCondition2_2::*)(::out_stream &)) &AttractingPlaneBoundaryCondition2_2::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetPointOnPlane",
            (void(AttractingPlaneBoundaryCondition2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &)) &AttractingPlaneBoundaryCondition2_2::SetPointOnPlane,
            " " , py::arg("rPoint") )
        .def(
            "SetAttractionThreshold",
            (void(AttractingPlaneBoundaryCondition2_2::*)(double)) &AttractingPlaneBoundaryCondition2_2::SetAttractionThreshold,
            " " , py::arg("attractionThreshold") )
    ;
}
