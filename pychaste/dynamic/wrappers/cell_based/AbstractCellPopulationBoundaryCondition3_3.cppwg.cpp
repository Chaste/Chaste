#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "AbstractCellPopulationBoundaryCondition3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellPopulationBoundaryCondition<3,3 > AbstractCellPopulationBoundaryCondition3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellPopulationBoundaryCondition3_3_Overrides : public AbstractCellPopulationBoundaryCondition3_3{
    public:
    using AbstractCellPopulationBoundaryCondition3_3::AbstractCellPopulationBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const & rOldLocations) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationBoundaryCondition3_3,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellPopulationBoundaryCondition3_3,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationBoundaryCondition3_3,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_AbstractCellPopulationBoundaryCondition3_3_class(py::module &m){
py::class_<AbstractCellPopulationBoundaryCondition3_3 , AbstractCellPopulationBoundaryCondition3_3_Overrides , boost::shared_ptr<AbstractCellPopulationBoundaryCondition3_3 >   >(m, "AbstractCellPopulationBoundaryCondition3_3")
        .def(py::init<::AbstractCellPopulation<3> * >(), py::arg("pCellPopulation"))
        .def(
            "ImposeBoundaryCondition",
            (void(AbstractCellPopulationBoundaryCondition3_3::*)(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const &)) &AbstractCellPopulationBoundaryCondition3_3::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(AbstractCellPopulationBoundaryCondition3_3::*)()) &AbstractCellPopulationBoundaryCondition3_3::VerifyBoundaryCondition,
            " "  )
        .def(
            "GetCellPopulation",
            (::AbstractCellPopulation<3> const *(AbstractCellPopulationBoundaryCondition3_3::*)() const ) &AbstractCellPopulationBoundaryCondition3_3::GetCellPopulation,
            " "  , py::return_value_policy::reference)
        .def(
            "OutputCellPopulationBoundaryConditionInfo",
            (void(AbstractCellPopulationBoundaryCondition3_3::*)(::out_stream &)) &AbstractCellPopulationBoundaryCondition3_3::OutputCellPopulationBoundaryConditionInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(AbstractCellPopulationBoundaryCondition3_3::*)(::out_stream &)) &AbstractCellPopulationBoundaryCondition3_3::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
