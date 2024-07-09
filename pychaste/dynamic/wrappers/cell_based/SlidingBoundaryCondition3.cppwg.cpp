#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PythonUblasObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SlidingBoundaryCondition.hpp"

#include "SlidingBoundaryCondition3.cppwg.hpp"

namespace py = pybind11;
typedef SlidingBoundaryCondition<3 > SlidingBoundaryCondition3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SlidingBoundaryCondition3_Overrides : public SlidingBoundaryCondition3{
    public:
    using SlidingBoundaryCondition3::SlidingBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            SlidingBoundaryCondition3,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            SlidingBoundaryCondition3,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SlidingBoundaryCondition3,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_SlidingBoundaryCondition3_class(py::module &m){
py::class_<SlidingBoundaryCondition3 , SlidingBoundaryCondition3_Overrides , boost::shared_ptr<SlidingBoundaryCondition3 >  , AbstractCellPopulationBoundaryCondition<3>  >(m, "SlidingBoundaryCondition3")
        .def(py::init<::AbstractCellPopulation<3> *, double >(), py::arg("pCellPopulation"), py::arg("threshold") = 0.80000000000000004)
        .def(
            "GetThreshold",
            (double(SlidingBoundaryCondition3::*)() const ) &SlidingBoundaryCondition3::GetThreshold,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(SlidingBoundaryCondition3::*)(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const &)) &SlidingBoundaryCondition3::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(SlidingBoundaryCondition3::*)()) &SlidingBoundaryCondition3::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(SlidingBoundaryCondition3::*)(::out_stream &)) &SlidingBoundaryCondition3::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
