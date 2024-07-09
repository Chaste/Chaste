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

#include "SlidingBoundaryCondition2.cppwg.hpp"

namespace py = pybind11;
typedef SlidingBoundaryCondition<2 > SlidingBoundaryCondition2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SlidingBoundaryCondition2_Overrides : public SlidingBoundaryCondition2{
    public:
    using SlidingBoundaryCondition2::SlidingBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            SlidingBoundaryCondition2,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            SlidingBoundaryCondition2,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SlidingBoundaryCondition2,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_SlidingBoundaryCondition2_class(py::module &m){
py::class_<SlidingBoundaryCondition2 , SlidingBoundaryCondition2_Overrides , boost::shared_ptr<SlidingBoundaryCondition2 >  , AbstractCellPopulationBoundaryCondition<2>  >(m, "SlidingBoundaryCondition2")
        .def(py::init<::AbstractCellPopulation<2> *, double >(), py::arg("pCellPopulation"), py::arg("threshold") = 0.80000000000000004)
        .def(
            "GetThreshold",
            (double(SlidingBoundaryCondition2::*)() const ) &SlidingBoundaryCondition2::GetThreshold,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(SlidingBoundaryCondition2::*)(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const &)) &SlidingBoundaryCondition2::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(SlidingBoundaryCondition2::*)()) &SlidingBoundaryCondition2::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(SlidingBoundaryCondition2::*)(::out_stream &)) &SlidingBoundaryCondition2::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
