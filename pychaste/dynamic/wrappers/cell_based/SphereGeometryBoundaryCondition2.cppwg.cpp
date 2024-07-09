#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PythonUblasObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SphereGeometryBoundaryCondition.hpp"

#include "SphereGeometryBoundaryCondition2.cppwg.hpp"

namespace py = pybind11;
typedef SphereGeometryBoundaryCondition<2 > SphereGeometryBoundaryCondition2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SphereGeometryBoundaryCondition2_Overrides : public SphereGeometryBoundaryCondition2{
    public:
    using SphereGeometryBoundaryCondition2::SphereGeometryBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition2,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            SphereGeometryBoundaryCondition2,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition2,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_SphereGeometryBoundaryCondition2_class(py::module &m){
py::class_<SphereGeometryBoundaryCondition2 , SphereGeometryBoundaryCondition2_Overrides , boost::shared_ptr<SphereGeometryBoundaryCondition2 >  , AbstractCellPopulationBoundaryCondition<2>  >(m, "SphereGeometryBoundaryCondition2")
        .def(py::init<::AbstractCellPopulation<2> *, ::boost::numeric::ublas::c_vector<double, 2>, double, double >(), py::arg("pCellPopulation"), py::arg("centre"), py::arg("radius"), py::arg("distance") = 1.0000000000000001E-5)
        .def(
            "rGetCentreOfSphere",
            (::boost::numeric::ublas::c_vector<double, 2> const &(SphereGeometryBoundaryCondition2::*)() const ) &SphereGeometryBoundaryCondition2::rGetCentreOfSphere,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetRadiusOfSphere",
            (double(SphereGeometryBoundaryCondition2::*)() const ) &SphereGeometryBoundaryCondition2::GetRadiusOfSphere,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(SphereGeometryBoundaryCondition2::*)(::std::map<Node<2> *, boost::numeric::ublas::c_vector<double, 2>> const &)) &SphereGeometryBoundaryCondition2::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(SphereGeometryBoundaryCondition2::*)()) &SphereGeometryBoundaryCondition2::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(SphereGeometryBoundaryCondition2::*)(::out_stream &)) &SphereGeometryBoundaryCondition2::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
