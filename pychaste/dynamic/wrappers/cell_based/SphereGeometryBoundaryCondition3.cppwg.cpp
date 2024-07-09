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

#include "SphereGeometryBoundaryCondition3.cppwg.hpp"

namespace py = pybind11;
typedef SphereGeometryBoundaryCondition<3 > SphereGeometryBoundaryCondition3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SphereGeometryBoundaryCondition3_Overrides : public SphereGeometryBoundaryCondition3{
    public:
    using SphereGeometryBoundaryCondition3::SphereGeometryBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition3,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            SphereGeometryBoundaryCondition3,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition3,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_SphereGeometryBoundaryCondition3_class(py::module &m){
py::class_<SphereGeometryBoundaryCondition3 , SphereGeometryBoundaryCondition3_Overrides , boost::shared_ptr<SphereGeometryBoundaryCondition3 >  , AbstractCellPopulationBoundaryCondition<3>  >(m, "SphereGeometryBoundaryCondition3")
        .def(py::init<::AbstractCellPopulation<3> *, ::boost::numeric::ublas::c_vector<double, 3>, double, double >(), py::arg("pCellPopulation"), py::arg("centre"), py::arg("radius"), py::arg("distance") = 1.0000000000000001E-5)
        .def(
            "rGetCentreOfSphere",
            (::boost::numeric::ublas::c_vector<double, 3> const &(SphereGeometryBoundaryCondition3::*)() const ) &SphereGeometryBoundaryCondition3::rGetCentreOfSphere,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetRadiusOfSphere",
            (double(SphereGeometryBoundaryCondition3::*)() const ) &SphereGeometryBoundaryCondition3::GetRadiusOfSphere,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(SphereGeometryBoundaryCondition3::*)(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const &)) &SphereGeometryBoundaryCondition3::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(SphereGeometryBoundaryCondition3::*)()) &SphereGeometryBoundaryCondition3::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(SphereGeometryBoundaryCondition3::*)(::out_stream &)) &SphereGeometryBoundaryCondition3::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
