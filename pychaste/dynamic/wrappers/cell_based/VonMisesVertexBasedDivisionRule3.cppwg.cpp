#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VonMisesVertexBasedDivisionRule.hpp"

#include "VonMisesVertexBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef VonMisesVertexBasedDivisionRule<3 > VonMisesVertexBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class VonMisesVertexBasedDivisionRule3_Overrides : public VonMisesVertexBasedDivisionRule3{
    public:
    using VonMisesVertexBasedDivisionRule3::VonMisesVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            VonMisesVertexBasedDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_VonMisesVertexBasedDivisionRule3_class(py::module &m){
py::class_<VonMisesVertexBasedDivisionRule3 , VonMisesVertexBasedDivisionRule3_Overrides , boost::shared_ptr<VonMisesVertexBasedDivisionRule3 >  , AbstractVertexBasedDivisionRule<3>  >(m, "VonMisesVertexBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "GetMeanParameter",
            (double(VonMisesVertexBasedDivisionRule3::*)()) &VonMisesVertexBasedDivisionRule3::GetMeanParameter,
            " "  )
        .def(
            "GetConcentrationParameter",
            (double(VonMisesVertexBasedDivisionRule3::*)()) &VonMisesVertexBasedDivisionRule3::GetConcentrationParameter,
            " "  )
        .def(
            "SetMeanParameter",
            (void(VonMisesVertexBasedDivisionRule3::*)(double)) &VonMisesVertexBasedDivisionRule3::SetMeanParameter,
            " " , py::arg("meanParameter") )
        .def(
            "SetConcentrationParameter",
            (void(VonMisesVertexBasedDivisionRule3::*)(double)) &VonMisesVertexBasedDivisionRule3::SetConcentrationParameter,
            " " , py::arg("concentrationParameter") )
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(VonMisesVertexBasedDivisionRule3::*)(::CellPtr, ::VertexBasedCellPopulation<3> &)) &VonMisesVertexBasedDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
