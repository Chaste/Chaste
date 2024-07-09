#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VonMisesVertexBasedDivisionRule.hpp"

#include "VonMisesVertexBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef VonMisesVertexBasedDivisionRule<2 > VonMisesVertexBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class VonMisesVertexBasedDivisionRule2_Overrides : public VonMisesVertexBasedDivisionRule2{
    public:
    using VonMisesVertexBasedDivisionRule2::VonMisesVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            VonMisesVertexBasedDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_VonMisesVertexBasedDivisionRule2_class(py::module &m){
py::class_<VonMisesVertexBasedDivisionRule2 , VonMisesVertexBasedDivisionRule2_Overrides , boost::shared_ptr<VonMisesVertexBasedDivisionRule2 >  , AbstractVertexBasedDivisionRule<2>  >(m, "VonMisesVertexBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "GetMeanParameter",
            (double(VonMisesVertexBasedDivisionRule2::*)()) &VonMisesVertexBasedDivisionRule2::GetMeanParameter,
            " "  )
        .def(
            "GetConcentrationParameter",
            (double(VonMisesVertexBasedDivisionRule2::*)()) &VonMisesVertexBasedDivisionRule2::GetConcentrationParameter,
            " "  )
        .def(
            "SetMeanParameter",
            (void(VonMisesVertexBasedDivisionRule2::*)(double)) &VonMisesVertexBasedDivisionRule2::SetMeanParameter,
            " " , py::arg("meanParameter") )
        .def(
            "SetConcentrationParameter",
            (void(VonMisesVertexBasedDivisionRule2::*)(double)) &VonMisesVertexBasedDivisionRule2::SetConcentrationParameter,
            " " , py::arg("concentrationParameter") )
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(VonMisesVertexBasedDivisionRule2::*)(::CellPtr, ::VertexBasedCellPopulation<2> &)) &VonMisesVertexBasedDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
