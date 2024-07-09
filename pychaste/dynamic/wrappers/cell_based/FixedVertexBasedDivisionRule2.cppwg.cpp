#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedVertexBasedDivisionRule.hpp"

#include "FixedVertexBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef FixedVertexBasedDivisionRule<2 > FixedVertexBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class FixedVertexBasedDivisionRule2_Overrides : public FixedVertexBasedDivisionRule2{
    public:
    using FixedVertexBasedDivisionRule2::FixedVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            FixedVertexBasedDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_FixedVertexBasedDivisionRule2_class(py::module &m){
py::class_<FixedVertexBasedDivisionRule2 , FixedVertexBasedDivisionRule2_Overrides , boost::shared_ptr<FixedVertexBasedDivisionRule2 >  , AbstractVertexBasedDivisionRule<2>  >(m, "FixedVertexBasedDivisionRule2")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2> & >(), py::arg("rDivisionVector"))
        .def(
            "rGetDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2> const &(FixedVertexBasedDivisionRule2::*)() const ) &FixedVertexBasedDivisionRule2::rGetDivisionVector,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(FixedVertexBasedDivisionRule2::*)(::CellPtr, ::VertexBasedCellPopulation<2> &)) &FixedVertexBasedDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
