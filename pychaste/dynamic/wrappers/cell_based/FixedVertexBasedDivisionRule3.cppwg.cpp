#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedVertexBasedDivisionRule.hpp"

#include "FixedVertexBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef FixedVertexBasedDivisionRule<3 > FixedVertexBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class FixedVertexBasedDivisionRule3_Overrides : public FixedVertexBasedDivisionRule3{
    public:
    using FixedVertexBasedDivisionRule3::FixedVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            FixedVertexBasedDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_FixedVertexBasedDivisionRule3_class(py::module &m){
py::class_<FixedVertexBasedDivisionRule3 , FixedVertexBasedDivisionRule3_Overrides , boost::shared_ptr<FixedVertexBasedDivisionRule3 >  , AbstractVertexBasedDivisionRule<3>  >(m, "FixedVertexBasedDivisionRule3")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 3> & >(), py::arg("rDivisionVector"))
        .def(
            "rGetDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3> const &(FixedVertexBasedDivisionRule3::*)() const ) &FixedVertexBasedDivisionRule3::rGetDivisionVector,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(FixedVertexBasedDivisionRule3::*)(::CellPtr, ::VertexBasedCellPopulation<3> &)) &FixedVertexBasedDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
