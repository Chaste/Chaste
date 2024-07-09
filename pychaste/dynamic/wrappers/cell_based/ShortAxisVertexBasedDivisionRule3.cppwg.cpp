#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"

#include "ShortAxisVertexBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef ShortAxisVertexBasedDivisionRule<3 > ShortAxisVertexBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class ShortAxisVertexBasedDivisionRule3_Overrides : public ShortAxisVertexBasedDivisionRule3{
    public:
    using ShortAxisVertexBasedDivisionRule3::ShortAxisVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            ShortAxisVertexBasedDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_ShortAxisVertexBasedDivisionRule3_class(py::module &m){
py::class_<ShortAxisVertexBasedDivisionRule3 , ShortAxisVertexBasedDivisionRule3_Overrides , boost::shared_ptr<ShortAxisVertexBasedDivisionRule3 >  , AbstractVertexBasedDivisionRule<3>  >(m, "ShortAxisVertexBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(ShortAxisVertexBasedDivisionRule3::*)(::CellPtr, ::VertexBasedCellPopulation<3> &)) &ShortAxisVertexBasedDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
