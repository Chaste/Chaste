#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"

#include "ShortAxisVertexBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef ShortAxisVertexBasedDivisionRule<2 > ShortAxisVertexBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class ShortAxisVertexBasedDivisionRule2_Overrides : public ShortAxisVertexBasedDivisionRule2{
    public:
    using ShortAxisVertexBasedDivisionRule2::ShortAxisVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ShortAxisVertexBasedDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_ShortAxisVertexBasedDivisionRule2_class(py::module &m){
py::class_<ShortAxisVertexBasedDivisionRule2 , ShortAxisVertexBasedDivisionRule2_Overrides , boost::shared_ptr<ShortAxisVertexBasedDivisionRule2 >  , AbstractVertexBasedDivisionRule<2>  >(m, "ShortAxisVertexBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(ShortAxisVertexBasedDivisionRule2::*)(::CellPtr, ::VertexBasedCellPopulation<2> &)) &ShortAxisVertexBasedDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
