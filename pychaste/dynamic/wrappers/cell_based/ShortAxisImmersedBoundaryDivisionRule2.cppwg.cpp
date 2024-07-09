#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule.hpp"

#include "ShortAxisImmersedBoundaryDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef ShortAxisImmersedBoundaryDivisionRule<2 > ShortAxisImmersedBoundaryDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class ShortAxisImmersedBoundaryDivisionRule2_Overrides : public ShortAxisImmersedBoundaryDivisionRule2{
    public:
    using ShortAxisImmersedBoundaryDivisionRule2::ShortAxisImmersedBoundaryDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ShortAxisImmersedBoundaryDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_ShortAxisImmersedBoundaryDivisionRule2_class(py::module &m){
py::class_<ShortAxisImmersedBoundaryDivisionRule2 , ShortAxisImmersedBoundaryDivisionRule2_Overrides , boost::shared_ptr<ShortAxisImmersedBoundaryDivisionRule2 >  , AbstractImmersedBoundaryDivisionRule<2>  >(m, "ShortAxisImmersedBoundaryDivisionRule2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(ShortAxisImmersedBoundaryDivisionRule2::*)(::CellPtr, ::ImmersedBoundaryCellPopulation<2> &)) &ShortAxisImmersedBoundaryDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
