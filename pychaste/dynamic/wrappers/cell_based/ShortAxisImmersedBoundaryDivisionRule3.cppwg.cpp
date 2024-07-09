#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule.hpp"

#include "ShortAxisImmersedBoundaryDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef ShortAxisImmersedBoundaryDivisionRule<3 > ShortAxisImmersedBoundaryDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class ShortAxisImmersedBoundaryDivisionRule3_Overrides : public ShortAxisImmersedBoundaryDivisionRule3{
    public:
    using ShortAxisImmersedBoundaryDivisionRule3::ShortAxisImmersedBoundaryDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            ShortAxisImmersedBoundaryDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_ShortAxisImmersedBoundaryDivisionRule3_class(py::module &m){
py::class_<ShortAxisImmersedBoundaryDivisionRule3 , ShortAxisImmersedBoundaryDivisionRule3_Overrides , boost::shared_ptr<ShortAxisImmersedBoundaryDivisionRule3 >  , AbstractImmersedBoundaryDivisionRule<3>  >(m, "ShortAxisImmersedBoundaryDivisionRule3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(ShortAxisImmersedBoundaryDivisionRule3::*)(::CellPtr, ::ImmersedBoundaryCellPopulation<3> &)) &ShortAxisImmersedBoundaryDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
