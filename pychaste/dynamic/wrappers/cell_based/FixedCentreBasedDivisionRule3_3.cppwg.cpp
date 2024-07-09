#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedCentreBasedDivisionRule.hpp"

#include "FixedCentreBasedDivisionRule3_3.cppwg.hpp"

namespace py = pybind11;
typedef FixedCentreBasedDivisionRule<3,3 > FixedCentreBasedDivisionRule3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_;

class FixedCentreBasedDivisionRule3_3_Overrides : public FixedCentreBasedDivisionRule3_3{
    public:
    using FixedCentreBasedDivisionRule3_3::FixedCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<3, 3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_,
            FixedCentreBasedDivisionRule3_3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_FixedCentreBasedDivisionRule3_3_class(py::module &m){
py::class_<FixedCentreBasedDivisionRule3_3 , FixedCentreBasedDivisionRule3_3_Overrides , boost::shared_ptr<FixedCentreBasedDivisionRule3_3 >  , AbstractCentreBasedDivisionRule<3>  >(m, "FixedCentreBasedDivisionRule3_3")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 3> & >(), py::arg("rDaughterLocation"))
        .def(
            "rGetDaughterLocation",
            (::boost::numeric::ublas::c_vector<double, 3> const &(FixedCentreBasedDivisionRule3_3::*)() const ) &FixedCentreBasedDivisionRule3_3::rGetDaughterLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>>(FixedCentreBasedDivisionRule3_3::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<3, 3> &)) &FixedCentreBasedDivisionRule3_3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
