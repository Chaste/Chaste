#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"

#include "RandomDirectionCentreBasedDivisionRule2_2.cppwg.hpp"

namespace py = pybind11;
typedef RandomDirectionCentreBasedDivisionRule<2,2 > RandomDirectionCentreBasedDivisionRule2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_;

class RandomDirectionCentreBasedDivisionRule2_2_Overrides : public RandomDirectionCentreBasedDivisionRule2_2{
    public:
    using RandomDirectionCentreBasedDivisionRule2_2::RandomDirectionCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<2, 2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_,
            RandomDirectionCentreBasedDivisionRule2_2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_RandomDirectionCentreBasedDivisionRule2_2_class(py::module &m){
py::class_<RandomDirectionCentreBasedDivisionRule2_2 , RandomDirectionCentreBasedDivisionRule2_2_Overrides , boost::shared_ptr<RandomDirectionCentreBasedDivisionRule2_2 >  , AbstractCentreBasedDivisionRule<2>  >(m, "RandomDirectionCentreBasedDivisionRule2_2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>>(RandomDirectionCentreBasedDivisionRule2_2::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<2, 2> &)) &RandomDirectionCentreBasedDivisionRule2_2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
