#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"

#include "RandomDirectionCentreBasedDivisionRule3_3.cppwg.hpp"

namespace py = pybind11;
typedef RandomDirectionCentreBasedDivisionRule<3,3 > RandomDirectionCentreBasedDivisionRule3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_;

class RandomDirectionCentreBasedDivisionRule3_3_Overrides : public RandomDirectionCentreBasedDivisionRule3_3{
    public:
    using RandomDirectionCentreBasedDivisionRule3_3::RandomDirectionCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<3, 3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_,
            RandomDirectionCentreBasedDivisionRule3_3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_RandomDirectionCentreBasedDivisionRule3_3_class(py::module &m){
py::class_<RandomDirectionCentreBasedDivisionRule3_3 , RandomDirectionCentreBasedDivisionRule3_3_Overrides , boost::shared_ptr<RandomDirectionCentreBasedDivisionRule3_3 >  , AbstractCentreBasedDivisionRule<3>  >(m, "RandomDirectionCentreBasedDivisionRule3_3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>>(RandomDirectionCentreBasedDivisionRule3_3::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<3, 3> &)) &RandomDirectionCentreBasedDivisionRule3_3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
