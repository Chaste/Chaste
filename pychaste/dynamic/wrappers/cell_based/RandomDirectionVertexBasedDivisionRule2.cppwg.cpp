#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"

#include "RandomDirectionVertexBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef RandomDirectionVertexBasedDivisionRule<2 > RandomDirectionVertexBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class RandomDirectionVertexBasedDivisionRule2_Overrides : public RandomDirectionVertexBasedDivisionRule2{
    public:
    using RandomDirectionVertexBasedDivisionRule2::RandomDirectionVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            RandomDirectionVertexBasedDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_RandomDirectionVertexBasedDivisionRule2_class(py::module &m){
py::class_<RandomDirectionVertexBasedDivisionRule2 , RandomDirectionVertexBasedDivisionRule2_Overrides , boost::shared_ptr<RandomDirectionVertexBasedDivisionRule2 >  , AbstractVertexBasedDivisionRule<2>  >(m, "RandomDirectionVertexBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(RandomDirectionVertexBasedDivisionRule2::*)(::CellPtr, ::VertexBasedCellPopulation<2> &)) &RandomDirectionVertexBasedDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
