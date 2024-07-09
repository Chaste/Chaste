#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"

#include "RandomDirectionVertexBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef RandomDirectionVertexBasedDivisionRule<3 > RandomDirectionVertexBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class RandomDirectionVertexBasedDivisionRule3_Overrides : public RandomDirectionVertexBasedDivisionRule3{
    public:
    using RandomDirectionVertexBasedDivisionRule3::RandomDirectionVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            RandomDirectionVertexBasedDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_RandomDirectionVertexBasedDivisionRule3_class(py::module &m){
py::class_<RandomDirectionVertexBasedDivisionRule3 , RandomDirectionVertexBasedDivisionRule3_Overrides , boost::shared_ptr<RandomDirectionVertexBasedDivisionRule3 >  , AbstractVertexBasedDivisionRule<3>  >(m, "RandomDirectionVertexBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(RandomDirectionVertexBasedDivisionRule3::*)(::CellPtr, ::VertexBasedCellPopulation<3> &)) &RandomDirectionVertexBasedDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
