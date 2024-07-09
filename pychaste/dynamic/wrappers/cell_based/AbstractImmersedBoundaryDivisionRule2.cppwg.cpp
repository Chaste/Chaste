#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractImmersedBoundaryDivisionRule.hpp"

#include "AbstractImmersedBoundaryDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractImmersedBoundaryDivisionRule<2 > AbstractImmersedBoundaryDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class AbstractImmersedBoundaryDivisionRule2_Overrides : public AbstractImmersedBoundaryDivisionRule2{
    public:
    using AbstractImmersedBoundaryDivisionRule2::AbstractImmersedBoundaryDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractImmersedBoundaryDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellImmersedBoundaryDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractImmersedBoundaryDivisionRule2,
            OutputCellImmersedBoundaryDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractImmersedBoundaryDivisionRule2_class(py::module &m){
py::class_<AbstractImmersedBoundaryDivisionRule2 , AbstractImmersedBoundaryDivisionRule2_Overrides , boost::shared_ptr<AbstractImmersedBoundaryDivisionRule2 >   >(m, "AbstractImmersedBoundaryDivisionRule2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractImmersedBoundaryDivisionRule2::*)(::CellPtr, ::ImmersedBoundaryCellPopulation<2> &)) &AbstractImmersedBoundaryDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellImmersedBoundaryDivisionRuleInfo",
            (void(AbstractImmersedBoundaryDivisionRule2::*)(::out_stream &)) &AbstractImmersedBoundaryDivisionRule2::OutputCellImmersedBoundaryDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
