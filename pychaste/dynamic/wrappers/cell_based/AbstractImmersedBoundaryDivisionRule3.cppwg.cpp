#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractImmersedBoundaryDivisionRule.hpp"

#include "AbstractImmersedBoundaryDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractImmersedBoundaryDivisionRule<3 > AbstractImmersedBoundaryDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class AbstractImmersedBoundaryDivisionRule3_Overrides : public AbstractImmersedBoundaryDivisionRule3{
    public:
    using AbstractImmersedBoundaryDivisionRule3::AbstractImmersedBoundaryDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractImmersedBoundaryDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellImmersedBoundaryDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractImmersedBoundaryDivisionRule3,
            OutputCellImmersedBoundaryDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractImmersedBoundaryDivisionRule3_class(py::module &m){
py::class_<AbstractImmersedBoundaryDivisionRule3 , AbstractImmersedBoundaryDivisionRule3_Overrides , boost::shared_ptr<AbstractImmersedBoundaryDivisionRule3 >   >(m, "AbstractImmersedBoundaryDivisionRule3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractImmersedBoundaryDivisionRule3::*)(::CellPtr, ::ImmersedBoundaryCellPopulation<3> &)) &AbstractImmersedBoundaryDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellImmersedBoundaryDivisionRuleInfo",
            (void(AbstractImmersedBoundaryDivisionRule3::*)(::out_stream &)) &AbstractImmersedBoundaryDivisionRule3::OutputCellImmersedBoundaryDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
