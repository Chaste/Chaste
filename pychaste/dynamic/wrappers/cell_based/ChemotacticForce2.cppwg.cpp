#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChemotacticForce.hpp"

#include "ChemotacticForce2.cppwg.hpp"

namespace py = pybind11;
typedef ChemotacticForce<2 > ChemotacticForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChemotacticForce2_Overrides : public ChemotacticForce2{
    public:
    using ChemotacticForce2::ChemotacticForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotacticForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotacticForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_ChemotacticForce2_class(py::module &m){
py::class_<ChemotacticForce2 , ChemotacticForce2_Overrides , boost::shared_ptr<ChemotacticForce2 >  , AbstractForce<2>  >(m, "ChemotacticForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(ChemotacticForce2::*)(::AbstractCellPopulation<2> &)) &ChemotacticForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(ChemotacticForce2::*)(::out_stream &)) &ChemotacticForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
