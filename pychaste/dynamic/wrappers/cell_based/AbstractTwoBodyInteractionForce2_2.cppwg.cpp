#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"

#include "AbstractTwoBodyInteractionForce2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTwoBodyInteractionForce<2,2 > AbstractTwoBodyInteractionForce2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class AbstractTwoBodyInteractionForce2_2_Overrides : public AbstractTwoBodyInteractionForce2_2{
    public:
    using AbstractTwoBodyInteractionForce2_2::AbstractTwoBodyInteractionForce;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractTwoBodyInteractionForce2_2,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce2_2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce2_2,
            OutputForceParameters,
                    rParamsFile);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce2_2,
            WriteDataToVisualizerSetupFile,
                    pVizSetupFile);
    }

};
void register_AbstractTwoBodyInteractionForce2_2_class(py::module &m){
py::class_<AbstractTwoBodyInteractionForce2_2 , AbstractTwoBodyInteractionForce2_2_Overrides , boost::shared_ptr<AbstractTwoBodyInteractionForce2_2 >  , AbstractForce<2>  >(m, "AbstractTwoBodyInteractionForce2_2")
        .def(
            "GetUseCutOffLength",
            (bool(AbstractTwoBodyInteractionForce2_2::*)()) &AbstractTwoBodyInteractionForce2_2::GetUseCutOffLength,
            " "  )
        .def(
            "SetCutOffLength",
            (void(AbstractTwoBodyInteractionForce2_2::*)(double)) &AbstractTwoBodyInteractionForce2_2::SetCutOffLength,
            " " , py::arg("cutOffLength") )
        .def(
            "GetCutOffLength",
            (double(AbstractTwoBodyInteractionForce2_2::*)()) &AbstractTwoBodyInteractionForce2_2::GetCutOffLength,
            " "  )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractTwoBodyInteractionForce2_2::*)(unsigned int, unsigned int, ::AbstractCellPopulation<2> &)) &AbstractTwoBodyInteractionForce2_2::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "AddForceContribution",
            (void(AbstractTwoBodyInteractionForce2_2::*)(::AbstractCellPopulation<2> &)) &AbstractTwoBodyInteractionForce2_2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(AbstractTwoBodyInteractionForce2_2::*)(::out_stream &)) &AbstractTwoBodyInteractionForce2_2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "WriteDataToVisualizerSetupFile",
            (void(AbstractTwoBodyInteractionForce2_2::*)(::out_stream &)) &AbstractTwoBodyInteractionForce2_2::WriteDataToVisualizerSetupFile,
            " " , py::arg("pVizSetupFile") )
    ;
}
