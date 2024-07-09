#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"

#include "AbstractTwoBodyInteractionForce3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTwoBodyInteractionForce<3,3 > AbstractTwoBodyInteractionForce3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class AbstractTwoBodyInteractionForce3_3_Overrides : public AbstractTwoBodyInteractionForce3_3{
    public:
    using AbstractTwoBodyInteractionForce3_3::AbstractTwoBodyInteractionForce;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractTwoBodyInteractionForce3_3,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce3_3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce3_3,
            OutputForceParameters,
                    rParamsFile);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTwoBodyInteractionForce3_3,
            WriteDataToVisualizerSetupFile,
                    pVizSetupFile);
    }

};
void register_AbstractTwoBodyInteractionForce3_3_class(py::module &m){
py::class_<AbstractTwoBodyInteractionForce3_3 , AbstractTwoBodyInteractionForce3_3_Overrides , boost::shared_ptr<AbstractTwoBodyInteractionForce3_3 >  , AbstractForce<3>  >(m, "AbstractTwoBodyInteractionForce3_3")
        .def(
            "GetUseCutOffLength",
            (bool(AbstractTwoBodyInteractionForce3_3::*)()) &AbstractTwoBodyInteractionForce3_3::GetUseCutOffLength,
            " "  )
        .def(
            "SetCutOffLength",
            (void(AbstractTwoBodyInteractionForce3_3::*)(double)) &AbstractTwoBodyInteractionForce3_3::SetCutOffLength,
            " " , py::arg("cutOffLength") )
        .def(
            "GetCutOffLength",
            (double(AbstractTwoBodyInteractionForce3_3::*)()) &AbstractTwoBodyInteractionForce3_3::GetCutOffLength,
            " "  )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractTwoBodyInteractionForce3_3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &)) &AbstractTwoBodyInteractionForce3_3::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "AddForceContribution",
            (void(AbstractTwoBodyInteractionForce3_3::*)(::AbstractCellPopulation<3> &)) &AbstractTwoBodyInteractionForce3_3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(AbstractTwoBodyInteractionForce3_3::*)(::out_stream &)) &AbstractTwoBodyInteractionForce3_3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "WriteDataToVisualizerSetupFile",
            (void(AbstractTwoBodyInteractionForce3_3::*)(::out_stream &)) &AbstractTwoBodyInteractionForce3_3::WriteDataToVisualizerSetupFile,
            " " , py::arg("pVizSetupFile") )
    ;
}
