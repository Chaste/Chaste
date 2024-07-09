#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeElasticForce.hpp"

#include "BuskeElasticForce2.cppwg.hpp"

namespace py = pybind11;
typedef BuskeElasticForce<2 > BuskeElasticForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class BuskeElasticForce2_Overrides : public BuskeElasticForce2{
    public:
    using BuskeElasticForce2::BuskeElasticForce;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            BuskeElasticForce2,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeElasticForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeElasticForce2_class(py::module &m){
py::class_<BuskeElasticForce2 , BuskeElasticForce2_Overrides , boost::shared_ptr<BuskeElasticForce2 >  , AbstractTwoBodyInteractionForce<2>  >(m, "BuskeElasticForce2")
        .def(py::init< >())
        .def(
            "GetDeformationEnergyParameter",
            (double(BuskeElasticForce2::*)()) &BuskeElasticForce2::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(BuskeElasticForce2::*)(double)) &BuskeElasticForce2::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 2>(BuskeElasticForce2::*)(unsigned int, unsigned int, ::AbstractCellPopulation<2> &)) &BuskeElasticForce2::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "GetMagnitudeOfForce",
            (double(BuskeElasticForce2::*)(double, double, double)) &BuskeElasticForce2::GetMagnitudeOfForce,
            " " , py::arg("distanceBetweenNodes"), py::arg("radiusOfCellOne"), py::arg("radiusOfCellTwo") )
        .def(
            "OutputForceParameters",
            (void(BuskeElasticForce2::*)(::out_stream &)) &BuskeElasticForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
