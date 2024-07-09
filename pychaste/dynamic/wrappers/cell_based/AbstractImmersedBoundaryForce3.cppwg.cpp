#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractImmersedBoundaryForce.hpp"

#include "AbstractImmersedBoundaryForce3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractImmersedBoundaryForce<3 > AbstractImmersedBoundaryForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractImmersedBoundaryForce3_Overrides : public AbstractImmersedBoundaryForce3{
    public:
    using AbstractImmersedBoundaryForce3::AbstractImmersedBoundaryForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractImmersedBoundaryForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractImmersedBoundaryForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_AbstractImmersedBoundaryForce3_class(py::module &m){
py::class_<AbstractImmersedBoundaryForce3 , AbstractImmersedBoundaryForce3_Overrides , boost::shared_ptr<AbstractImmersedBoundaryForce3 >   >(m, "AbstractImmersedBoundaryForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(AbstractImmersedBoundaryForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &AbstractImmersedBoundaryForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(AbstractImmersedBoundaryForce3::*)(::out_stream &)) &AbstractImmersedBoundaryForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetAdditiveNormalNoise",
            (bool(AbstractImmersedBoundaryForce3::*)() const ) &AbstractImmersedBoundaryForce3::GetAdditiveNormalNoise,
            " "  )
        .def(
            "SetAdditiveNormalNoise",
            (void(AbstractImmersedBoundaryForce3::*)(bool)) &AbstractImmersedBoundaryForce3::SetAdditiveNormalNoise,
            " " , py::arg("additiveNormalNoise") )
        .def(
            "GetNormalNoiseMean",
            (double(AbstractImmersedBoundaryForce3::*)() const ) &AbstractImmersedBoundaryForce3::GetNormalNoiseMean,
            " "  )
        .def(
            "SetNormalNoiseMean",
            (void(AbstractImmersedBoundaryForce3::*)(double)) &AbstractImmersedBoundaryForce3::SetNormalNoiseMean,
            " " , py::arg("normalNoiseMean") )
        .def(
            "GetNormalNoiseStdDev",
            (double(AbstractImmersedBoundaryForce3::*)() const ) &AbstractImmersedBoundaryForce3::GetNormalNoiseStdDev,
            " "  )
        .def(
            "SetNormalNoiseStdDev",
            (void(AbstractImmersedBoundaryForce3::*)(double)) &AbstractImmersedBoundaryForce3::SetNormalNoiseStdDev,
            " " , py::arg("normalNoiseStdDev") )
    ;
}
