#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractImmersedBoundaryForce.hpp"

#include "AbstractImmersedBoundaryForce2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractImmersedBoundaryForce<2 > AbstractImmersedBoundaryForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractImmersedBoundaryForce2_Overrides : public AbstractImmersedBoundaryForce2{
    public:
    using AbstractImmersedBoundaryForce2::AbstractImmersedBoundaryForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractImmersedBoundaryForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractImmersedBoundaryForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_AbstractImmersedBoundaryForce2_class(py::module &m){
py::class_<AbstractImmersedBoundaryForce2 , AbstractImmersedBoundaryForce2_Overrides , boost::shared_ptr<AbstractImmersedBoundaryForce2 >   >(m, "AbstractImmersedBoundaryForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(AbstractImmersedBoundaryForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &AbstractImmersedBoundaryForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(AbstractImmersedBoundaryForce2::*)(::out_stream &)) &AbstractImmersedBoundaryForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetAdditiveNormalNoise",
            (bool(AbstractImmersedBoundaryForce2::*)() const ) &AbstractImmersedBoundaryForce2::GetAdditiveNormalNoise,
            " "  )
        .def(
            "SetAdditiveNormalNoise",
            (void(AbstractImmersedBoundaryForce2::*)(bool)) &AbstractImmersedBoundaryForce2::SetAdditiveNormalNoise,
            " " , py::arg("additiveNormalNoise") )
        .def(
            "GetNormalNoiseMean",
            (double(AbstractImmersedBoundaryForce2::*)() const ) &AbstractImmersedBoundaryForce2::GetNormalNoiseMean,
            " "  )
        .def(
            "SetNormalNoiseMean",
            (void(AbstractImmersedBoundaryForce2::*)(double)) &AbstractImmersedBoundaryForce2::SetNormalNoiseMean,
            " " , py::arg("normalNoiseMean") )
        .def(
            "GetNormalNoiseStdDev",
            (double(AbstractImmersedBoundaryForce2::*)() const ) &AbstractImmersedBoundaryForce2::GetNormalNoiseStdDev,
            " "  )
        .def(
            "SetNormalNoiseStdDev",
            (void(AbstractImmersedBoundaryForce2::*)(double)) &AbstractImmersedBoundaryForce2::SetNormalNoiseStdDev,
            " " , py::arg("normalNoiseStdDev") )
    ;
}
