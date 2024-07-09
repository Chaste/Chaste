#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce.hpp"

#include "ImmersedBoundaryKinematicFeedbackForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryKinematicFeedbackForce<2 > ImmersedBoundaryKinematicFeedbackForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryKinematicFeedbackForce2_Overrides : public ImmersedBoundaryKinematicFeedbackForce2{
    public:
    using ImmersedBoundaryKinematicFeedbackForce2::ImmersedBoundaryKinematicFeedbackForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryKinematicFeedbackForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryKinematicFeedbackForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryKinematicFeedbackForce2_class(py::module &m){
py::class_<ImmersedBoundaryKinematicFeedbackForce2 , ImmersedBoundaryKinematicFeedbackForce2_Overrides , boost::shared_ptr<ImmersedBoundaryKinematicFeedbackForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryKinematicFeedbackForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryKinematicFeedbackForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryKinematicFeedbackForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryKinematicFeedbackForce2::*)(::out_stream &)) &ImmersedBoundaryKinematicFeedbackForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetSpringConst",
            (double(ImmersedBoundaryKinematicFeedbackForce2::*)() const ) &ImmersedBoundaryKinematicFeedbackForce2::GetSpringConst,
            " "  )
        .def(
            "SetSpringConst",
            (void(ImmersedBoundaryKinematicFeedbackForce2::*)(double)) &ImmersedBoundaryKinematicFeedbackForce2::SetSpringConst,
            " " , py::arg("springConst") )
    ;
}
