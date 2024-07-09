#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce.hpp"

#include "ImmersedBoundaryKinematicFeedbackForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryKinematicFeedbackForce<3 > ImmersedBoundaryKinematicFeedbackForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryKinematicFeedbackForce3_Overrides : public ImmersedBoundaryKinematicFeedbackForce3{
    public:
    using ImmersedBoundaryKinematicFeedbackForce3::ImmersedBoundaryKinematicFeedbackForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryKinematicFeedbackForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryKinematicFeedbackForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryKinematicFeedbackForce3_class(py::module &m){
py::class_<ImmersedBoundaryKinematicFeedbackForce3 , ImmersedBoundaryKinematicFeedbackForce3_Overrides , boost::shared_ptr<ImmersedBoundaryKinematicFeedbackForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryKinematicFeedbackForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryKinematicFeedbackForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryKinematicFeedbackForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryKinematicFeedbackForce3::*)(::out_stream &)) &ImmersedBoundaryKinematicFeedbackForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetSpringConst",
            (double(ImmersedBoundaryKinematicFeedbackForce3::*)() const ) &ImmersedBoundaryKinematicFeedbackForce3::GetSpringConst,
            " "  )
        .def(
            "SetSpringConst",
            (void(ImmersedBoundaryKinematicFeedbackForce3::*)(double)) &ImmersedBoundaryKinematicFeedbackForce3::SetSpringConst,
            " " , py::arg("springConst") )
    ;
}
