#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VtkScene.hpp"

#include "VtkScene3.cppwg.hpp"

namespace py = pybind11;
typedef VtkScene<3 > VtkScene3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VtkScene3_Overrides : public VtkScene3{
    public:
    using VtkScene3::VtkScene;
    void ResetRenderer(unsigned int timeStep) override {
        PYBIND11_OVERRIDE(
            void,
            VtkScene3,
            ResetRenderer,
                    timeStep);
    }

};
void register_VtkScene3_class(py::module &m){
py::class_<VtkScene3 , VtkScene3_Overrides , boost::shared_ptr<VtkScene3 >   >(m, "VtkScene3")
        .def(py::init< >())
        .def(
            "End",
            (void(VtkScene3::*)()) &VtkScene3::End,
            " "  )
        .def(
            "GetSceneAsCharBuffer",
            (::vtkSmartPointer<vtkUnsignedCharArray>(VtkScene3::*)()) &VtkScene3::GetSceneAsCharBuffer,
            " "  )
        .def(
            "GetRenderer",
            (::vtkSmartPointer<vtkRenderer>(VtkScene3::*)()) &VtkScene3::GetRenderer,
            " "  )
        .def(
            "GetCellPopulationActorGenerator",
            (::boost::shared_ptr<CellPopulationPyChasteActorGenerator<3>>(VtkScene3::*)()) &VtkScene3::GetCellPopulationActorGenerator,
            " "  )
        .def(
            "ResetRenderer",
            (void(VtkScene3::*)(unsigned int)) &VtkScene3::ResetRenderer,
            " " , py::arg("timeStep") = 0 )
        .def(
            "Start",
            (void(VtkScene3::*)()) &VtkScene3::Start,
            " "  )
        .def(
            "SetCellPopulation",
            (void(VtkScene3::*)(::boost::shared_ptr<AbstractCellPopulation<3, 3>>)) &VtkScene3::SetCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetOutputFilePath",
            (void(VtkScene3::*)(::std::string const &)) &VtkScene3::SetOutputFilePath,
            " " , py::arg("rPath") )
        .def(
            "SetIsInteractive",
            (void(VtkScene3::*)(bool)) &VtkScene3::SetIsInteractive,
            " " , py::arg("isInteractive") )
        .def(
            "SetSaveAsAnimation",
            (void(VtkScene3::*)(bool)) &VtkScene3::SetSaveAsAnimation,
            " " , py::arg("saveAsAnimation") )
        .def(
            "SetSaveAsImages",
            (void(VtkScene3::*)(bool)) &VtkScene3::SetSaveAsImages,
            " " , py::arg("saveAsImages") )
        .def(
            "StartInteractiveEventHandler",
            (void(VtkScene3::*)()) &VtkScene3::StartInteractiveEventHandler,
            " "  )
    ;
}
