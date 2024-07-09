#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VtkScene.hpp"

#include "VtkScene2.cppwg.hpp"

namespace py = pybind11;
typedef VtkScene<2 > VtkScene2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VtkScene2_Overrides : public VtkScene2{
    public:
    using VtkScene2::VtkScene;
    void ResetRenderer(unsigned int timeStep) override {
        PYBIND11_OVERRIDE(
            void,
            VtkScene2,
            ResetRenderer,
                    timeStep);
    }

};
void register_VtkScene2_class(py::module &m){
py::class_<VtkScene2 , VtkScene2_Overrides , boost::shared_ptr<VtkScene2 >   >(m, "VtkScene2")
        .def(py::init< >())
        .def(
            "End",
            (void(VtkScene2::*)()) &VtkScene2::End,
            " "  )
        .def(
            "GetSceneAsCharBuffer",
            (::vtkSmartPointer<vtkUnsignedCharArray>(VtkScene2::*)()) &VtkScene2::GetSceneAsCharBuffer,
            " "  )
        .def(
            "GetRenderer",
            (::vtkSmartPointer<vtkRenderer>(VtkScene2::*)()) &VtkScene2::GetRenderer,
            " "  )
        .def(
            "GetCellPopulationActorGenerator",
            (::boost::shared_ptr<CellPopulationPyChasteActorGenerator<2>>(VtkScene2::*)()) &VtkScene2::GetCellPopulationActorGenerator,
            " "  )
        .def(
            "ResetRenderer",
            (void(VtkScene2::*)(unsigned int)) &VtkScene2::ResetRenderer,
            " " , py::arg("timeStep") = 0 )
        .def(
            "Start",
            (void(VtkScene2::*)()) &VtkScene2::Start,
            " "  )
        .def(
            "SetCellPopulation",
            (void(VtkScene2::*)(::boost::shared_ptr<AbstractCellPopulation<2, 2>>)) &VtkScene2::SetCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetOutputFilePath",
            (void(VtkScene2::*)(::std::string const &)) &VtkScene2::SetOutputFilePath,
            " " , py::arg("rPath") )
        .def(
            "SetIsInteractive",
            (void(VtkScene2::*)(bool)) &VtkScene2::SetIsInteractive,
            " " , py::arg("isInteractive") )
        .def(
            "SetSaveAsAnimation",
            (void(VtkScene2::*)(bool)) &VtkScene2::SetSaveAsAnimation,
            " " , py::arg("saveAsAnimation") )
        .def(
            "SetSaveAsImages",
            (void(VtkScene2::*)(bool)) &VtkScene2::SetSaveAsImages,
            " " , py::arg("saveAsImages") )
        .def(
            "StartInteractiveEventHandler",
            (void(VtkScene2::*)()) &VtkScene2::StartInteractiveEventHandler,
            " "  )
    ;
}
