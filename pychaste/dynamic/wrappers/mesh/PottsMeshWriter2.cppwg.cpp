#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PottsMeshWriter.hpp"

#include "PottsMeshWriter2.cppwg.hpp"

namespace py = pybind11;
typedef PottsMeshWriter<2 > PottsMeshWriter2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::vector<double> _std_vector_lt_double_gt_;
typedef ::ElementData _ElementData;

class PottsMeshWriter2_Overrides : public PottsMeshWriter2{
    public:
    using PottsMeshWriter2::PottsMeshWriter;
    ::std::vector<double> GetNextNode() override {
        PYBIND11_OVERRIDE(
            _std_vector_lt_double_gt_,
            PottsMeshWriter2,
            GetNextNode,
            );
    }
    ::ElementData GetNextElement() override {
        PYBIND11_OVERRIDE(
            _ElementData,
            PottsMeshWriter2,
            GetNextElement,
            );
    }
    void WriteFiles() override {
        PYBIND11_OVERRIDE(
            void,
            PottsMeshWriter2,
            WriteFiles,
            );
    }

};
void register_PottsMeshWriter2_class(py::module &m){
py::class_<PottsMeshWriter2 , PottsMeshWriter2_Overrides , boost::shared_ptr<PottsMeshWriter2 >   >(m, "PottsMeshWriter2")
        .def(py::init<::std::string const &, ::std::string const &, bool const >(), py::arg("rDirectory"), py::arg("rBaseName"), py::arg("clearOutputDir") = true)
        .def(
            "WriteFilesUsingMesh",
            (void(PottsMeshWriter2::*)(::PottsMesh<2> &)) &PottsMeshWriter2::WriteFilesUsingMesh,
            " " , py::arg("rMesh") )
        .def(
            "GetNextNode",
            (::std::vector<double>(PottsMeshWriter2::*)()) &PottsMeshWriter2::GetNextNode,
            " "  )
        .def(
            "GetNextElement",
            (::ElementData(PottsMeshWriter2::*)()) &PottsMeshWriter2::GetNextElement,
            " "  )
        .def(
            "WriteFiles",
            (void(PottsMeshWriter2::*)()) &PottsMeshWriter2::WriteFiles,
            " "  )
    ;
}
