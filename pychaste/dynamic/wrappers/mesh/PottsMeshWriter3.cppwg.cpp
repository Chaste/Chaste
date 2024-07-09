#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PottsMeshWriter.hpp"

#include "PottsMeshWriter3.cppwg.hpp"

namespace py = pybind11;
typedef PottsMeshWriter<3 > PottsMeshWriter3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::vector<double> _std_vector_lt_double_gt_;
typedef ::ElementData _ElementData;

class PottsMeshWriter3_Overrides : public PottsMeshWriter3{
    public:
    using PottsMeshWriter3::PottsMeshWriter;
    ::std::vector<double> GetNextNode() override {
        PYBIND11_OVERRIDE(
            _std_vector_lt_double_gt_,
            PottsMeshWriter3,
            GetNextNode,
            );
    }
    ::ElementData GetNextElement() override {
        PYBIND11_OVERRIDE(
            _ElementData,
            PottsMeshWriter3,
            GetNextElement,
            );
    }
    void WriteFiles() override {
        PYBIND11_OVERRIDE(
            void,
            PottsMeshWriter3,
            WriteFiles,
            );
    }

};
void register_PottsMeshWriter3_class(py::module &m){
py::class_<PottsMeshWriter3 , PottsMeshWriter3_Overrides , boost::shared_ptr<PottsMeshWriter3 >   >(m, "PottsMeshWriter3")
        .def(py::init<::std::string const &, ::std::string const &, bool const >(), py::arg("rDirectory"), py::arg("rBaseName"), py::arg("clearOutputDir") = true)
        .def(
            "WriteFilesUsingMesh",
            (void(PottsMeshWriter3::*)(::PottsMesh<3> &)) &PottsMeshWriter3::WriteFilesUsingMesh,
            " " , py::arg("rMesh") )
        .def(
            "GetNextNode",
            (::std::vector<double>(PottsMeshWriter3::*)()) &PottsMeshWriter3::GetNextNode,
            " "  )
        .def(
            "GetNextElement",
            (::ElementData(PottsMeshWriter3::*)()) &PottsMeshWriter3::GetNextElement,
            " "  )
        .def(
            "WriteFiles",
            (void(PottsMeshWriter3::*)()) &PottsMeshWriter3::WriteFiles,
            " "  )
    ;
}
