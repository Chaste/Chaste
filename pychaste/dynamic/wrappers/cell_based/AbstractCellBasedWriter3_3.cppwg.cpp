#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedWriter.hpp"

#include "AbstractCellBasedWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedWriter<3,3 > AbstractCellBasedWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellBasedWriter3_3_Overrides : public AbstractCellBasedWriter3_3{
    public:
    using AbstractCellBasedWriter3_3::AbstractCellBasedWriter;
    void OpenOutputFile(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter3_3,
            OpenOutputFile,
                    rOutputFileHandler);
    }
    void WriteTimeStamp() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter3_3,
            WriteTimeStamp,
            );
    }
    void WriteNewline() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter3_3,
            WriteNewline,
            );
    }

};
void register_AbstractCellBasedWriter3_3_class(py::module &m){
py::class_<AbstractCellBasedWriter3_3 , AbstractCellBasedWriter3_3_Overrides , boost::shared_ptr<AbstractCellBasedWriter3_3 >   >(m, "AbstractCellBasedWriter3_3")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "CloseFile",
            (void(AbstractCellBasedWriter3_3::*)()) &AbstractCellBasedWriter3_3::CloseFile,
            " "  )
        .def(
            "OpenOutputFile",
            (void(AbstractCellBasedWriter3_3::*)(::OutputFileHandler &)) &AbstractCellBasedWriter3_3::OpenOutputFile,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "OpenOutputFileForAppend",
            (void(AbstractCellBasedWriter3_3::*)(::OutputFileHandler &)) &AbstractCellBasedWriter3_3::OpenOutputFileForAppend,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "WriteTimeStamp",
            (void(AbstractCellBasedWriter3_3::*)()) &AbstractCellBasedWriter3_3::WriteTimeStamp,
            " "  )
        .def(
            "WriteNewline",
            (void(AbstractCellBasedWriter3_3::*)()) &AbstractCellBasedWriter3_3::WriteNewline,
            " "  )
        .def(
            "SetFileName",
            (void(AbstractCellBasedWriter3_3::*)(::std::string)) &AbstractCellBasedWriter3_3::SetFileName,
            " " , py::arg("fileName") )
        .def(
            "GetFileName",
            (::std::string(AbstractCellBasedWriter3_3::*)()) &AbstractCellBasedWriter3_3::GetFileName,
            " "  )
    ;
}
