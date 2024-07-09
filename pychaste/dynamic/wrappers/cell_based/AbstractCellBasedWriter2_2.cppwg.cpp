#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedWriter.hpp"

#include "AbstractCellBasedWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedWriter<2,2 > AbstractCellBasedWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellBasedWriter2_2_Overrides : public AbstractCellBasedWriter2_2{
    public:
    using AbstractCellBasedWriter2_2::AbstractCellBasedWriter;
    void OpenOutputFile(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter2_2,
            OpenOutputFile,
                    rOutputFileHandler);
    }
    void WriteTimeStamp() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter2_2,
            WriteTimeStamp,
            );
    }
    void WriteNewline() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter2_2,
            WriteNewline,
            );
    }

};
void register_AbstractCellBasedWriter2_2_class(py::module &m){
py::class_<AbstractCellBasedWriter2_2 , AbstractCellBasedWriter2_2_Overrides , boost::shared_ptr<AbstractCellBasedWriter2_2 >   >(m, "AbstractCellBasedWriter2_2")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "CloseFile",
            (void(AbstractCellBasedWriter2_2::*)()) &AbstractCellBasedWriter2_2::CloseFile,
            " "  )
        .def(
            "OpenOutputFile",
            (void(AbstractCellBasedWriter2_2::*)(::OutputFileHandler &)) &AbstractCellBasedWriter2_2::OpenOutputFile,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "OpenOutputFileForAppend",
            (void(AbstractCellBasedWriter2_2::*)(::OutputFileHandler &)) &AbstractCellBasedWriter2_2::OpenOutputFileForAppend,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "WriteTimeStamp",
            (void(AbstractCellBasedWriter2_2::*)()) &AbstractCellBasedWriter2_2::WriteTimeStamp,
            " "  )
        .def(
            "WriteNewline",
            (void(AbstractCellBasedWriter2_2::*)()) &AbstractCellBasedWriter2_2::WriteNewline,
            " "  )
        .def(
            "SetFileName",
            (void(AbstractCellBasedWriter2_2::*)(::std::string)) &AbstractCellBasedWriter2_2::SetFileName,
            " " , py::arg("fileName") )
        .def(
            "GetFileName",
            (::std::string(AbstractCellBasedWriter2_2::*)()) &AbstractCellBasedWriter2_2::GetFileName,
            " "  )
    ;
}
