import cppwg.templates.custom


class MutableElementCustomTemplate(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_def_code(self, class_name):
        """
        Adds MutableElement constructor not implemented in 1D template.
        """

        ELEMENT_DIM, SPACE_DIM = class_name.replace("MutableElement", "").split("_")

        if ELEMENT_DIM == "1":
            return ""

        code = f"""\
        .def(py::init<unsigned int, ::std::vector<Node<{SPACE_DIM}> *> const & >(), py::arg("index"), py::arg("rNodes"))
"""
        return code
