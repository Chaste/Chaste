import cppwg.templates.custom


class PopulationWriterCustomTemplate(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_def_code(self, class_name):
        """
        Creates custom wrapper code for adding population and cell writers.

        Parameters
        ----------
        class_name : str
            The name of the class being wrapped.

        Returns
        -------
        str
            The custom code to add to the class wrappers.
        """
        population_writers = ["VoronoiDataWriter"]

        cell_writers = ["CellLabelWriter"]

        population_writer_template = """\
        .def("AddPopulationWriter{writer}", &{class_name}::AddPopulationWriter<{writer}>)
"""

        cell_writer_template = """\
        .def("AddCellWriter{writer}", &{class_name}::AddCellWriter<{writer}>)
"""

        code = ""

        for writer in population_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += population_writer_template.format(**replacements)

        for writer in cell_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += cell_writer_template.format(**replacements)

        return code
