
import cppwg.templates.custom


class PopulationWriter_custom(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_def_code(self, class_name):

        code = ""
        population_writers = ["VoronoiDataWriter"]
        for eachWriter in population_writers:
            local_reps = {"class_name": class_name,
                          "writer": eachWriter}
            local_code = """\
        .def("AddPopulationWriter{writer}", &{class_name}::AddPopulationWriter<{writer}>)
"""
            code += local_code.format(**local_reps)

        cell_writers = ["CellLabelWriter"]
        for eachWriter in cell_writers:
            local_reps = {"class_name": class_name,
                          "writer": eachWriter}
            local_code = """\
        .def("AddCellWriter{writer}", &{class_name}::AddCellWriter<{writer}>)
"""
            code += local_code.format(**local_reps)
        return code
