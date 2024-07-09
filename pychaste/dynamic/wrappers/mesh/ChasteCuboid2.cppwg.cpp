#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChasteCuboid.hpp"

#include "ChasteCuboid2.cppwg.hpp"

namespace py = pybind11;
typedef ChasteCuboid<2 > ChasteCuboid2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChasteCuboid2_Overrides : public ChasteCuboid2{
    public:
    using ChasteCuboid2::ChasteCuboid;
    bool DoesContain(::ChastePoint<2> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE(
            bool,
            ChasteCuboid2,
            DoesContain,
                    rPointToCheck);
    }

};
void register_ChasteCuboid2_class(py::module &m){
py::class_<ChasteCuboid2 , ChasteCuboid2_Overrides , boost::shared_ptr<ChasteCuboid2 >  , AbstractChasteRegion<2>  >(m, "ChasteCuboid2")
        .def(py::init<::ChastePoint<2> &, ::ChastePoint<2> & >(), py::arg("rLowerPoint"), py::arg("rUpperPoint"))
        .def(
            "DoesContain",
            (bool(ChasteCuboid2::*)(::ChastePoint<2> const &) const ) &ChasteCuboid2::DoesContain,
            " " , py::arg("rPointToCheck") )
        .def(
            "rGetUpperCorner",
            (::ChastePoint<2> const &(ChasteCuboid2::*)() const ) &ChasteCuboid2::rGetUpperCorner,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetLowerCorner",
            (::ChastePoint<2> const &(ChasteCuboid2::*)() const ) &ChasteCuboid2::rGetLowerCorner,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetWidth",
            (double(ChasteCuboid2::*)(unsigned int) const ) &ChasteCuboid2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetLongestAxis",
            (unsigned int(ChasteCuboid2::*)() const ) &ChasteCuboid2::GetLongestAxis,
            " "  )
    ;
}
