#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChasteCuboid.hpp"

#include "ChasteCuboid3.cppwg.hpp"

namespace py = pybind11;
typedef ChasteCuboid<3 > ChasteCuboid3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChasteCuboid3_Overrides : public ChasteCuboid3{
    public:
    using ChasteCuboid3::ChasteCuboid;
    bool DoesContain(::ChastePoint<3> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE(
            bool,
            ChasteCuboid3,
            DoesContain,
                    rPointToCheck);
    }

};
void register_ChasteCuboid3_class(py::module &m){
py::class_<ChasteCuboid3 , ChasteCuboid3_Overrides , boost::shared_ptr<ChasteCuboid3 >  , AbstractChasteRegion<3>  >(m, "ChasteCuboid3")
        .def(py::init<::ChastePoint<3> &, ::ChastePoint<3> & >(), py::arg("rLowerPoint"), py::arg("rUpperPoint"))
        .def(
            "DoesContain",
            (bool(ChasteCuboid3::*)(::ChastePoint<3> const &) const ) &ChasteCuboid3::DoesContain,
            " " , py::arg("rPointToCheck") )
        .def(
            "rGetUpperCorner",
            (::ChastePoint<3> const &(ChasteCuboid3::*)() const ) &ChasteCuboid3::rGetUpperCorner,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetLowerCorner",
            (::ChastePoint<3> const &(ChasteCuboid3::*)() const ) &ChasteCuboid3::rGetLowerCorner,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetWidth",
            (double(ChasteCuboid3::*)(unsigned int) const ) &ChasteCuboid3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetLongestAxis",
            (unsigned int(ChasteCuboid3::*)() const ) &ChasteCuboid3::GetLongestAxis,
            " "  )
    ;
}
