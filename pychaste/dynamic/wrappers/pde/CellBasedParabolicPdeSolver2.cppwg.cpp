#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellBasedParabolicPdeSolver.hpp"

#include "CellBasedParabolicPdeSolver2.cppwg.hpp"

namespace py = pybind11;
typedef CellBasedParabolicPdeSolver<2 > CellBasedParabolicPdeSolver2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class CellBasedParabolicPdeSolver2_Overrides : public CellBasedParabolicPdeSolver2{
    public:
    using CellBasedParabolicPdeSolver2::CellBasedParabolicPdeSolver;
    ::boost::numeric::ublas::c_vector<double, 3> ComputeVectorTerm(::boost::numeric::ublas::c_vector<double, 3> & rPhi, ::boost::numeric::ublas::c_matrix<double, 2, 3> & rGradPhi, ::ChastePoint<2> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 2> & rGradU, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            CellBasedParabolicPdeSolver2,
            ComputeVectorTerm,
                    rPhi,
        rGradPhi,
        rX,
        rU,
        rGradU,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeMatrixTerm(::boost::numeric::ublas::c_vector<double, 3> & rPhi, ::boost::numeric::ublas::c_matrix<double, 2, 3> & rGradPhi, ::ChastePoint<2> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 2> & rGradU, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            CellBasedParabolicPdeSolver2,
            ComputeMatrixTerm,
                    rPhi,
        rGradPhi,
        rX,
        rU,
        rGradU,
        pElement);
    }
    void ResetInterpolatedQuantities() override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedParabolicPdeSolver2,
            ResetInterpolatedQuantities,
            );
    }
    void IncrementInterpolatedQuantities(double phiI, ::Node<2> const * arg1) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedParabolicPdeSolver2,
            IncrementInterpolatedQuantities,
                    phiI,
        arg1);
    }

};
void register_CellBasedParabolicPdeSolver2_class(py::module &m){
py::class_<CellBasedParabolicPdeSolver2 , CellBasedParabolicPdeSolver2_Overrides , boost::shared_ptr<CellBasedParabolicPdeSolver2 >   >(m, "CellBasedParabolicPdeSolver2")
        .def(py::init<::TetrahedralMesh<2, 2> *, ::AbstractLinearParabolicPde<2> *, ::BoundaryConditionsContainer<2, 2, 1> * >(), py::arg("pMesh"), py::arg("pPde"), py::arg("pBoundaryConditions"))
    ;
}
