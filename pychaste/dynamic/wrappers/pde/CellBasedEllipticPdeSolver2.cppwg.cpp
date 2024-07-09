#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellBasedEllipticPdeSolver.hpp"

#include "CellBasedEllipticPdeSolver2.cppwg.hpp"

namespace py = pybind11;
typedef CellBasedEllipticPdeSolver<2 > CellBasedEllipticPdeSolver2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class CellBasedEllipticPdeSolver2_Overrides : public CellBasedEllipticPdeSolver2{
    public:
    using CellBasedEllipticPdeSolver2::CellBasedEllipticPdeSolver;
    ::boost::numeric::ublas::c_vector<double, 3> ComputeVectorTerm(::boost::numeric::ublas::c_vector<double, 3> & rPhi, ::boost::numeric::ublas::c_matrix<double, 2, 3> & rGradPhi, ::ChastePoint<2> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 2> & rGradU, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            CellBasedEllipticPdeSolver2,
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
            CellBasedEllipticPdeSolver2,
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
            CellBasedEllipticPdeSolver2,
            ResetInterpolatedQuantities,
            );
    }
    void IncrementInterpolatedQuantities(double phiI, ::Node<2> const * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver2,
            IncrementInterpolatedQuantities,
                    phiI,
        pNode);
    }
    void InitialiseForSolve(::Vec initialSolution) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver2,
            InitialiseForSolve,
                    initialSolution);
    }

};
void register_CellBasedEllipticPdeSolver2_class(py::module &m){
py::class_<CellBasedEllipticPdeSolver2 , CellBasedEllipticPdeSolver2_Overrides , boost::shared_ptr<CellBasedEllipticPdeSolver2 >   >(m, "CellBasedEllipticPdeSolver2")
        .def(py::init<::TetrahedralMesh<2, 2> *, ::AbstractLinearEllipticPde<2, 2> *, ::BoundaryConditionsContainer<2, 2, 1> * >(), py::arg("pMesh"), py::arg("pPde"), py::arg("pBoundaryConditions"))
    ;
}
