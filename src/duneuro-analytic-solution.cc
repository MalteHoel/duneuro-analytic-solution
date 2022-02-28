// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// Class for creating bindings for analytic MEG forward solution in sphere models.
// Inspired by duneuro-py
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/operators.h>                                           // include for easy binding of +=, *=, etc.
#include <dune/duneuro-analytic-solution/duneuro-analytic-solution.hh>                // include for analytic MEG solution in sphere models
#include <dune/common/fvector.hh>
#include <duneuro/common/dipole.hh>
#include <iostream>
#include <algorithm>

namespace py = pybind11;
using Scalar = double;
enum {dim = 3};
using CoordinateType = Dune::FieldVector<Scalar, dim>;
using Dipole = duneuro::Dipole<Scalar, dim>;

///////////////////////////////////////////////////////////
// Bindings for the AnalyticSolutionMEG class
///////////////////////////////////////////////////////////
void register_analytic_solution_meg(py::module& m) {
  py::class_<duneuro::AnalyticSolutionMEG<Scalar>>(m, "AnalyticSolutionMEG", "class implementing the analytic solution of the MEG forward problem in multilayer sphere models")
    .def(py::init<const CoordinateType&, Scalar>(), "create analytic solver using the sphere center and the scaling factor", py::arg("sphere_center"), py::arg("scaling_factor") = 1.0)
    .def("bind", &duneuro::AnalyticSolutionMEG<Scalar>::bind, "bind the dipole we want to solve for")
    .def("totalField", py::overload_cast<const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::totalField), "compute the total magnetic field vector at the specified position")
    .def("totalField", py::overload_cast<const CoordinateType&, const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::totalField), "compute the total magnetic field at the specified position in the specified direction")
    .def("primaryField", py::overload_cast<const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::primaryField), "compute the primary magnetic field vector at the specified position")
    .def("primaryField", py::overload_cast<const CoordinateType&, const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::primaryField), "compute the primary magnetic field at the specified position in the specified direction")
    .def("secondaryField", py::overload_cast<const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::secondaryField), "compute the secondary magnetic field vector at the specified position")
    .def("secondaryField", py::overload_cast<const CoordinateType&, const CoordinateType&>(&duneuro::AnalyticSolutionMEG<Scalar>::secondaryField), "compute the secondary magnetic field at the specified position in the specified direction")
    ; // end definition of class
} // end register_analytic_solution_meg

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Create bindings
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
PYBIND11_MODULE(duneuroAnalyticSolutionPy, m) {
  register_analytic_solution_meg(m);
}
