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

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>                                           // include for easy binding of +=, *=, etc.
#include <dune/duneuro-analytic-solution/duneuro-analytic-solution.hh>    // include for analytic MEG solution in sphere models
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
// Bindings for the FieldVector class
///////////////////////////////////////////////////////////

// this is a smaller version of the corresponding code in duneuro-py
void register_coordinate_vector(py::module& m) {
  py::class_<CoordinateType>(m, 
                         "Coordinate", 
                         "3-dimensional Vector containing the cartesian coordinates of a point",
                         py::buffer_protocol())
    .def_buffer([] (CoordinateType& vector) -> py::buffer_info {
        return py::buffer_info(
          &vector[0],                                        // pointer to buffer
          sizeof(Scalar),                                    // size of one entry in bytes
          py::format_descriptor<Scalar>::format(),           // python struct for data type
          1,                                                 // number of dimensions
          {dim},                                             // dimension sizes
          {sizeof(Scalar)}                                   // stride sizes in bytes
        ); // end buffer_info constructor
      } // end lambda definition
    ) // end def_buffer
    .def(py::init<Scalar>(), "create Coordinate from scalar")
    .def(py::init(
      [](py::buffer buffer) {
        // check buffer format  
        py::buffer_info info = buffer.request();
        if(info.format != py::format_descriptor<Scalar>::format()) {
          DUNE_THROW(Dune::Exception, "buffer entries are of the wrong type");
        }
        if(info.ndim != 1) {
          DUNE_THROW(Dune::Exception, "buffer has to consist of 1 dimension, but consists of " << info.ndim << " dimensions");
        }
        if(info.shape[0] != dim) {
          DUNE_THROW(Dune::Exception, "buffer has to contain 3 entries, but contains " << info.shape[0] << " entries");
        }
        
        Scalar* data_ptr = static_cast<Scalar*>(info.ptr);
        return CoordinateType({data_ptr[0], data_ptr[1], data_ptr[2]});
      }) // end definition of lambda
      , "create Coordinate from Python buffer"
    ) // end definition of py::init
    .def(py::init(
      [](const py::list& value_list) {

        // validate list
        if(value_list.size() != dim) {
          DUNE_THROW(Dune::Exception, "list has to contain 3 entries, but contains" << value_list.size() << " entries");
        }
        // list validated
        
        // copy values to coordinate vector
        CoordinateType coordinate;
        std::transform(value_list.begin(), value_list.end(), coordinate.begin(), [] (const py::handle& handle) -> Scalar {return handle.cast<Scalar>();});
        return coordinate;
      }), // end definition of lambda
      "create coordinate vector from list"
    ) // end definition of py::init
    .def("__len__", [] (const CoordinateType& coordinate) {return coordinate.size();})
    .def("__getitem__",
      [](const CoordinateType& coordinate, size_t index) {
        return coordinate[index];
      } // end definition of lambda
    ) // end definition of __getitem__
    .def("__setitem__",
      [] (CoordinateType& coordinate, size_t index, Scalar value) {
        coordinate[index] = value;
      } // end definition of lambda
    ) // end definition of __setitem_
    .def("__str__",
      [](const CoordinateType& coordinate) {
        std::stringstream sstr;
        sstr << " Coordinate with entries [" << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << "]";
        return sstr.str();
      } // end definition of lambda
    ) // end definition of __str__
    // bind arithmetic operations
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self += Scalar())
    .def(py::self -= Scalar())
    .def(py::self *= Scalar())
    .def(py::self /= Scalar())
  ; // end definition of class
} // end register_coordinate_vector


///////////////////////////////////////////////////////////
// Bindings for the Dipole class
///////////////////////////////////////////////////////////

void register_dipole(py::module& m) {
  py::class_<Dipole>(m, "Dipole", "Class representing a mathematical point dipole, consisting of a position vector and a moment vector")
    // define constructors
    .def(py::init<const CoordinateType&, const CoordinateType&>(), "create dipole from given position and moment", py::arg("position"), py::arg("moment"))
    .def(py::init(
      [] (py::buffer position_buffer, py::buffer moment_buffer) {

        // validate buffer
        py::buffer_info info_position_buffer = position_buffer.request();
        py::buffer_info info_moment_buffer = moment_buffer.request();
        
        if(info_position_buffer.format != py::format_descriptor<Scalar>::format() || info_moment_buffer.format != py::format_descriptor<Scalar>::format()) {
          DUNE_THROW(Dune::Exception, "buffer entries are of the wrong type");
        }
        if(info_position_buffer.ndim != 1 || info_moment_buffer.ndim != 1) {
          DUNE_THROW(Dune::Exception, "buffers have to consist of 1 dimension, but consist of " << info_position_buffer.ndim 
                                      << "and " << info_moment_buffer.ndim << " dimensions");
        }
        if(info_position_buffer.shape[0] != dim || info_moment_buffer.shape[0] != dim) {
          DUNE_THROW(Dune::Exception, "buffers have to contain 3 entries each, but contain " << info_position_buffer.shape[0] 
                                      << " and " << info_moment_buffer.shape[0] << " entries");
        }
        // buffer validated
        
        // copy data into Coordinate vectors
        CoordinateType position;
        CoordinateType moment;
        
        Scalar* position_ptr = static_cast<Scalar*>(info_position_buffer.ptr);
        Scalar* moment_ptr = static_cast<Scalar*>(info_moment_buffer.ptr);
        
        std::copy(position_ptr, position_ptr + dim, position.begin());
        std::copy(moment_ptr, moment_ptr + dim, moment.begin());
        
        return Dipole(position, moment);
      }), // end definition of lambda
      "create dipole from a position buffer and a moment buffer"
    ) // end definition of constructor
    .def(py::init(
      [] (py::buffer combined_position_moment_buffer) {
        
        // validate buffer
        py::buffer_info info_buffer = combined_position_moment_buffer.request();
        
        if(info_buffer.format != py::format_descriptor<Scalar>::format()) {
          DUNE_THROW(Dune::Exception, "buffer entries are of the wrong type");
        }
        if(info_buffer.ndim != 1) {
          DUNE_THROW(Dune::Exception, "buffer has to consist of 1 dimension, but consists of " << info_buffer.ndim << " dimensions");
        }
        if(info_buffer.shape[0] != 2 * dim) {
          DUNE_THROW(Dune::Exception, "buffer has to contain 6 entries, but contains " << info_buffer.shape[0] << " entries");
        }
        // buffer validated
        
        // copy data from buffer into position and moment vector
        CoordinateType position;
        CoordinateType moment;
        
        Scalar* data_ptr = static_cast<Scalar*>(info_buffer.ptr);
        std::copy(data_ptr, data_ptr + dim, position.begin());
        std::copy(data_ptr + dim, data_ptr + 2 * dim, moment.begin());
        
        return Dipole(position, moment);
      }), // end definition of lambda
      "create dipole from a single buffer containing a position and a moment"
    ) // end definition of constructor
    .def(py::init(
      [] (const py::list& pos_list, const py::list& mom_list) {
        
        // validate lists
        if(pos_list.size() != dim || mom_list.size() != dim) {
          DUNE_THROW(Dune::Exception, "lists have to be of size 3, but are of sizes " << pos_list.size() << " and " << mom_list.size());
        }
        // lists validated
        
        // copy lists to vectors
        CoordinateType position;
        CoordinateType moment;
        
        std::transform(pos_list.begin(), pos_list.end(), position.begin(), [] (const py::handle& handle) -> Scalar {return handle.cast<Scalar>();});
        std::transform(mom_list.begin(), mom_list.end(), moment.begin(), [] (const py::handle& handle) -> Scalar {return handle.cast<Scalar>();});
        
        return Dipole(position, moment);
      }), // end definition of lambda
      "create dipole from two lists containing the position and the moment"
    ) // end definition of constructor
    .def(py::init(
      [] (const py::list& combined_list) {

        // validate list
        if(combined_list.size() != 2 * dim) {
          DUNE_THROW(Dune::Exception, "list has to be of size 6, but is of size " << combined_list.size());
        }
        // list validated
        
        // copy list to vectors
        CoordinateType position;
        CoordinateType moment;
        
        auto pos_iterator = position.begin();
        auto moment_iterator = moment.begin();
        auto pos_end = position.end();
        for(const py::handle& handle : combined_list) {
          if(pos_iterator != pos_end) {
            *pos_iterator = handle.cast<Scalar>();
            ++pos_iterator;
          }
          else {
            *moment_iterator = handle.cast<Scalar>();
            ++moment_iterator;
          }
        }
        
        return Dipole(position, moment);
      }), // end definition of lambda
      "create dipole from a single list containing the position and the moment"
    ) // end definition of constructor
    .def("position", &Dipole::position, "position of the dipole", py::return_value_policy::reference_internal)
    .def("moment", &Dipole::moment, "moment of the dipole", py::return_value_policy::reference_internal)
    .def("__str__",
      [] (const Dipole& dipole) {
        const auto& pos = dipole.position();
        const auto& mom = dipole.moment();
        std::stringstream sstr;
        sstr << "Dipole with position [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "] and moment ["
             << mom[0] << ", " << mom[1] << ", " << mom[2] << "]";
        return sstr.str();
      } // end definition of lambda
    ) // end definitio of __str__
  ; // end definition of class
} // end register_dipole

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
  register_coordinate_vector(m);
  register_dipole(m);
  register_analytic_solution_meg(m);
}
