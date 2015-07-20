// ============================================================================
// A simple python interface for Torus code
// ============================================================================

#include <Python.h>
#include <iostream>
#include <fstream>
#include <utils/Vector.h>
#include <string>
#include <Torus.h>
#include "PJM_utils.h"
#include "pot/falPot.h"
#include "pot/LogPot.h"
#include "pot/MiyamotoNagaiPot.h"
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/list.hpp>
#include <numpy/arrayobject.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

using namespace boost::python;

// ==========================================================================
// Define converters for Vector <==> numpy array

template <class T, int N>
struct WDvector_to_ndarray {
  static PyObject* convert(const Vector<T, N> v) {
    list l;
    for (int p = 0; p < N; ++p) {
      l.append(object(v[p]));
    }
    return incref(numeric::array(l).ptr());
  }
};

template <typename T, int N>
struct WDvector_from_ndarray {
  WDvector_from_ndarray() {
    converter::registry::push_back(&WDvector_from_ndarray<T, N>::convertible,
                                   &WDvector_from_ndarray<T, N>::construct,
                                   type_id<Vector<T, N>>());
  }

  // Determine if obj_ptr can be converted in a std::vector<T>
  static void* convertible(PyObject* obj_ptr) {
    if (!PyArray_Check(obj_ptr)) {
      std::cerr << "You have passed a non-numpy array" << std::endl;
      return 0;
    }
    return obj_ptr;
  }

  // Convert obj_ptr into a std::vector<T>
  static void construct(PyObject* obj_ptr,
                        converter::rvalue_from_python_stage1_data* data) {
    list l(handle<>(borrowed(obj_ptr)));
    // Grab pointer to memory into which to construct the new std::vector<T>
    void* storage = ((converter::rvalue_from_python_storage<Vector<T, N>>*)data)
                        ->storage.bytes;
    // in-place construct the new std::vector<T> using the character data
    // extraced from the python object
    Vector<T, N>& v = *(new (storage) Vector<T, N>());
    // populate the vector from list contains !!!
    for (int i = 0; i < N; ++i) {
      v[i] = extract<T>(l[i]);
    }
    // Stash the memory chunk pointer for later use by boost.python
    data->convertible = storage;
  }
};

// ==========================================================================
// Define converters for Vector <==> numpy array

struct PSPT_to_ndarray {
  static PyObject* convert(const PSPT v) {
    list l;
    for (int p = 0; p < 6; ++p) {
      l.append(object(v(p)));
    }
    return incref(numeric::array(l).ptr());
  }
};

struct PSPT_from_ndarray {
  PSPT_from_ndarray() {
    converter::registry::push_back(&PSPT_from_ndarray::convertible,
                                   &PSPT_from_ndarray::construct,
                                   type_id<PSPT>());
  }

  // Determine if obj_ptr can be converted in a std::vector<T>
  static void* convertible(PyObject* obj_ptr) {
    if (!PyArray_Check(obj_ptr)) {
      std::cerr << "You have passed a non-numpy array" << std::endl;
      return 0;
    }
    return obj_ptr;
  }

  // Convert obj_ptr into a std::vector<T>
  static void construct(PyObject* obj_ptr,
                        converter::rvalue_from_python_stage1_data* data) {
    list l(handle<>(borrowed(obj_ptr)));
    void* storage =
        ((converter::rvalue_from_python_storage<PSPT>*)data)->storage.bytes;
    PSPT& v = *(new (storage) PSPT);
    for (int i = 0; i < 6; ++i) v[i] = extract<double>(l[i]);
    data->convertible = storage;
  }
};

// ==========================================================================
// For creating potentials in python

struct Potential_PythonCallback : Potential {
  Potential_PythonCallback(PyObject* p) : self(p) {}
  double operator()(const double R, const double z) const {
    return call_method<double>(self, "Phi",
                               static_cast<numeric::array>(make_tuple(R, z)));
  }
  double operator()(const double R, const double z, double& dR,
                    double& dz) const {
    auto g = call_method<numeric::array>(self, "Forces",
                                         numeric::array(make_tuple(R, z)));
    dR = -extract<double>(g[0]);
    dz = -extract<double>(g[1]);
    return (*this)(R, z);
  }

  double LfromRc(const double R, double* dR) const {
    double dPR, dPz;
    (*this)(R, 0., dPR, dPz);
    return sqrt(R * R * R * dPR);
  }

  double RfromLc(const double L, double* dR) const {
    bool more = false;
    double R, lR = 0., dlR = 0.001, dPR, dPz, LcR, oldL;
    R = exp(lR);
    (*this)(R, 0., dPR, dPz);
    LcR = pow(R * R * R * dPR, 0.5);
    if (LcR == L) return R;
    if (L > LcR) more = true;
    oldL = LcR;

    for (;;) {
      lR += (more) ? dlR : -dlR;
      R = exp(lR);
      (*this)(R, 0., dPR, dPz);
      LcR = pow(R * R * R * dPR, 0.5);
      if (LcR == L) return R;
      if ((L < LcR && L > oldL) || (L > LcR && L < oldL)) {
        R = (more) ? exp(lR - 0.5 * dlR) : exp(lR + 0.5 * dlR);
        return R;
      }
      oldL = LcR;
    }
  }
  Frequencies KapNuOm(const double R) const {
    Frequencies epi;
    std::cerr<<"Not yet implemented"<<std::endl;
    return epi;
  }

  PyObject* self;
};

// We can't pass ifstream from python so we have this horrible hack

class GalPot_ForPython: public Potential{
private:
  GalaxyPotential *GP;
public:
  GalPot_ForPython(std::string filename){
    std::ifstream file(filename.c_str());
    GP = new GalaxyPotential(file);
    file.close();
  }
  double operator()(const double R, const double z) const {
    return (*GP)(R,z);
  }
  double operator()(const double R, const double z, double& dR,
                    double& dz) const {return (*GP)(R,z,dR,dz);}
  double LfromRc(const double R, double* dR) const{return GP->LfromRc(R,dR);}
  double RfromLc(const double L, double* dR) const{return GP->RfromLc(L,dR);}
  Frequencies KapNuOm(const double R) const {return GP->KapNuOm(R);}
};

// ==========================================================================
// Now define members of python library

BOOST_PYTHON_MODULE_INIT(Torus_py) {
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  class_<Potential, boost::noncopyable,
         boost::shared_ptr<Potential_PythonCallback> >("WDPotential", init<>())
      .def("__call__",
           static_cast<double (Potential::*)(const double, const double) const>(
               &Potential::operator()));
  class_<LogPotential, bases<Potential> >(
      "LogarithmicPotential", init<double, double, double, double>());
  class_<MiyamotoNagaiPotential, bases<Potential> >(
      "MiyamotoNagaiPotential", init<double, double, double>());
  class_<GalPot_ForPython, bases<Potential> >("GalaxyPotential",
    init<std::string>());

  class_<Torus, boost::noncopyable>("Torus", init<>())
      .def("energy",&Torus::energy)
      .def("fsample",&Torus::fsample)
      .def("actions",&Torus::actions).def("action",&Torus::action)
      // .def("omegas",&Torus::omega)
      .def("errors",&Torus::errors).def("error",&Torus::error)
      // .def("toymap_params",&Torus::TP)
      // .def("SN",&Torus::SN)
      .def("AngMap",&Torus::AP)
      .def("minR",&Torus::minR)
      .def("maxR",&Torus::maxR)
      .def("maxz",&Torus::maxz)
      .def("dS1",&Torus::dS1)
      .def("dS2",&Torus::dS2)
      .def("dS3",&Torus::dS3)
      .def("n1",&Torus::n1)
      .def("n2",&Torus::n2)
      .def("make_torus", &Torus::AutoTorus)
      .def("autofit", &Torus::AutoFit)
      .def("map3D", &Torus::Map3D);

  // For action, angle, frequency vectors
  to_python_converter<Vector<double, 3>, WDvector_to_ndarray<double, 3> >();
  WDvector_from_ndarray<double, 3>();

  to_python_converter<PSPT, PSPT_to_ndarray>();
  PSPT_from_ndarray();
  import_array();
}

// ==========================================================================
