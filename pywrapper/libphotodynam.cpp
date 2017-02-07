/*
  Wrappers for Josh Carter's libphotodynam routines

  Need to install for Python.h header:
    apt-get install python-dev

  Author:
    Kyle Conroy,
    Martin Horvat, October 2016
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <iostream>
#include <Python.h>
#include <numpy/arrayobject.h>


// Porting to Python 3
// Ref: http://python3porting.com/cextensions.html
#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
        static struct PyModuleDef moduledef = { \
          PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
        ob = PyModule_Create(&moduledef);

  // adding missing declarations and functions
  #define PyString_Type PyBytes_Type
  #define PyString_AsString PyBytes_AsString
  #define PyString_Check PyBytes_Check

#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
        ob = Py_InitModule3(name, methods, doc);
#endif

#include "n_body_state.h"
#include "n_body.h"

static PyObject *kep2cartesian(PyObject *self, PyObject *args) {

  // printf("entered: kep2cartesian\n");
  // parse input arguments
  PyObject *mass_o, *a_o, *e_o, *inc_o, *om_o, *ln_o, *ma_o;
  double *mass, *a, *e, *inc, *om, *ln, *ma, *mj, t0;
  int Nobjects, i, j;

  // printf("parsing arguments...\n");
  if (!PyArg_ParseTuple(args, "OOOOOOOd", &mass_o, &a_o, &e_o, &inc_o, &om_o, &ln_o, &ma_o, &t0))
      return NULL;
  // printf("arguments parsed!\n");

  // get lengths of lists
  Nobjects = PyTuple_Size(mass_o);

  // convert objects to C arrays
  //~ printf("converting to C arrays...\n");
  mass = new double[Nobjects];
  a = new double[Nobjects-1];
  e = new double[Nobjects-1];
  inc = new double[Nobjects-1];
  om = new double[Nobjects-1];
  ln = new double[Nobjects-1];
  ma = new double[Nobjects-1];
  mj = new double[Nobjects-1];  // NOT REALLY SURE WHAT THIS DOES

  for (i = 0; i < Nobjects; i++){
      mass[i] = PyFloat_AsDouble(PyTuple_GetItem(mass_o, i));
  }
  for (i = 0; i < Nobjects-1; i++){
    a[i] = PyFloat_AsDouble(PyTuple_GetItem(a_o, i));
    e[i] = PyFloat_AsDouble(PyTuple_GetItem(e_o, i));
    inc[i] = PyFloat_AsDouble(PyTuple_GetItem(inc_o, i));
    om[i] = PyFloat_AsDouble(PyTuple_GetItem(om_o, i));
    ln[i] = PyFloat_AsDouble(PyTuple_GetItem(ln_o, i));
    ma[i] = PyFloat_AsDouble(PyTuple_GetItem(ma_o, i));
  }

  // initialize output
  //~ printf("initializing output...\n");
  PyObject *dict = PyDict_New();

  PyObject *x = PyTuple_New(Nobjects);
  PyObject *y = PyTuple_New(Nobjects);
  PyObject *z = PyTuple_New(Nobjects);
  PyObject *vx = PyTuple_New(Nobjects);
  PyObject *vy = PyTuple_New(Nobjects);
  PyObject *vz = PyTuple_New(Nobjects);

  // create Nbody state
  NBodyState state(mass, a, e, inc, om, ln, ma, Nobjects, t0);

  state.kep_elements(mj,a,e,inc,om,ln,ma);

  for (j = 0; j < Nobjects; j++){
    PyTuple_SetItem(x, j, Py_BuildValue("d", state.X_B(j)));
    PyTuple_SetItem(y, j, Py_BuildValue("d", state.Y_B(j)));
    PyTuple_SetItem(z, j, Py_BuildValue("d", state.Z_B(j)));
    PyTuple_SetItem(vx, j, Py_BuildValue("d", state.V_X_B(j)));
    PyTuple_SetItem(vy, j, Py_BuildValue("d", state.V_Y_B(j)));
    PyTuple_SetItem(vz, j, Py_BuildValue("d", state.V_Z_B(j)));
  }

  PyDict_SetItem(dict, Py_BuildValue("s", "x"), x);
  PyDict_SetItem(dict, Py_BuildValue("s", "y"), y);
  PyDict_SetItem(dict, Py_BuildValue("s", "z"), z);
  PyDict_SetItem(dict, Py_BuildValue("s", "vx"), vx);
  PyDict_SetItem(dict, Py_BuildValue("s", "vy"), vy);
  PyDict_SetItem(dict, Py_BuildValue("s", "vz"), vz);

  return dict;
}


static PyObject *do_dynamics(PyObject *self, PyObject *args) {

  //~ printf("entered: do_dynamics\n");
  // parse input arguments

  PyObject *mass_o, *a_o, *e_o, *inc_o, *om_o, *ln_o, *ma_o, *times_o;

  double *mass, *a, *e, *inc, *om, *ln, *ma, *times, *mj;

  int Nobjects, Ntimes, i, j, ltte, return_keplerian;

  double t0, maxh, orbit_error;

  //~ printf("parsing arguments...\n");
  if (!PyArg_ParseTuple(args, "OOOOOOOOdddii",
        &times_o, &mass_o, &a_o, &e_o, &inc_o, &om_o, &ln_o,
        &ma_o, &t0, &maxh, &orbit_error, &ltte, &return_keplerian)
      )
      return NULL;

  // get lengths of lists
  Nobjects = PyTuple_Size(mass_o);
  Ntimes = PyTuple_Size(times_o);

  // convert objects to C arrays
  //~ printf("converting to C arrays...\n");
  mass = new double[Nobjects];
  a = new double[Nobjects-1];
  e = new double[Nobjects-1];
  inc = new double[Nobjects-1];
  om = new double[Nobjects-1];
  ln = new double[Nobjects-1];
  ma = new double[Nobjects-1];
  mj = new double[Nobjects-1];  // NOT REALLY SURE WHAT THIS DOES
  times = new double[Ntimes];

  for (i = 0; i < Nobjects; i++)
    mass[i] = PyFloat_AsDouble(PyTuple_GetItem(mass_o, i));

  for (i = 0; i < Nobjects-1; i++){
    a[i] = PyFloat_AsDouble(PyTuple_GetItem(a_o, i));
    e[i] = PyFloat_AsDouble(PyTuple_GetItem(e_o, i));
    inc[i] = PyFloat_AsDouble(PyTuple_GetItem(inc_o, i));
    om[i] = PyFloat_AsDouble(PyTuple_GetItem(om_o, i));
    ln[i] = PyFloat_AsDouble(PyTuple_GetItem(ln_o, i));
    ma[i] = PyFloat_AsDouble(PyTuple_GetItem(ma_o, i));
  }

  for (i = 0; i < Ntimes; i++)
    times[i] = PyFloat_AsDouble(PyTuple_GetItem(times_o, i));


  // initialize output
  //~ printf("initializing output...\n");
  PyObject *dict = PyDict_New();

  PyObject *x = PyTuple_New(Ntimes);
  PyObject *y = PyTuple_New(Ntimes);
  PyObject *z = PyTuple_New(Ntimes);
  PyObject *vx = PyTuple_New(Ntimes);
  PyObject *vy = PyTuple_New(Ntimes);
  PyObject *vz = PyTuple_New(Ntimes);

  // if (return_keplerian > 0) {
  PyObject *kepl_a = PyTuple_New(Ntimes);
  PyObject *kepl_e = PyTuple_New(Ntimes);
  PyObject *kepl_in = PyTuple_New(Ntimes);
  PyObject *kepl_o = PyTuple_New(Ntimes);
  PyObject *kepl_ln = PyTuple_New(Ntimes);
  PyObject *kepl_m = PyTuple_New(Ntimes);
  // }

  // create Nbody state
  NBodyState state(mass, a, e, inc, om, ln, ma, Nobjects, t0);

  //state.kep_elements(mj,a,e,inc,om,ln,ma);


  for (i = 0; i < Ntimes; i++) {
    // integrate to this time

    state.kep_elements(mj,a,e,inc,om,ln,ma);

    state(times[i],maxh,orbit_error,1e-16);
    //~ state.kep_elements(mj,a,e,inc,om,ln,ma);

    // need to (re)create tuples to hold positions for each object at this time
    PyObject *xi = PyTuple_New(Nobjects);
    PyObject *yi = PyTuple_New(Nobjects);
    PyObject *zi = PyTuple_New(Nobjects);
    PyObject *vxi = PyTuple_New(Nobjects);
    PyObject *vyi = PyTuple_New(Nobjects);
    PyObject *vzi = PyTuple_New(Nobjects);

    PyObject *kepl_ai = PyTuple_New(Nobjects);
    PyObject *kepl_ei = PyTuple_New(Nobjects);
    PyObject *kepl_ini = PyTuple_New(Nobjects);
    PyObject *kepl_oi = PyTuple_New(Nobjects);
    PyObject *kepl_lni = PyTuple_New(Nobjects);
    PyObject *kepl_mi = PyTuple_New(Nobjects);

    for (j = 0; j < Nobjects; j++){
      if (ltte > 0) {
        PyTuple_SetItem(xi, j, Py_BuildValue("d", state.X_LT(j)));
        PyTuple_SetItem(yi, j, Py_BuildValue("d", state.Y_LT(j)));
        PyTuple_SetItem(zi, j, Py_BuildValue("d", state.Z_LT(j)));
        PyTuple_SetItem(vxi, j, Py_BuildValue("d", state.V_X_LT(j)));
        PyTuple_SetItem(vyi, j, Py_BuildValue("d", state.V_Y_LT(j)));
        PyTuple_SetItem(vzi, j, Py_BuildValue("d", state.V_Z_LT(j)));
      } else {
        PyTuple_SetItem(xi, j, Py_BuildValue("d", state.X_B(j)));
        PyTuple_SetItem(yi, j, Py_BuildValue("d", state.Y_B(j)));
        PyTuple_SetItem(zi, j, Py_BuildValue("d", state.Z_B(j)));
        PyTuple_SetItem(vxi, j, Py_BuildValue("d", state.V_X_B(j)));
        PyTuple_SetItem(vyi, j, Py_BuildValue("d", state.V_Y_B(j)));
        PyTuple_SetItem(vzi, j, Py_BuildValue("d", state.V_Z_B(j)));
      }

      // if (return_keplerian > 0) {
          // state.kep_elements(mji, ai, ei, ini, oi, lni, mi);
      PyTuple_SetItem(kepl_ai, j, Py_BuildValue("d", a));
      PyTuple_SetItem(kepl_ei, j, Py_BuildValue("d", e));
      PyTuple_SetItem(kepl_ini, j, Py_BuildValue("d", inc));
      PyTuple_SetItem(kepl_oi, j, Py_BuildValue("d", om));
      PyTuple_SetItem(kepl_lni, j, Py_BuildValue("d", ln));
      PyTuple_SetItem(kepl_mi, j, Py_BuildValue("d", ma));
      // }
    }

    PyTuple_SetItem(x, i, xi);
    PyTuple_SetItem(y, i, yi);
    PyTuple_SetItem(z, i, zi);
    PyTuple_SetItem(vx, i, vxi);
    PyTuple_SetItem(vy, i, vyi);
    PyTuple_SetItem(vz, i, vzi);

    // if (return_keplerian > 0) {
    PyTuple_SetItem(kepl_a, i, kepl_ai);
    PyTuple_SetItem(kepl_e, i, kepl_ei);
    PyTuple_SetItem(kepl_in, i, kepl_ini);
    PyTuple_SetItem(kepl_o, i, kepl_oi);
    PyTuple_SetItem(kepl_ln, i, kepl_lni);
    PyTuple_SetItem(kepl_m, i, kepl_mi);
    // }
  }

  // prepare dictionary for returning output
  PyDict_SetItem(dict, Py_BuildValue("s", "t"), times_o);
  PyDict_SetItem(dict, Py_BuildValue("s", "x"), x);
  PyDict_SetItem(dict, Py_BuildValue("s", "y"), y);
  PyDict_SetItem(dict, Py_BuildValue("s", "z"), z);
  PyDict_SetItem(dict, Py_BuildValue("s", "vx"), vx);
  PyDict_SetItem(dict, Py_BuildValue("s", "vy"), vy);
  PyDict_SetItem(dict, Py_BuildValue("s", "vz"), vz);

  // if (return_keplerian > 0) {
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_a"), kepl_a);
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_e"), kepl_e);
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_in"), kepl_in);
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_o"), kepl_o);
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_ln"), kepl_ln);
  PyDict_SetItem(dict, Py_BuildValue("s", "kepl_m"), kepl_m);
  // }

  // return format (for a triple) will be something like
  // {'t': [...], 'x': [(t1, t2, t3), (t1, t2, t3), ...], 'y': ..., 'z': ...}

  return dict;

}

int PyDict_SetItemStringMatrix(PyObject *p, const char *key, const int &n, const int & m, double *a){
 
  npy_intp dims[2] = {n, m};
  
  int size = n*m;

  double *data = new double [size];
  
  for (int i = 0; i < size; ++i) data[i] = a[i];
    
  PyObject *val = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data); 

  PyArray_ENABLEFLAGS((PyArrayObject *)val, NPY_ARRAY_OWNDATA);
   
  int status = PyDict_SetItemString(p, key, val);
 
  Py_XDECREF(val);
 
  return status;
}

static PyObject *do_dynamics_numpy(PyObject *self, PyObject *args) {

  PyArrayObject *mass_o, *a_o, *e_o, *inc_o, *om_o, *ln_o, *ma_o, *times_o;

  int ltte, return_keplerian;

  double t0, maxh, orbit_error;

  //~ printf("parsing arguments...\n");
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!dddii",
        &PyArray_Type, &times_o, 
        &PyArray_Type, &mass_o, 
        &PyArray_Type, &a_o, 
        &PyArray_Type, &e_o, 
        &PyArray_Type, &inc_o,
        &PyArray_Type, &om_o,
        &PyArray_Type, &ln_o,
        &PyArray_Type, &ma_o, 
        &t0, 
        &maxh, 
        &orbit_error, 
        &ltte, 
        &return_keplerian)
      )
      return NULL;

  // get lengths of lists
  int 
    Nobjects = PyArray_DIM(mass_o, 0),
    Nobjects_ = Nobjects - 1,
    Ntimes = PyArray_DIM(times_o, 0);

  double
    *times = (double*)PyArray_DATA(times_o),
    *mass = (double*)PyArray_DATA(mass_o),
    *a =    (double*)PyArray_DATA(a_o),
    *e =    (double*)PyArray_DATA(e_o),
    *inc =  (double*)PyArray_DATA(inc_o),
    *om =   (double*)PyArray_DATA(om_o),
    *ln =   (double*)PyArray_DATA(ln_o),
    *ma =   (double*)PyArray_DATA(ma_o);

  int size = Ntimes*Nobjects;
  
  double 
    *buf = new double [12*size],
    *x = buf,
    *y = buf + size,
    *z = buf + 2*size,
    *vx = buf + 3*size,
    *vy = buf + 4*size,
    *vz = buf + 5*size,
    *kepl_a = buf + 6*size,
    *kepl_e = buf + 7*size,
    *kepl_inc = buf + 8*size,
    *kepl_om = buf + 9*size,
    *kepl_ln = buf + 10*size,
    *kepl_ma = buf + 11*size,

    *tmp = new double [Nobjects];
  
  // create Nbody state
  NBodyState state(mass, a, e, inc, om, ln, ma, Nobjects, t0);
  
  for (int i = 0; i < Ntimes; ++i) {
    
    // integrate to this time
    state.kep_elements(tmp, kepl_a, kepl_e, kepl_inc, kepl_om, kepl_ln, kepl_ma);
    
    state(times[i], maxh, orbit_error, 1e-16);
    
    for (int j = 0; j < Nobjects; ++j)
      if (ltte > 0) {
        *(x++) = state.X_LT(j);
        *(y++) = state.Y_LT(j);
        *(z++) = state.Z_LT(j);
        *(vx++) = state.V_X_LT(j);
        *(vy++) = state.V_Y_LT(j);
        *(vz++) = state.V_Z_LT(j);
      } else {
        *(x++) = state.X_B(j);
        *(y++) = state.Y_B(j);
        *(z++) = state.Z_B(j);
        *(vx++) = state.V_X_B(j);
        *(vy++) = state.V_Y_B(j);
        *(vz++) = state.V_Z_B(j);
     }

    kepl_a += Nobjects_;
    kepl_e += Nobjects_;
    kepl_inc += Nobjects_;
    kepl_om += Nobjects_;
    kepl_ln += Nobjects_;
    kepl_ma += Nobjects_;
  }

  // prepare dictionary for returning output
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "t", (PyObject*)times_o);
  PyDict_SetItemStringMatrix(dict, "x",  Ntimes, Nobjects, buf);
  PyDict_SetItemStringMatrix(dict, "y",  Ntimes, Nobjects, buf + size);
  PyDict_SetItemStringMatrix(dict, "z",  Ntimes, Nobjects, buf + 2*size);
  PyDict_SetItemStringMatrix(dict, "vx", Ntimes, Nobjects, buf + 3*size);
  PyDict_SetItemStringMatrix(dict, "vy", Ntimes, Nobjects, buf + 4*size);
  PyDict_SetItemStringMatrix(dict, "vz", Ntimes, Nobjects, buf + 5*size);
  
  PyDict_SetItemStringMatrix(dict, "kepl_a",  Ntimes, Nobjects_, buf + 6*size);
  PyDict_SetItemStringMatrix(dict, "kepl_e",  Ntimes, Nobjects_, buf + 7*size);
  PyDict_SetItemStringMatrix(dict, "kepl_in", Ntimes, Nobjects_, buf + 8*size);
  PyDict_SetItemStringMatrix(dict, "kepl_o",  Ntimes, Nobjects_, buf + 9*size);
  PyDict_SetItemStringMatrix(dict, "kepl_ln", Ntimes, Nobjects_, buf + 10*size);
  PyDict_SetItemStringMatrix(dict, "kepl_m",  Ntimes, Nobjects_, buf + 11*size);
   
  delete [] tmp;
  delete [] buf;
  
  return dict;
}


static PyMethodDef Methods[] = {
    {"do_dynamics",      do_dynamics,   METH_VARARGS, "Do N-body dynamics to get positions as a function of time"},
    {"do_dynamics_numpy",      do_dynamics_numpy,   METH_VARARGS, "Do N-body dynamics to get positions as a function of time"},
    {"kep2cartesian",    kep2cartesian, METH_VARARGS, "Convert from keplerian to cartesian initial conditions"},
    {NULL,               NULL,             0,            NULL}
};


static char const *Docstring =
  "Module wraping John Carters photodynam routines.";


/* module initialization */
MOD_INIT(photodynam) {

  PyObject *backend;

  MOD_DEF(backend, "photodynam", Docstring, Methods)

  if (!backend) return MOD_ERROR_VAL;

  // Added to handle Numpy arrays
  // Ref:
  // * http://docs.scipy.org/doc/numpy-1.10.1/user/c-info.how-to-extend.html
  import_array();

  return MOD_SUCCESS_VAL(backend);
}
