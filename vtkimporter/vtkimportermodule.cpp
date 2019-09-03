// VisualSPHysics
// Copyright (C) 2019 Orlando Garcia-Feal orlando@uvigo.es

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <Python.h>
#include <iostream>
#include <array>


/**
   \file vtkimportermodule.cpp
   \brief C++ module for Python to load vtk files.
   This Python module allows loading VTK geometry files into Blender. 
   Due to its implementation, it is faster than importing STL or PLY files natively into Blender.
   Moreover, VTK is the format used by the DualSPHysics tools so the workflow is faster.
 */
extern "C"
{

  std::array<double, 15> loVertices = {-1.299038, -0.750000, 0.009955,
				                               0.000000, 0.000000, -1.492738,
				                               1.299038, -0.750000, 0.009955,
				                               0.000000, 1.500000, 0.009955,
				                               0.000000, 0.000000, 1.512648};
  
  std::array<int, 18> loFaces = {0, 1, 2,
				                         3, 1, 0,
				                         2, 1, 3,
				                         3, 4, 2,
				                         2, 4, 0,
				                         0, 4, 3};
  
  /** 
     Load a file corresponding to diffuse particles data. Diffuse particles are fixed points in the space,
     it is necessary to create an associated geometry to each particle, in this case, a 6-triangle polyedron.
     The size of the object will be pointed in the file.
     \param self Pointer to the associated Python object.
     \param args Pointer to the Python object that contains the parameters of the function, in this case, the
     string with the file name.
     \return Pointer to a Python object that contains a tuple with the list of vertices, list of polygons and velocity vectors.
   */
  static PyObject * vtkimporter_loaddiffuse(PyObject *self, PyObject *args){
    const char * FILE_NAME;

    double step = 0.1;

    if(!PyArg_ParseTuple(args, "s", &FILE_NAME))
      return NULL;

    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(FILE_NAME);
    reader->Update();
    
    vtkPolyData * output = reader->GetOutput();
    
    vtkDataArray * points = output->GetPoints()->GetData();
    vtkDataArray * pvel = output->GetPointData()->GetArray("Velocity");

    vtkPointData * pointData = output->GetPointData();
    vtkDataArray * size = pointData->GetArray("Size");

    PyObject *vertices = PyList_New(output->GetPoints()->GetNumberOfPoints()*5),
      *velocity = PyList_New(output->GetPoints()->GetNumberOfPoints()*5),
      *faces = PyList_New(output->GetPoints()->GetNumberOfPoints()*6);
    
    for(long i=0; i < output->GetPoints()->GetNumberOfPoints(); i++){
      double *p = points->GetTuple(i),
	      *v = pvel->GetTuple(i),
	      s = size->GetTuple(i)[0];

      /* 1.- Write the vertices */
      for(int iv=0; iv<5; iv++){
	      double px = loVertices[iv*3] *  s + p[0],
	      py = loVertices[iv * 3 + 1] * s + p[1],
	      pz = loVertices[iv * 3 + 2] * s + p[2];

	      PyList_SET_ITEM(vertices, i*5+iv, Py_BuildValue("(ddd)", px, py, pz));
	      PyList_SET_ITEM(velocity, i*5+iv, Py_BuildValue("(ddd)", px + step * v[0],
				          			py + step * v[1],
				          			pz + step * v[2]));
      }
      
      /* 2.- Write the faces */
      for(int iv=0; iv<6; iv++){
	      PyList_SET_ITEM(faces, i * 6 + iv, Py_BuildValue("(iii)",
                        loFaces[iv * 3 ] + i * 5,
                        loFaces[iv * 3 + 1] + i * 5,
                        loFaces[iv * 3 + 2] + i * 5));
      }
    }

    PyObject *ret = PyTuple_New(3);
    
    PyTuple_SET_ITEM(ret, 0, vertices);
    PyTuple_SET_ITEM(ret, 1, faces);
    PyTuple_SET_ITEM(ret, 2, velocity);
    
    return ret;
  }

  /**
     Load a VTK file with a geometry. Normally this function is employed to load water or floating bodies mesh.
     \param self Pointer to the associated Python object.
     \param args Pointer to the Python object that contains the parameters of the function, in this case, the string with the file name.
     \return Pointer to a Python object that contains a tuple with the list of vertices, list of polygons and velocity vectors.
   */
  static PyObject * vtkimporter_loadvel(PyObject *self, PyObject *args){
    std::cerr<<"DEBUG: VtkImporter v1.5\n";
    const char * FILE_NAME;

    double step = 0.1;

    if(!PyArg_ParseTuple(args, "s", &FILE_NAME))
      return NULL;
    
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(FILE_NAME);
    reader->Update();

    vtkPolyData * output = reader->GetOutput();
    vtkDataArray * points = output->GetPoints()->GetData();
    vtkDataArray * pvel = output->GetPointData()->GetArray("Vel");

    PyObject *vertices = PyList_New(output->GetPoints()->GetNumberOfPoints()),
      *velocity = PyList_New(output->GetPoints()->GetNumberOfPoints());

    std::cerr<<"DEBUG: data read ok";

    for(long i=0; i < output->GetPoints()->GetNumberOfPoints(); i++){
      double *p = points->GetTuple(i),
	      *v = pvel->GetTuple(i);
      PyList_SET_ITEM(vertices, i, Py_BuildValue("(ddd)", p[0], p[1], p[2]));

      PyList_SET_ITEM(velocity, i, Py_BuildValue("(ddd)",
                                                 p[0] + step * v[0],
       						                               p[1] + step * v[1],
       						                               p[2] + step * v[2]));
    }

    std::cerr<<"DEBUG: vertex read ok";

    vtkCellArray * polys = output->GetPolys();
    vtkIdType *indices; // This should be vtkSmartPointer or just regular pointer?
    vtkIdType numberOfPoints;
    long polyCount = 0;

    PyObject *faces = PyList_New(polys->GetNumberOfCells());

    for (polys->InitTraversal(); polys->GetNextCell(numberOfPoints, indices); polyCount++) {
      if(numberOfPoints == 3){
	      PyList_SET_ITEM(faces, polyCount, Py_BuildValue("(iii)",indices[0],indices[1],indices[2]));
      }else{ /* Should be 4 ;-) */
	      PyList_SET_ITEM(faces, polyCount, Py_BuildValue("(iiii)",indices[0],indices[1],indices[2],indices[3]));
      }
    }

    std::cerr<<"DEBUG: faces read ok";

    PyObject *ret = PyTuple_New(3);
    
    PyTuple_SET_ITEM(ret, 0, vertices);
    PyTuple_SET_ITEM(ret, 1, faces);
    PyTuple_SET_ITEM(ret, 2, velocity);
    
    return ret;
  }

  /**
     Load a rope-type object from a vtk file.
     \param self Pointer to the associated Python object.
     \param args Pointer to the Python object that contains the parameters of the function, in this case, the string with the file name.
     \return Pointer to a Python object that contains a tuple with the list of vertices and list of polygons.
   */
  static PyObject * vtkimporter_loadrope(PyObject *self, PyObject *args){
    const char * FILE_NAME;

    if(!PyArg_ParseTuple(args, "s", &FILE_NAME))
      return NULL;
    
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(FILE_NAME);
    reader->Update();

    vtkPolyData * output = reader->GetOutput();
    vtkDataArray * points = output->GetPoints()->GetData();

    PyObject *vertices = PyList_New(output->GetPoints()->GetNumberOfPoints());

    for(long i=0; i < output->GetPoints()->GetNumberOfPoints(); i++){
      double *p = points->GetTuple(i);
      PyList_SET_ITEM(vertices, i, Py_BuildValue("(ddd)",p[0],p[1],p[2]));
    }

    vtkCellArray * lines = output->GetLines();
    int nlines = lines->GetNumberOfCells(),
      nl=0;
    vtkSmartPointer<vtkIdList> clids = vtkSmartPointer<vtkIdList>::New();

    PyObject *vlines = PyList_New(nlines);
    
    lines->InitTraversal();

    while(lines->GetNextCell(clids) != 0){
      int lnids = clids->GetNumberOfIds();
      PyObject * tline = PyList_New(lnids);

      for(long i=0; i < lnids; i++){
	      PyList_SET_ITEM(tline, i, Py_BuildValue("i",clids->GetId(i)));
      }
      
      PyList_SET_ITEM(vlines, nl, tline);
      nl++;
    }    

    PyObject *ret = PyTuple_New(2);
    
    PyTuple_SET_ITEM(ret, 0, vertices);
    PyTuple_SET_ITEM(ret, 1, vlines);
    
    return ret;
  }

  /**
     Load a VTK file with a geometry. Normally this function is employed to load water or floating bodies mesh.
     \param self Pointer to the associated Python object.
     \param args Pointer to the Python object that contains the parameters of the function, in this case, the string with the file name.
     \return Pointer to a Python object that contains a tuple with the list of vertices and list of polygons.
   */
  static PyObject * vtkimporter_load(PyObject *self, PyObject *args){
    std::cerr<<"DEBUG: VtkImporter v1.5\n";
    
    const char * FILE_NAME;

    if(!PyArg_ParseTuple(args, "s", &FILE_NAME))
      return NULL;
    
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(FILE_NAME);
    reader->Update();

    vtkPolyData * output = reader->GetOutput();
    vtkDataArray * points = output->GetPoints()->GetData();

    PyObject *vertices = PyList_New(output->GetPoints()->GetNumberOfPoints());

    for(long i=0; i < output->GetPoints()->GetNumberOfPoints(); i++){
      double *p = points->GetTuple(i);
      PyList_SET_ITEM(vertices, i, Py_BuildValue("(ddd)",p[0],p[1],p[2]));
    }

    vtkCellArray * polys = output->GetPolys();
    vtkIdType *indices;
    vtkIdType numberOfPoints;
    long polyCount = 0;

    PyObject *faces = PyList_New(polys->GetNumberOfCells());
    
    for (polys->InitTraversal(); polys->GetNextCell(numberOfPoints, indices); polyCount++) {
      if(numberOfPoints == 3){
	      PyList_SET_ITEM(faces, polyCount, Py_BuildValue("(iii)",indices[0],indices[1],indices[2]));
      }else{ /* Should be 4 ;-) */
	      PyList_SET_ITEM(faces, polyCount, Py_BuildValue("(iiii)",indices[0],indices[1],indices[2],indices[3]));
      }
    }

    PyObject *ret = PyTuple_New(2);
    
    PyTuple_SET_ITEM(ret, 0, vertices);
    PyTuple_SET_ITEM(ret, 1, faces);
    
    return ret;
  }


  static PyMethodDef VtkImporterMethods[] = {
    {"load", vtkimporter_load, METH_VARARGS, "Load a vtk file."},
    {"loadvel", vtkimporter_loadvel, METH_VARARGS, "Load a vtk file with velocity vectors."},
    {"loadrope", vtkimporter_loadrope, METH_VARARGS, "Load a vtk file with rope data."},
    {"loaddiffuse", vtkimporter_loaddiffuse, METH_VARARGS, "Load a vtk file with diffuse particles data."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
  };

  static struct PyModuleDef vtkimportermodule = {
    PyModuleDef_HEAD_INIT,
    "vtkimpoter", // Name of the module
    NULL,         // Module of documentation 
    -1,           // size of per-interpreter state of the module, or -1 if the module keeps state in global variables.
    VtkImporterMethods    
  };

  /**
     This function initilises the Python module.
     \return Pointer to a Python object.
   */
  PyMODINIT_FUNC PyInit_vtkimporter(void) {
    return PyModule_Create(&vtkimportermodule);
  }
  
}

