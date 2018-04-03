// VisualSPHysics
// Copyright (C) 2018 Orlando Garcia-Feal orlando@uvigo.es

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

#include <Python.h>
#include <iostream>
#include "SimulationParams.h"
#include "DiffuseCalculator.h"

/*
 * This is a simple Python module to run the foam simulation.
 */

extern "C"
{

  static PyObject * diffuseparticles_run(PyObject *self, PyObject *args){
    const char * dataPath, * filePrefix, * outputPath, * outputPreffix, * exclusionZoneFile;
    SimulationParams sp;

    if(!PyArg_ParseTuple(args, "sssssiiipppppdddddddddddddddddddddd",
			 &dataPath,
			 &filePrefix,
			 &outputPath,
			 &outputPreffix,
			 &exclusionZoneFile,
			 &sp.nstart,
			 &sp.nend,
			 &sp.nzeros,
			 &sp.text_files,
			 &sp.ply_files,
			 &sp.vtk_files,
			 &sp.vtk_diffuse_data,
			 &sp.vtk_fluid_data,
			 &sp.h, &sp.mass, &sp.TIMESTEP,
			 &sp.MINX, &sp.MINY, &sp.MINZ, &sp.MAXX, &sp.MAXY, &sp.MAXZ,
			 &sp.MINTA, &sp.MAXTA, &sp.MINWC, &sp.MAXWC,
			 &sp.MINK, &sp.MAXK, &sp.KTA, &sp.KWC,
			 &sp.SPRAY, &sp.BUBBLES, &sp.LIFEFIME, &sp.KB, &sp.KD
			 )){
      return NULL;
    }

    
    sp.dataPath = dataPath;
    sp.filePrefix = filePrefix;
    sp.outputPath = outputPath;
    sp.outputPreffix = outputPreffix;
    sp.exclusionZoneFile = exclusionZoneFile;
      
    DiffuseCalculator dc(sp);
    dc.runSimulation();
    
    return Py_True;
  }

  static PyMethodDef DiffuseParticlesMethods[] = {
    {"run", diffuseparticles_run, METH_VARARGS, "Run Diffuse Particles Simulation"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
  };

  static struct PyModuleDef diffuseparticlesmodule = {
    PyModuleDef_HEAD_INIT,
    "diffuseparticles", /* Name of the module */
    NULL, /* Module of documentation */
    -1, /* size of per-interpreter state of the module,
	   or -1 if the module keeps state in global variables. */
    DiffuseParticlesMethods
  };

  PyMODINIT_FUNC PyInit_diffuseparticles(void) {
    return PyModule_Create(&diffuseparticlesmodule);
  }
  
}
