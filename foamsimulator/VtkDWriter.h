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

#ifndef VTKWRITER_H
#define VTKWRITER_H

#include "FileWriter.h"
#include "BucketContainer.h"

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>

struct oparticle {
  std::array<double,3> pos; ///< Position vector of a particle.
  std::array<double,3> vel; ///< Velocity vector of a particle.
};

/**
   \brief Implement writing to vtk format files.
   This class is responsible for storing particles in vtk files. A clustering phase 
   is performed in which very close particles are grouped, creating larger particles. 
   In this way the randomness is increased and the number of polygons used is reduced. 
   In the vtk files, only the position of the particles and their size are stored.
*/
class VtkDWriter : public FileWriter {

 private:
  std::string fileName;

  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkDoubleArray> velocity;
  vtkSmartPointer<vtkDoubleArray> size;
  
  //TODO: check memory consumption...
  BucketContainer<oparticle> bc;

  double h;

 public:
  /**
    Class constructor.
    \param name File name to write.
    \param xmin Domain limits: min x value.
    \param xmax Domain limits: max x value.
    \param ymin Domain limits: min y value.
    \param ymax Domain limits: max y value.
    \param zmin Domain limits: min z value.
    \param zmax Domain limits: max z value.
    \param h Cell size.
 */
  VtkDWriter(std::string const& name,
	    double xmin, double xmax,
	    double ymin, double ymax,
	    double zmin, double zmax, double h);

  /**
    Set the data to write including velocity vectors.
    \param d Position vectors.
    \param v Velocity vectors.
  */
  void setData(std::vector<std::array<double,3>> *d, std::vector<std::array<double,3>> *v);
  
  /**
    Dumps the data to the file.
    \return Zero if it fails. Any other value otherwise.
  */
  virtual int write();
  
};


#endif
