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

#ifndef FLUIDDATA_H
#define FLUIDDATA_H

#include <string>
#include <functional>
#include <vector>
#include <tuple>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "BucketContainer.h"

/**
   This structure stores all the data of a fluid particle.
 */
struct particle {
  long id;          ///< Particle id
  std::array<double,3> pos; ///< Postion
  std::array<double,3> vel; ///< Velocity
  double rhop;              ///< Density
};

/**
   \brief This class loads a fluid data vtk file.
   This class loads a vtk file with the information of the fluid particles from DualSPHysics,
   storing it in a BucketContainer class container.
   \see BucketContainer
 */
class FluidData{
 private:
  int exclude;
  std::string exFile;

  BucketContainer<particle> bc;

 public:
    /**
     Class constructor. Creates an empty data structure of the given size.
     \param xmin Domain limits: min x value.
     \param xmax Domain limits: max x value.
     \param ymin Domain limits: min y value.
     \param ymax Domain limits: max y value.
     \param zmin Domain limits: min z value.
     \param zmax Domain limits: max z value.
     \param h Cell size.
   */
  FluidData(double xmin, double xmax,
	    double ymin, double ymax,
	    double zmin, double zmax, double h);

  /**
     Load a vtk file with data of fluid particles.
     \param fileName File name.
     \return True if the file was correctly loaded.
   */
  bool loadFile(std::string const& fileName);

  /**
     Loads a file with thte geometry of an exclusion zone.
     Warning: very slow.
     \param fileName File name.
   */
  void setExclusionZone(std::string const& fileName);

  /**
     Return the particle container.
     \return Pointer to an object of type BucketContainer<particle>
   */
  BucketContainer<particle> * getBucketContainer();  
};

#endif
