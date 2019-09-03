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

#ifndef PLYWRITER_H
#define PLYWRITER_H

#include <string>
#include <vector>
#include <array>

#include "FileWriter.h"
#include "rply/rply.h"
#include "BucketContainer.h"

/**
   \brief Implement writing to ply format files.
   This class is responsible for generating a geometry associated with the 
   foam particles and saving said geometry in PLY format files. To generate 
   the geometry polyhedra with 6 triangles are used. In addition, a clustering 
   phase is performed, in which very close particles are grouped, generating
   larger particles. In this way the randomness is increased and the number 
   of polygons used is reduced.
 */
class PlyWriter : public FileWriter {

private:

  std::string fileName;
  std::string comment;

  std::vector<std::array<double,3>> data;
  std::vector<double> size;

  BucketContainer<std::array<double,3>> bc;

  double h;

  //Reference low-poly
  std::array<double, 15> loVertices; 

  // Vertices of each face
  std::array<int, 18> loFaces; 

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
 PlyWriter(std::string const& name,
	  double xmin, double xmax,
	  double ymin, double ymax,
	  double zmin, double zmax, double h);

 /**
    Adds a comment to the PLY file.
    \param c Text string with the comment.
  */
 void setComment(std::string const& c);

 /**
    Sets the data to write.
    \param d Pointer to vector of arrays of doubles of three components.
  */
 void setData(std::vector<std::array<double,3>> *d);

 /**
    Dumps the data to the file.
    \return Zero if it fails. Any other value otherwise.
  */
 virtual int write();

};

#endif
