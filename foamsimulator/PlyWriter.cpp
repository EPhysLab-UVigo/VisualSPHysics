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

#include <iostream>
#include <cmath>

#include "FileWriter.h"
#include "PlyWriter.h"
#include "Ops.h"


PlyWriter::PlyWriter(std::string const& name,
		     double xmin, double xmax,
		     double ymin, double ymax,
		     double zmin, double zmax, double h):
  bc(xmin, xmax, ymin, ymax, zmin, zmax, h), h(h), fileName(name){

  loVertices = {-1.299038, -0.750000, 0.009955,
		0.000000, 0.000000, -1.492738,
		1.299038, -0.750000, 0.009955,
		0.000000, 1.500000, 0.009955,
		0.000000, 0.000000, 1.512648};

  loFaces = {0, 1, 2,
	     3, 1, 0,
	     2, 1, 3,
	     3, 4, 2,
	     2, 4, 0,
	     0, 4, 3};   
}

void PlyWriter::setComment(std::string const& c){
  comment = c;
}
void PlyWriter::setData(std::vector<std::array<double,3>> *d){
  std::cerr<<" PlyWriter: "<<d->size()<<" input particles"<<std::endl;

  for(unsigned long i=0; i < d->size(); i++){    
    auto p = d->at(i);
    if(p[0] != 0 && p[1] != 0 && p[2] != 0){
      bc.addElement(p, p[0], p[1], p[2]);
    }
  }

  std::cerr<<" Particles non zero: "<< bc.getNElements() << std::endl;
  
  for(auto &bucket : bc.getBuckets()){
    for(unsigned long i=0; i<bucket.size(); i++){
      std::array<double,3> newp = bucket[i];
      int psize = 1;
      for(unsigned long j=0; j<bucket.size(); j++){
	if(j!=i && (ops::magnitude(ops::substract(bucket[i],bucket[j])) < (h/5))) {
	  newp[0] = (newp[0]+bucket[j][0])/2.;
	  newp[1] = (newp[1]+bucket[j][1])/2.;
	  newp[2] = (newp[2]+bucket[j][2])/2.;
	  bucket.erase(bucket.begin()+j);
	  j--;
	  psize++;
	}	
      }
      bucket.erase(bucket.begin()+i);
      i--;
      data.push_back(newp);
      size.push_back(psize);      
    }
  }
  std::cerr<<" PlyWriter: "<<data.size()<<" output spheres"<<std::endl;  
}

int PlyWriter::write(){
  
  p_ply oply = ply_create( fileName.c_str(), PLY_LITTLE_ENDIAN, NULL, 0, NULL );
  if ( !oply ) {
    std::cerr<<"ERROR: Could not create »new.ply«\n";
    return 0;
  }

  /* Add object information to PLY. */
  if ( !ply_add_obj_info( oply, "Diffuse particle geometry." ) ) {
    std::cerr<<"ERROR: Could not add object info.\n";
    return 0;
  }

  /* Add a comment, too. */
  if ( !ply_add_comment( oply, comment.c_str() ) ) {
    std::cerr<<"ERROR: Could not add comment.\n";
    return 0;
  }

  if(data.size() > 0){

    /* Add vertex element. */
    if ( !ply_add_element( oply, "vertex", data.size() * 5 ) ) {
      std::cerr<<"ERROR: Could not add element.\n";
      return 0;
    }

    /* Add vertex properties: x, y, z */
    if ( !ply_add_property( oply, "x", PLY_FLOAT, PLY_UINT, PLY_UINT ) ) {
      std::cerr<<"ERROR: Could not add property x.\n";
      return 0;
    }

    if ( !ply_add_property( oply, "y", PLY_FLOAT, PLY_UINT, PLY_UINT ) ) {
      std::cerr<<"ERROR: Could not add property x.\n";
      return 0;
    }

    if ( !ply_add_property( oply, "z", PLY_FLOAT, PLY_UINT, PLY_UINT ) ) {
      std::cerr<<"ERROR: Could not add property x.\n";
      return 0;
    }

    /* Add face element. */
    if ( !ply_add_element( oply, "face", data.size() * 6 ) ) {
      std::cerr<<"ERROR: Could not add element.\n";
      return 0;
    }

    if ( !ply_add_list_property( oply, "vertex_indices", PLY_UCHAR, PLY_UINT ) ) {
      std::cerr<<"ERROR: Could not add property vertex indices.\n";
      return 0;
    }

    /* Write header to file */
    if ( !ply_write_header( oply ) ) {
      std::cerr<<"ERROR: Could not write header.\n";
      return 0;
    }

    /* 1.- Write the vertices */

    for(unsigned long i=0; i<data.size(); i++){
      for(int iv=0; iv<5; iv++){
	double f = h/(10/std::pow(size[i],1/3));
	ply_write( oply, loVertices[iv*3] *  f + data[i][0]);
	ply_write( oply, loVertices[iv*3+1] * f + data[i][1]);
	ply_write( oply, loVertices[iv*3+2] * f + data[i][2]);
      }
    }

    /* 2.- Write the faces */

    for(unsigned long i=0; i<data.size(); i++){
      for(int iv=0; iv<6; iv++){
	ply_write( oply, 3 );
	ply_write( oply, loFaces[iv*3] + i*5 );
	ply_write( oply, loFaces[iv*3+1] + i*5 );
	ply_write( oply, loFaces[iv*3+2] + i*5 );
      }
    }
  }

  /* Close the file */

  if ( !ply_close( oply ) ) {
    std::cerr<<"ERROR: Could not close file.\n";
    return 0;
  }
	
  return 1;
}
