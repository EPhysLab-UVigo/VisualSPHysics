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

#ifndef BUCKETCONTAINER_H
#define BUCKETCONTAINER_H

#include <vector>
#include <tuple>
#include <cmath>
#include <functional>
#include <array>

/**
   \brief Template class that implements a particle container.
   This template implements a data structure to accelerate the search for neighbour particles.
   In order to find the neighbours of a given point in a point vector, it is necessary to iterate 
   over all that vector. This is too complex to address cases with a large number of particles.
   To improve the performance, the spatial domain is divided into cubes with an edge of size "h". 
   All the particles are inserted into these buckets. In this way, when performing the search of 
   neighbour particles of a given point, only the particles in the bucket in which this point is 
   placed and the 26 surrounding buckets are evaluated.
 */
template <class T>
class BucketContainer {

 private:
  //long nPoints;
  double xmin, xmax, ymin, ymax, zmin, zmax; // Maximum and minimum coordinates in space
  long nx,ny,nz; // Number of buckets on each dimension
  double width; // Size of each bucket

  const int nneig = 27;
  const int addvals[27][3] = {{ 0, 0, 0},
			      { 1, 0, 0},
			      {-1, 0, 0},				
			      { 0, 1, 0},
			      { 1, 1, 0},
			      {-1, 1, 0},				
			      { 0,-1, 0},
			      { 1,-1, 0},
			      {-1,-1, 0},				
			      { 0, 0, 1},
			      { 1, 0, 1},
			      {-1, 0, 1},				
			      { 0, 1, 1},
			      { 1, 1, 1},
			      {-1, 1, 1},				
			      { 0,-1, 1},
			      { 1,-1, 1},
			      {-1,-1, 1},				
			      { 0, 0,-1},
			      { 1, 0,-1},
			      {-1, 0,-1},				
			      { 0, 1,-1},
			      { 1, 1,-1},
			      {-1, 1,-1},				
			      { 0,-1,-1},
			      { 1,-1,-1},
			      {-1,-1,-1}};
  
  std::vector<std::vector<T>> buckets; // Buckets
  std::vector<std::pair<long, std::vector<T> &>> nebuckets; // Not empty buckets

 public:
  /**
     Class constructor. It creates an empty data structure of the given size.
     \param xmin Domain limits: min x value.
     \param xmax Domain limits: max x value.
     \param ymin Domain limits: min y value.
     \param ymax Domain limits: max y value.
     \param zmin Domain limits: min z value.
     \param zmax Domain limits: max z value.
     \param h Cell size.
   */
  BucketContainer(double xmin, double xmax,
		  double ymin, double ymax,
		  double zmin, double zmax, double h);

  /**
     Adds and element given its spatial coordinates.
     \param e Element to be inserted.
     \param x Coordinate x.
     \param y Coordinate y.
     \param z Coordinate z.
     \return Returns true if the operation is executed correctly.
   */
  bool addElement(T e, double  x, double  y, double  z);

  /**
     Given some spatial coordinates, return the index of the bucket to which it belongs.
     \param x Coordinate x.
     \param y Coordinate y.
     \param z Coordinate z.
     \return Bucket index.
   */
  long getBucketNumber(double  x, double  y, double  z) const;

  /**
     Given some spatial coordinates, return the coordinates of the bucket inside the bucket matrix.
     \param x Coordinate x.
     \param y Coordinate y.
     \param z Coordinate z.
     \return Bucket coordinates.
   */
  std::array<long,3> getBucketCoords(double  x, double  y, double  z) const;

  /**
     Given the index of a bucket, returns an array with the coordinates of the bucket inside the bucket matrix.
     \param n Bucket index.
     \return Bucket coordinates.
   */
  std::array<long,3> getBucketCoords(long n) const;

  /**
     Returns the total number of elements in the structure. 
     Note: this operation loops over the whole structure, so it is relatively slow.
     \return Total number of elements.
   */
  long getNElements() const;

  /**
     \return A reference to the bucket vector that contains all the elements. Including those empty buckets.
   */  
  std::vector<std::vector<T>> & getBuckets();

  /**
     \return A reference to a vector of pairs id-bucket. This function only returns those non-empty buckets.
  */
  std::vector<std::pair<long, std::vector<T> &>> & getNoEmptyBuckets();

  /**
     Given a bucket index, return a vector with the elements in that bucket and the 26 surrounding buckets.
     \param nbucket Bucket index.
     \return Element vector.
   */
  std::vector<T> getSurroundingElements(long nbucket);

  /**
     Given a bucket index, returns a vector of pointers to the pointed bucket and the 26 surrounding buckets.
     \param nbucket Bucket index.
     \return Vector of pointers to buckets.
   */
  std::vector<std::vector<T> *> getSurroundingBuckets(long nbucket);

  /**
     Given the coordinates of a bucket, returns a vector of pointers to the pointed bucket and the 26 surrounding buckets.
     \param bp Array with the coordinates of the bucket.
     \return Vector of pointers to buckets.
   */
  std::vector<std::vector<T> *> getSurroundingBuckets(std::array<long,3> bp);

  /**
     Given the coordinates of an element, returns a vector of pointers to the pointed bucket and the 26 surrounding buckets.
     \param px Coordinate x.
     \param py Coordinate y.
     \param pz Coordinate z.
     \return Vector of pointers to buckets.
   */
  std::vector<std::vector<T> *> getSurroundingBuckets(double px, double py, double pz);

  /**
     Given the coordinates of an element, returns a vector of pointers to the pointed bucket and the 26 surrounding buckets.
     \param pos Array with the coordinates of the element.
     \return Vector of pointers to buckets.
   */  
  std::vector<std::vector<T> *> getSurroundingBuckets(std::array<double,3> pos);

};

/* Method definitions */

template <class T>
BucketContainer<T>::BucketContainer(double xmin, double xmax,
				 double ymin, double ymax,
				 double zmin, double zmax, double h):
  xmin(xmin), xmax(xmax),
  ymin(ymin), ymax(ymax),
  zmin(zmin), zmax(zmax), width(h){

  // Make subdivisions
  nx = std::floor((xmax - xmin) / h) + 1;
  ny = std::floor((ymax - ymin) / h) + 1;
  nz = std::floor((zmax - zmin) / h) + 1;

  buckets.resize(nx*ny*nz);
}

template <class T>
bool BucketContainer<T>::addElement(T e, double  x, double  y, double  z){
  if(x <= xmin || x >= xmax ||
     y <= ymin || y >= ymax ||
     z <= zmin || z >= zmax ){
    return false; // Out of bounds
  }
  int nbucket = getBucketNumber(x, y, z);
  buckets.at(nbucket).push_back(e);
  return true;
}

template <class T>
long BucketContainer<T>::getBucketNumber(double  x, double  y, double  z) const{
  auto coords = getBucketCoords(x,y,z);
  return coords[0] + nx * coords[1] + nx * ny * coords[2];
}

template <class T>
std::array<long,3> BucketContainer<T>::getBucketCoords(double  x, double  y, double  z) const {
  return std::array<long,3>{(long)((x - xmin)/width), 
      (long)((y - ymin)/width), 
      (long)((z - zmin)/width)};
}

template <class T>
std::array<long,3> BucketContainer<T>::getBucketCoords(long n) const {
  return std::array<long,3>{ ((n % (nx*ny)) % nx),
    (n % (nx*ny)) / nx,
    n / (nx*ny)
  };
}

template <class T>
long BucketContainer<T>::getNElements() const {
  long ne = 0;
  for(auto &i : buckets)
    ne += i.size();
  return ne;
}

template <class T>
std::vector<std::vector<T>> & BucketContainer<T>::getBuckets(){
  return buckets;
}

template <class T>
std::vector<std::pair<long, std::vector<T> &>> & BucketContainer<T>::getNoEmptyBuckets(){
  if(nebuckets.size() == 0){ 
    //Compact bucket list
    for(long nbucket=0; nbucket < buckets.size(); nbucket++){
      if(buckets[nbucket].size() > 0){
	nebuckets.push_back(std::make_pair(nbucket, std::ref(buckets[nbucket])));
      }
    }
  }  
  return nebuckets;
}

template <class T>
std::vector<T> BucketContainer<T>::getSurroundingElements(long nbucket){
  auto buckets = getSurroundingBuckets(nbucket);
  std::vector<T> ret;
  for(auto bucket : buckets){
    for(auto p : *bucket){
      ret.push_back(p);
    }
  }
  return ret;
}

template <class T>
std::vector<std::vector<T> *> BucketContainer<T>::getSurroundingBuckets(long nbucket){
  return getSurroundingBuckets(getBucketCoords(nbucket));
}

template <class T>
std::vector<std::vector<T> *> BucketContainer<T>::getSurroundingBuckets(std::array<long,3> bp){
  std::vector<std::vector<T> *> retvec;
  
  for(long i=0; i<nneig; i++){
    long vx = bp[0] + addvals[i][0], 
      vy = bp[1] + addvals[i][1], 
      vz = bp[2] + addvals[i][2];
    if(vx>=0 && vx<nx &&
       vy>=0 && vy<ny &&
       vz>=0 && vz<nz){
      retvec.push_back(&buckets[vx + nx * vy + nx * ny * vz]);
    }
  }
  return retvec;
}

template <class T>
std::vector<std::vector<T> *> BucketContainer<T>::getSurroundingBuckets(double px, double py, double pz){
  auto bp = getBucketCoords(px,py,pz);
  return getSurroundingBuckets(bp);  
}

template <class T>
std::vector<std::vector<T> *> BucketContainer<T>::getSurroundingBuckets(std::array<double,3> pos){
  return getSurroundingBuckets(pos[0],pos[1],pos[2]);
}

#endif
