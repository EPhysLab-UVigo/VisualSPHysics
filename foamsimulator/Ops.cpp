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

#include "Ops.h"
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <cstdio>
#include <limits>
#include <cmath>
#include <iomanip>
#include <math.h>

namespace ops {

  void printArray(std::array<double,3> x){
    std::cout<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<"["<<x[0]<<" "<<x[1]<<" "<<x[2]<<"]"<<std::endl;
  }

  void printArray(std::array<int,3> x){
    std::cout<<"["<<x[0]<<" "<<x[1]<<" "<<x[2]<<"]"<<std::endl;
  }

  template<typename T>
  T dfmin(T a, T b)
  {
	  if (a < b)
		  return a;
	  return b;
  }

  template<typename T>
  T dfmax(T a, T b)
  {
	  if (a > b)
		  return a;
	  return b;
  }

  void printVectorStats(std::vector<double> v){
    std::cout << "Mean: " << std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size() << std::endl
	      << "Min: " << std::accumulate(v.begin(), v.end(), std::numeric_limits<double>::max(), dfmin<double>) << std::endl
	      << "Max: " << std::accumulate(v.begin(), v.end(), std::numeric_limits<double>::min(), dfmax<double>) << std::endl
	      << "Min non-zero: " << std::accumulate(v.begin(), v.end(), std::numeric_limits<double>::max(),
						     [](double a, double b){
						       if (a == 0.0)
							 return b;
						       else if (b == 0.0)
							 return a;
						       return fmin(a,b);
						     }) << std::endl;
  }

  // Velocity difference between two particles
  // double * substract(double *  xi, double *  xj){
  //   double * rval = new double(3);
  //   rval[0] = xi[0] - xj[0];
  //   rval[1] = xi[1] - xj[1];
  //   rval[2] = xi[2] - xj[2];
  //   return rval;
  // }
  
  std::array<double,3> substract(std::array<double,3> xi, std::array<double,3> xj){
    std::array<double,3> rval;
    rval[0] = xi[0] - xj[0];
    rval[1] = xi[1] - xj[1];
    rval[2] = xi[2] - xj[2];
    return rval;
  }

  double magnitude(std::array<double,3> x){
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);  
  }

  std::array<double,3> normalize(std::array<double,3> x){
    double mag = magnitude(x);
    return std::array<double,3>{x[0]/mag, x[1]/mag, x[2]/mag};
  }

  double dotProduct(std::array<double,3> xi, std::array<double,3> xj){
    return xi[0]*xj[0] + xi[1]*xj[1] + xi[2]*xj[2];
  }

  std::array<double,3> distanceVector(std::array<double,3> xi, std::array<double,3> xj){
    std::array<double,3> d1 = substract(xi,xj);
    double d2 = magnitude(d1);
    return std::array<double,3>{d1[0]/d2, d1[1]/d2, d1[2]/d2};
  }

}
