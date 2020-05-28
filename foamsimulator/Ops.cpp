// VisualSPHysics
// Copyright (C) 2020 Orlando Garcia-Feal orlando@uvigo.es

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
#include <numeric>
#include <functional>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace ops {

  std::string vectorStats(std::vector<double> &vec){
    auto v(vec);
    auto q1 = v.size() / 4,
         q2 = v.size() / 2,
         q3 = 3 * v.size() / 4;

    std::ostringstream s;

    std::nth_element(v.begin(), v.begin(), v.end());
    s << "[Min: " << std::setw(11) << *v.begin() << " ] ";

    std::nth_element(v.begin(), v.begin() + q1, v.end());
    s << "[Q1: " << std::setw(11) << v[q1] << " ] ";
 
    std::nth_element(v.begin(), v.begin() + q2, v.end());
    s << "[Q2: " << std::setw(11) << v[q2] << " ] ";
    
    std::nth_element(v.begin(), v.begin() + q3, v.end());
    s << "[Q3: " << std::setw(11) << v[q3] << " ] ";    
 
    std::nth_element(v.begin(), v.end()-1, v.end());
    s << "[Max: " << std::setw(11) << *(v.end()-1) << " ] ";

    return s.str();
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
