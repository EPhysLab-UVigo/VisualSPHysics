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

#ifndef OPS_H
#define OPS_H

#include <array>
#include <vector>
#include <string>

/**
   \brief This namespace groups some common math functions.
   This namespace includes several math operations used in some stages of the simulation.
 */
namespace ops {

  /**
     Compute several statistics of a vector.
     The statistics computed are the minimum, maximum and Q1, Q2 and Q3 quartiles.
     \param vec Vector of type double.
     \param std::string with the stats
  */
  std::string vectorStats(std::vector<double> &vec);

  /**
     Substraction component by component of two double arrays of three components.
     \param xi Double array of three components.
     \param xj Double array of three components.
     \return Double array of three components.
   */
  std::array<double,3> substract(std::array<double,3> xi, std::array<double,3> xj);

  /**
     Compute the magnitude of a vector.
     \param x Double array of three components.
     \return Magnitude.
  */
  double magnitude(std::array<double,3> x);

  /**
     Normalizes a vector.
     \param x Double array of three components.
     \return Normalized vector.
  */
  std::array<double,3> normalize(std::array<double,3> x);

  /**
     Computes the dot product.
     \param xi Double array of three components. 
     \param xj Double array of three components.
     \return Dot product.
   */
  double dotProduct(std::array<double,3> xi, std::array<double,3> xj);

  /**
     Computes the distance vector. This is, the substraction divided by the modulus.
     \param xi Double array of three components. 
     \param xj Double array of three components.
     \return Distance vector.
   */
  std::array<double,3> distanceVector(std::array<double,3> xi, std::array<double,3> xj);
}

#endif
