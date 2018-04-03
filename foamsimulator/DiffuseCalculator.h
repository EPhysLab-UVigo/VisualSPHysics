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

#ifndef DIFFUSECALCULATOR_H
#define DIFFUSECALCULATOR_H

#include <string>
#include <array>
#include "SimulationParams.h"

/**
   \brief This class provides the main functionality to compute the foam simulation.

   This class simulates the foam generation from a DualSPHysics simulation and some
   parameters given by the user that adjust the foam generation thresholds. The output files
   can be written in several formats.
 */
class DiffuseCalculator {

 public:
  /**
     Class constructor.
     \param p Simulation parameters.
     \see SimulationParams
   */
  DiffuseCalculator(SimulationParams p);

  /**
     Starts the simulation.
   */
  void runSimulation();
  
 private:
  SimulationParams sp;

  /**
     Maps a value I between zero and one according to a max and min thresolds.
     \param I Value to map.
     \param tmin Min threshold.
     \param tmax Max threshold.
     \return Mapped value.
   */
  double phi(double I, double tmin, double tmax);

  /**
     Implementation of radially symmetrical kernel. This function weighs the given value according to the distance.
     It is simpler than the cubic spline or Wendland, but provides good results some cases according to Ihmsen et al.
     \param xij Value to weigh.
     \param h Smoothing length.
     \return Weighted value.
   */
  double W(std::array<double,3> xij, double h);

  /**
     Implementation of a Wendland kernel. This function weighs the given value according to the distance.
     Usually provides better results than the classic cubic spline.
     \param xij Value to weigh.
     \param h Smoothing length.
     \return Weighted value.
   */
  double Wwendland(std::array<double,3> xij, double h);

  /**
     Implementation of a "poly6" kernel. This function weighs the given value according to the distance. Created by Müller et al.
     \param xij Value to weigh.
     \param h Smoothing length.
     \return Weighted value.
  */
  double Wpoly6(std::array<double,3> xij, double h);

  /**
     Computes the scaled velocity difference between two particles, necessary to obtain the trapped air potential.
     \param vi Velocity vector of the particle i.
     \param vj Velocity vector of the particle j.
     \param xi Position of the particle i.
     \param xj Position of the particle j.
     \param h Smoothing length.
     \return Scaled velocity difference between two particles.
   */
  double vdiff2p(std::array<double,3> vi, std::array<double,3> vj, 
		 std::array<double,3> xi, std::array<double,3> xj, double h);
  
  /**
     Computes the smoothed color field for two particles.
     \param xi Position of the particle i.
     \param xj Position of the particle j.
     \param h Smoothing length.
     \param mj Mass of the particle j.
     \param pj Density of the particle j.
     \return Smoothed colorfield for two particles.
   */
  double colorField2p(std::array<double,3> xi, std::array<double,3> xj, double h, double mj, double pj);

  /**
     Computes the gradient of the smoothed color field between two particles.
     \param xi Position of the particle i.
     \param xj Position of the particle j.
     \param h Smoothing length.
     \param csi Color field for the particle i.
     \param csj Color field for the particle j.
     \return Gradient of the smoothed color field between two particles.
  */
  std::array<double,3> gradient2p(std::array<double,3> xi, std::array<double,3> xj, double h, double csi, double csj);

  /**
     Computes the curvature factor between two particles.
     \param xi Position of the particle i.
     \param xj Position of the particle j.
     \param ni Normalized surface normal of the particle i.
     \param nj Normalized surface normal of the particle j.
     \param h Smoothing length.
     \return Curvature factor between two particles.
   */
  double curvature2p(std::array<double,3> xi, std::array<double,3> xj,
		     std::array<double,3> ni, std::array<double,3> nj, double h);

  /**
     Computes an estimation of wave crest between two particles.
     \param xi Position of the particle i.
     \param xj Position of the particle j.
     \param vi Velocity vector of the particle i.
     \param ni Normalized surface normal of the particle i.
     \param nj Normalized surface normal of the particle j.
     \param h Smoothing length.
     \return Wave crest factor between two particles.
   */
  double crests2p(std::array<double,3> xi, std::array<double,3> xj, std::array<double,3> vi, 
		  std::array<double,3> ni, std::array<double,3> nj, double h);

  /**
     Computes the third component of an orthogonal vector to a vector "v" that crosses a point p.
     \param px Component x of the point p.
     \param py Component y of the point p.
     \param pz Component z of the point p.
     \param vx Component x of the vector v.
     \param vy Component y of the vector v.
     \param vz Component z of the vector v.
     \param x First component of the normal vector.
     \param y Second component of the normal vector.
     \return Third component of a vector.
  */
  double solveEq(double px, double py, double pz,
		 double vx, double vy, double vz, 
		 double x, double y);
};

#endif
