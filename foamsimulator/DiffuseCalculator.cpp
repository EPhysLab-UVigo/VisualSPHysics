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

#define _USE_MATH_DEFINES

#include "DiffuseCalculator.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include <omp.h>

// VTK stuff

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include <vtkSmartPointer.h>

#include <vtkIntArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>

#include "FluidData.h"
#include "PlyWriter.h"
#include "VtkDWriter.h"

#include "BucketContainer.h"

#include "Ops.h"

#define SURFACE 0.75

// Clamping function
#ifndef _MSVC
#pragma omp declare simd
#endif
double DiffuseCalculator::phi(double I, double tmin, double tmax) {
  return (fmin(I, tmax) - fmin(I, tmin)) / (tmax - tmin);
}

//(3) Kernel
double DiffuseCalculator::W(std::array<double, 3> xij, double h) {
  double rval = 0;
  double mxij = ops::magnitude(xij);
  if (mxij <= h)
    rval = 1. - (mxij / h);
  return rval;
}

// Wendland kernel
double DiffuseCalculator::Wwendland(std::array<double, 3> xij, double h) {
  double rval = 0;
  double mxij = ops::magnitude(xij), q = mxij / h;

  if (q >= 0. && q <= 2) {
    double ad = 21. / (16. * M_PI * h * h * h), e1 = (1. - (q / 2.0));
    rval = ad * e1 * e1 * e1 * e1 * (2 * q + 1.);
  }
  return rval;
}

// Poly6 kernel
double DiffuseCalculator::Wpoly6(std::array<double, 3> xij, double h) {
  double rval = 0;
  double mxij = ops::magnitude(xij);
  if (mxij >= 0. && mxij <= h) {
    rval = (315. / (64. * M_PI * pow(h, 9))) * pow((h * h - mxij * mxij), 3);
  }
  return rval;
}

double DiffuseCalculator::vdiff2p(std::array<double, 3> vi,
                                  std::array<double, 3> vj,
                                  std::array<double, 3> xi,
                                  std::array<double, 3> xj, double h) {
  double e1 = ops::magnitude(ops::substract(vi, vj));
  double e2 = 1 - ops::dotProduct(ops::distanceVector(vi, vj),
                                  ops::distanceVector(xi, xj));
  double e3 = W(ops::substract(xi, xj), h);

  return e1 * e2 * e3;
}

// Colorfield between two particles
double DiffuseCalculator::colorField2p(std::array<double, 3> xi,
                                       std::array<double, 3> xj, double h,
                                       double mj, double pj) {
  return (mj / pj) * Wwendland(ops::substract(xi, xj), h);
}

// Smoothed gradient field of the smoothed color field
std::array<double, 3> DiffuseCalculator::gradient2p(std::array<double, 3> xi,
                                                    std::array<double, 3> xj,
                                                    double h, double csi,
                                                    double csj) {
  double wval = Wwendland(ops::substract(xi, xj), h);
  return std::array<double, 3>{{wval * csj * (xi[0] - xj[0]),
                                wval * csj * (xi[1] - xj[1]),
                                wval * csj * (xi[2] - xj[2])}};
}

// Surface curvature for 2 particles
double DiffuseCalculator::curvature2p(std::array<double, 3> xi,
                                      std::array<double, 3> xj,
                                      std::array<double, 3> ni,
                                      std::array<double, 3> nj, double h) {
  double e1 = 1 - ops::dotProduct(ops::normalize(ni), ops::normalize(nj));
  double e2 = W(ops::substract(xi, xj), h);
  return e1 * e2;
}

// Wave crests
double DiffuseCalculator::crests2p(std::array<double, 3> xi,
                                   std::array<double, 3> xj,
                                   std::array<double, 3> vi,
                                   std::array<double, 3> ni,
                                   std::array<double, 3> nj, double h) {
  std::array<double, 3> xji = ops::distanceVector(xj, xi),
                        nni = ops::normalize(ni), nvi = ops::normalize(vi);
  double kij = 0;
  if (ops::dotProduct(xji, nni) < 0 && ops::dotProduct(nvi, nni) >= 0.6)
    kij = curvature2p(xi, xj, ni, nj, h);
  return kij;
}

double DiffuseCalculator::solveEq(double px, double py, double pz, double vx,
                                  double vy, double vz, double x, double y) {
  return ((-(x - px) * vx - (y - py) * vy) / vz) + pz;
}

// Constructor
DiffuseCalculator::DiffuseCalculator(SimulationParams p) : sp(p) {}

void DiffuseCalculator::runSimulation() {
  std::string seqnum(sp.nzeros, '0'),
      formats = std::string("%.") + std::to_string(sp.nzeros) + std::string("d");

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> xunif(0, 1);

  long difId = 0;

  // Persistent particle vector
  std::vector<std::array<double, 3>> ppPosit, ppVel;
  std::vector<int> ppIds, ppTTL;
  std::vector<double> ppDensity;

  // Let's loop!
  for (int nstep = sp.nstart; nstep <= sp.nend; nstep++) {

    std::sprintf(&seqnum[0], formats.c_str(), nstep);
    std::string fileName = fs::path(sp.dataPath) / (sp.filePrefix + seqnum + ".vtk");

    std::cout << "\n\n=====================================================================\n";
    std::cout << "Opening: " << fileName << std::endl;

    FluidData file(sp.MINX, sp.MAXX, sp.MINY, sp.MAXY, sp.MINZ, sp.MAXZ, sp.h);

    if (std::string(sp.exclusionZoneFile) != "")
      file.setExclusionZone(sp.exclusionZoneFile);
    file.loadFile(fileName);

    BucketContainer<particle> &f = *(file.getBucketContainer());

    long npoints =
        f.getNElements(); // output->GetPoints()->GetNumberOfPoints();

    // Create a vector with the scaled velocity difference for each particle
    std::vector<double> Ita(npoints, 0.0);
    std::vector<double> colorField(npoints, 0.0);
    std::vector<double> waveCrest(npoints, 0.0);
    std::vector<double> energy(npoints, 0.0);
    std::vector<int> ndiffuse(npoints, 0);
    std::vector<std::array<double, 3>> gradient(npoints, std::array<double, 3>{{0, 0, 0}});

    std::cout << "Total fluid particles: " << npoints << std::endl;

    std::cerr << "\n[Stage 1] trapped air potential, Energy and colorfield. ";

    double wall_timer = omp_get_wtime();
    auto &buckets = f.getNoEmptyBuckets();

    /*
     * First pass: trapped air potential, Energy and colorfield
     */
    {
#pragma omp parallel for schedule(guided)
      for (long nebucket = 0; nebucket < buckets.size();
           nebucket++) { // Iterate over all buckets
        auto &bucket = buckets[nebucket].second;
        auto sbuckets = f.getSurroundingBuckets(buckets[nebucket].first);

        for (auto &pi : bucket) { // Iterate over each particle in the bucket
          long i = pi.id;
          auto vi = pi.vel, xi = pi.pos;

          for (auto sb : sbuckets) { // Iterate over surrounding buckets
            for (auto &pj : *sb) {   // Iterate over each particle in the bucket

              if (pi.id != pj.id) {
                auto vj = pj.vel, xj = pj.pos;

                // Substract position
                double spx = xi[0] - xj[0], spy = xi[1] - xj[1],
                       spz = xi[2] - xj[2];

                double mp = sqrt(spx * spx + spy * spy + spz * spz);
                double q = mp / sp.h;

                if (mp <= sp.h) {
                  // Substract velocity
                  double svx = vi[0] - vj[0], svy = vi[1] - vj[1],
                         svz = vi[2] - vj[2];

                  // Magnitude
                  double mv = sqrt(svx * svx + svy * svy + svz * svz);

                  // Distance vector
                  double dvx = svx / mv, dvy = svy / mv, dvz = svz / mv;

                  double dpx = spx / mp, dpy = spy / mp, dpz = spz / mp;

                  double e = 1 - (dvx * dpx + dvy * dpy + dvz * dpz);

                  double w = 1. - q;

                  Ita[i] += mv * e * w;
                }

                if (q >= 0 && q <= 2) {
                  double ad = 21. / (16. * M_PI * sp.h * sp.h * sp.h),
                         e1 = (1. - (q / 2.0));
                  colorField[i] += (sp.mass / pj.rhop) * ad * e1 * e1 * e1 *
                                   e1 * (2 * q + 1.);
                }
              }
            }
          }

          energy[i] =
              0.5 * sp.mass * (vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]);
        }
      }
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    std::cerr << "[Stage 2] gradient. ";
    /*
     * Second pass: gradient
     */
    {
      wall_timer = omp_get_wtime();

#pragma omp parallel for schedule(guided)
      for (long nebucket = 0; nebucket < buckets.size();
           nebucket++) { // Iterate over all buckets
        auto &bucket = buckets[nebucket].second;
        auto sbuckets = f.getSurroundingBuckets(buckets[nebucket].first);

        for (auto &pi : bucket) { // Iterate over each particle in the bucket
          long i = pi.id;

          for (auto sb : sbuckets) { // Iterate over surrounding buckets
            for (auto &pj : *sb) {   // Iterate over each particle in the bucket
              auto xij = ops::substract(pi.pos, pj.pos);
              double mxij = ops::magnitude(xij), q = mxij / sp.h;

              if (q >= 0 && q <= 2) {
                double ad = 21. / (16. * M_PI * sp.h * sp.h * sp.h),
                       e1 = (1. - (q / 2.0)),
                       rval = colorField[pj.id] * ad * e1 * e1 * e1 * e1 *
                              (2 * q + 1.);
                gradient[i][0] += rval * xij[0];
                gradient[i][1] += rval * xij[1];
                gradient[i][2] += rval * xij[2];
              }
            }
          }
        }
      }
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    std::cerr << "[Stage 3] wave crests. ";
    /*
     * Third pass: wave crests
     */
    {

      wall_timer = omp_get_wtime();

#pragma omp parallel for schedule(guided)
      for (long nebucket = 0; nebucket < buckets.size(); nebucket++) { // Iterate over all buckets
        auto &bucket = buckets[nebucket].second;
        std::vector<std::vector<particle> *> sbuckets;

        for (auto &pi : bucket) { // Iterate over each particle in the bucket
          long i = pi.id;
          if (colorField[i] < SURFACE) {
            if (sbuckets.size() == 0)
              sbuckets = f.getSurroundingBuckets(buckets[nebucket].first);
            for (auto sb : sbuckets) { // Iterate over surrounding buckets
              for (auto &pj : *sb) { // Iterate over each particle in the bucket
                waveCrest[i] += crests2p(pi.pos, pj.pos, pi.vel, gradient[i],
                                         gradient[pj.id], sp.h);
              }
            }
          }
        }
      }
      std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;
    }

    std::cout << "\nStatistics " << std::endl << "===========" << std::endl;

    std::cout << "*** Wave crests ***" << std::endl;
    ops::printVectorStats(waveCrest);
    std::cout << "*** Trapped air ***" << std::endl;
    ops::printVectorStats(Ita);
    std::cout << "*** Energy ***" << std::endl;
    ops::printVectorStats(energy);

    std::cerr << "\n[Stage 4] clamping function. ";
    /*
     * Fourth pass: clamping function
     */
    wall_timer = omp_get_wtime();

#ifndef _MSVC
#pragma omp parallel for simd
#else
#pragma omp parallel for schedule(static)
#endif
    for (long i = 0; i < npoints; i++) {
      waveCrest[i] = phi(waveCrest[i], sp.MINWC, sp.MAXWC);
      Ita[i] = phi(Ita[i], sp.MINTA, sp.MAXTA);
      energy[i] = phi(energy[i], sp.MINK, sp.MAXK);
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    long npdiffuse = 0;

    std::cerr << "[Stage 5] number of diffuse particles generated: ";
    /*
     * Fifth pass: number of diffuse particles generated
     */
    wall_timer = omp_get_wtime();
#ifndef _MSVC
#pragma omp parallel for simd reduction(+ : npdiffuse)
#endif
    for (long i = 0; i < npoints; i++) {
      ndiffuse[i] = std::floor(
          energy[i] * (sp.KTA * Ita[i] + sp.KWC * waveCrest[i]) * sp.TIMESTEP);
      npdiffuse += ndiffuse[i];
    }

    std::cout << npdiffuse << " Time: " << omp_get_wtime() - wall_timer
              << std::endl;

    std::cerr << "[Stage 6] calculate diffuse particle positions. ";
    /*
     * Sixth pass: calculate diffuse particle positions
     */

    // Diffuse particle vector!
    std::vector<std::array<double, 3>> diffusePosit(npdiffuse);
    std::vector<std::array<double, 3>> diffuseVel(npdiffuse);
    std::vector<int> diffuseIds(npdiffuse), diffuseTTL(npdiffuse);
    std::vector<double> diffuseDensity(npdiffuse, 0.0);

		// Generate random numbers: this is done out of the loop because it is not thread safe
		std::vector<double> tempRand(npdiffuse * 3);
		for (auto &x : tempRand)
			x = xunif(gen);

    wall_timer = omp_get_wtime();

    {
      long idif = 0;

#pragma omp parallel for schedule(guided)
      for (long nebucket = 0; nebucket < buckets.size(); nebucket++) { // Iterate over all buckets
        auto &bucket = buckets[nebucket].second;

        for (auto &pi : bucket) { // Iterate over each particle in the bucket
          long i = pi.id;

          if (ndiffuse[i] >= 1) {
            std::array<double, 3> pos = pi.pos, vel = pi.vel;

            // Obtain orthogonal vectors to velocity vector
            std::array<double, 3> e1, e2;

            // Find non-zero component of velocity vector in order to avoid
            // division by 0 and calculate e1
            if (vel[0] != 0) { // x non zero

              e1 = ops::normalize({{solveEq(pos[2], pos[1], pos[0], vel[2],
                                            vel[1], vel[0], 0, 1),
                                    1, 0}});
            } else if (vel[1] != 0) { // y non zero
              e1 = ops::normalize({{1,
                                    solveEq(pos[0], pos[2], pos[1], vel[0],
                                            vel[2], vel[1], 1, 0),
                                    0}});
            } else { // z non zero
              e1 = ops::normalize({{1, 0,
                                    solveEq(pos[0], pos[1], pos[2], vel[0],
                                            vel[1], vel[2], 1, 0)}});
            }

            // Cross product of two orthogonal vectors generate a vector
            // orthogonal to them
            e2 = ops::normalize({{e1[1] * vel[2] - vel[1] * e1[2],
                                  e1[0] * vel[2] - vel[0] * e1[2],
                                  e1[0] * vel[1] - vel[0] * e1[1]}});

            std::array<double, 3> nvel =
                ops::normalize(std::array<double, 3>{{vel[0], vel[1], vel[2]}});

            for (int j = 0; j < ndiffuse[i]; j++) {
              double h = tempRand[idif * 3] *
                         (ops::magnitude(std::array<double, 3>{{vel[0], vel[1], vel[2]}}) * sp.TIMESTEP) *
												 .5,
                     r = sp.h * sqrt(tempRand[idif * 3 + 1]), theta = tempRand[idif * 3 + 2] * 2 * M_PI;

              // Position of newly created diffuse particle
              diffusePosit[idif] = {{pos[0] + r * cos(theta) * e1[0] +
                                         r * sin(theta) * e2[0] + h * nvel[0],
                                     pos[1] + r * cos(theta) * e1[1] +
                                         r * sin(theta) * e2[1] + h * nvel[1],
                                     pos[2] + r * cos(theta) * e1[2] +
                                         r * sin(theta) * e2[2] + h * nvel[2]}};

              // Velocity of newly created diffuse particle
              diffuseVel[idif] = {
                  {r * cos(theta) * e1[0] + r * sin(theta) * e2[0] + vel[0],
                   r * cos(theta) * e1[1] + r * sin(theta) * e2[1] + vel[1],
                   r * cos(theta) * e1[2] + r * sin(theta) * e2[2] + vel[2]}};

              // Particle ID
              diffuseIds[idif] = difId;

              // Particle lifetime
              diffuseTTL[idif] = ndiffuse[i] * sp.LIFEFIME;

              difId++;
              idif++;
            }
          }
        }
      }
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    // Seventh pass: classify particles
    //[0-6]Spray [6-20]Foam [20..]Bubbles ¿?
    std::cerr << "[Stage 7] classify particles. ";
    wall_timer = omp_get_wtime();

#pragma omp parallel for schedule(guided)
    for (long i = 0; i < npdiffuse; i++) {
      auto pxd = diffusePosit[i];
      auto sbuckets = f.getSurroundingBuckets(pxd);
      for (auto sb : sbuckets) { // Iterate over surrounding buckets
        for (auto &pj : *sb) {   // Iterate over each particle in the bucket
          if (ops::magnitude(ops::substract(pxd, pj.pos)) <= sp.h) {
            diffuseDensity[i]++;
          }
        }
      }
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    // Update particles

    std::cerr << "[Stage 8] Update particles! ";

    wall_timer = omp_get_wtime();

#pragma omp parallel for schedule(guided)
    for (long i = 0; i < ppIds.size(); i++) {
      auto &pxd = ppPosit[i];
			auto temppos = ppPosit[i];
      std::array<double, 3> num{{0, 0, 0}};
      double den = 0;

			// Recalculate density: should be placed before the new position calculation.
      ppDensity[i] = 0;
      for (auto sb : f.getSurroundingBuckets(pxd)) { // Iterate over surrounding buckets
        for (auto &pj : *sb) { // Iterate over each particle in the bucket
          if (ops::magnitude(ops::substract(pxd, pj.pos)) <= sp.h) {
            ppDensity[i]++;
          }
        }
      }

			if (ppDensity[i] >= sp.SPRAY) { // This is not needed for spray particles.
				std::vector<std::vector<particle> *> sbuckets = f.getSurroundingBuckets(pxd);
				for (auto sb : sbuckets) { // Iterate over surrounding buckets
					for (auto &pj : *sb) {   // Iterate over each particle in the bucket
						double tval = Wwendland(ops::substract(pxd, pj.pos), sp.h);
						num = {{num[0] + pj.vel[0] * tval,
										num[1] + pj.vel[1] * tval,
										num[2] + pj.vel[2] * tval}};
						den += tval;
					}
				}
			}

      // Now we can re-clasify and calculate new positions
      if (ppDensity[i] < sp.SPRAY) { // It's spray!
        // TODO: we are avoiding external forces (like wind)
        ppVel[i] = {{ppVel[i][0], ppVel[i][1], ppVel[i][2] + -9.81 * sp.TIMESTEP}};
        pxd = {{pxd[0] + sp.TIMESTEP * ppVel[i][0],
                pxd[1] + sp.TIMESTEP * ppVel[i][1],
                pxd[2] + sp.TIMESTEP * ppVel[i][2]}};

      } else if (ppDensity[i] > sp.BUBBLES) { // It's a bubble!
        num = {{num[0] / den, num[1] / den, num[2] / den}};
        ppVel[i] = {{ppVel[i][0] + sp.TIMESTEP * (sp.KD * (num[0] - ppVel[i][0]) / sp.TIMESTEP),
		    						 ppVel[i][1] + sp.TIMESTEP * (sp.KD * (num[1] - ppVel[i][1]) / sp.TIMESTEP),
		     						 ppVel[i][2] + sp.TIMESTEP * (-sp.KB * -9.81 + sp.KD * (num[2] - ppVel[i][2]) / sp.TIMESTEP)}};
        pxd = {{pxd[0] + sp.TIMESTEP * ppVel[i][0],
                pxd[1] + sp.TIMESTEP * ppVel[i][1],
                pxd[2] + sp.TIMESTEP * ppVel[i][2]}};

      } else { // It's foam!
        num = {{num[0] / den, num[1] / den, num[2] / den}};
        ppVel[i] = {{num[0], num[1], num[2]}};
	
        pxd = {{pxd[0] + sp.TIMESTEP * num[0],
								pxd[1] + sp.TIMESTEP * num[1],
                pxd[2] + sp.TIMESTEP * num[2]}};
      }
    }

    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;

    // Delete particles
    std::cerr << "[Stage 9] Delete particles! ";

    std::vector<std::array<double, 3>> tempPosit, tempVel;
    std::vector<int> tempIds, tempTTL;
    std::vector<double> tempDensity;

    wall_timer = omp_get_wtime();

    // #pragma omp parallel for schedule(auto)
    for (long i = 0; i < ppIds.size(); i++) {
      // Decrease TTL for foam particles
      if (ppDensity[i] > sp.SPRAY && ppDensity[i] < sp.BUBBLES)
        ppTTL[i]--;

      // If TTL is less than zero delete particle
      // OR If particle is out of the domain, delete particle
      if (!(ppTTL[i] < 0 || 
						ppPosit[i][0] <= sp.MINX || ppPosit[i][1] <= sp.MINY || ppPosit[i][2] <= sp.MINZ || 
						ppPosit[i][0] >= sp.MAXX || ppPosit[i][1] >= sp.MAXY || ppPosit[i][2] >= sp.MAXZ)) {
        tempPosit.push_back(ppPosit[i]);
        tempVel.push_back(ppVel[i]);
        tempIds.push_back(ppIds[i]);
        tempTTL.push_back(ppTTL[i]);
        tempDensity.push_back(ppDensity[i]);
      }
    }

		ppIds = std::move(tempIds);
		ppPosit = std::move(tempPosit);
		ppVel = std::move(tempVel);
		ppDensity = std::move(tempDensity);
		ppTTL = std::move(tempTTL);

    std::cout << "Deleted: " << ppIds.size() - tempIds.size();
    std::cout << " Time: " << omp_get_wtime() - wall_timer << std::endl;

    // Append new particles
    std::cerr << "[Stage 10] Append new particles! Total diffuse particles: "
              << ppIds.size() << std::endl;

    if (npdiffuse > 0) {
      std::copy(diffuseIds.begin(), diffuseIds.end(), std::back_inserter(ppIds));
      std::copy(diffusePosit.begin(), diffusePosit.end(), std::back_inserter(ppPosit));
      std::copy(diffuseVel.begin(), diffuseVel.end(), std::back_inserter(ppVel));
      std::copy(diffuseDensity.begin(), diffuseDensity.end(), std::back_inserter(ppDensity));
      std::copy(diffuseTTL.begin(), diffuseTTL.end(), std::back_inserter(ppTTL));
    }

    std::cout << "\n*** Density ***" << std::endl;
    ops::printVectorStats(ppDensity);

    /*
     * Write diffuse particle files
     */
    std::cout << "\n[Stage 11] Save to file! ";

    wall_timer = omp_get_wtime();

#ifndef _MSVC
#pragma omp parallel sections
#endif
    {

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.text_files) {
        // Save diffuse particles to simple text files
        std::string outFilename = fs::path(sp.outputPath) / (sp.outputPreffix + seqnum + ".txt");
        std::ofstream tfile(outFilename, std::ios::trunc);
        tfile.setf(std::ios::scientific);

        for (long i = 0; i < ppIds.size(); i++) {
          tfile << ppPosit[i][0] << " " << ppPosit[i][1] << " "
                << ppPosit[i][2];
          int ttype = 1;
          if (ppDensity[i] < sp.SPRAY) {
            ttype = 0;
          } else if (ppDensity[i] > sp.BUBBLES) {
            ttype = 2;
          }
          tfile << " " << ttype << std::endl;
        }
        tfile.close();
      }

      // Generate small icosphere per each diffuse particle and save to ply file
      // 3 ply files per step will be generated one for foam, one spray and a
      // third one for bubbles

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.ply_files) {
        std::vector<std::array<double, 3>> material;
        std::string plyFilename = fs::path(sp.outputPath) / (sp.outputPreffix + seqnum +
                                  ".spray.ply"),
                    comment = "Case: " + sp.outputPreffix +
                              " Frame: " + seqnum + " Spray";
        PlyWriter output(plyFilename, sp.MINX, sp.MAXX, sp.MINY, sp.MAXY,
                         sp.MINZ, sp.MAXZ, sp.h);
        output.setComment(comment);

        for (long i = 0; i < ppIds.size(); i++)
          if (ppDensity[i] < sp.SPRAY)
            material.push_back(ppPosit[i]);

        output.setData(&material);
        output.write();
      }

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.ply_files) {
        std::vector<std::array<double, 3>> material;
        std::string plyFilename =
                        fs::path(sp.outputPath) / (sp.outputPreffix + seqnum + ".foam.ply"),
                    comment = "Case: " + sp.outputPreffix +
                              " Frame: " + seqnum + " Foam";
        PlyWriter output(plyFilename, sp.MINX, sp.MAXX, sp.MINY, sp.MAXY,
                         sp.MINZ, sp.MAXZ, sp.h);
        output.setComment(comment);

        for (long i = 0; i < ppIds.size(); i++)
          if (ppDensity[i] > sp.SPRAY && ppDensity[i] < sp.BUBBLES)
            material.push_back(ppPosit[i]);

        output.setData(&material);
        output.write();
      }

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.ply_files) {
        std::vector<std::array<double, 3>> material;
        std::string plyFilename = fs::path(sp.outputPath) / (sp.outputPreffix + seqnum +
                                  ".bubbles.ply"),
                    comment = "Case: " + sp.outputPreffix +
                              " Frame: " + seqnum + " Bubbles";
        PlyWriter output(plyFilename, sp.MINX, sp.MAXX, sp.MINY, sp.MAXY,
                         sp.MINZ, sp.MAXZ, sp.h);
        output.setComment(comment);

        for (long i = 0; i < ppIds.size(); i++)
          if (ppDensity[i] > sp.BUBBLES)
            material.push_back(ppPosit[i]);

        output.setData(&material);
        output.write();
      }

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.vtk_files) {
        std::string vtkFilename =
            fs::path(sp.outputPath) / (sp.outputPreffix + seqnum + ".vtk");
        VtkDWriter output(vtkFilename, sp.MINX, sp.MAXX, sp.MINY, sp.MAXY,
                          sp.MINZ, sp.MAXZ, sp.h);
        output.setData(&ppPosit, &ppVel);
        output.write();
      }

#ifndef _MSVC
#pragma omp section
#endif
      if (sp.vtk_diffuse_data) {

        // Save diffuse particle data
        vtkSmartPointer<vtkPoints> dpoints = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkDoubleArray> dvels =
            vtkSmartPointer<vtkDoubleArray>::New();
        dvels->SetName("Velocity");
        dvels->SetNumberOfComponents(3);

        vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
        ids->SetName("id");

        vtkSmartPointer<vtkIntArray> ptype =
            vtkSmartPointer<vtkIntArray>::New();
        ptype->SetName("ParticleType");

        vtkSmartPointer<vtkDoubleArray> density =
            vtkSmartPointer<vtkDoubleArray>::New();
        density->SetName("Density");

        for (long i = 0; i < ppIds.size(); i++) {
          ids->InsertNextValue(ppIds[i]);
          int ttype = 1;
          if (ppDensity[i] < sp.SPRAY) {
            ttype = 0;
          } else if (ppDensity[i] > sp.BUBBLES) {
            ttype = 2;
          }
          ptype->InsertNextValue(ttype);
          density->InsertNextValue(ppDensity[i]);
          dpoints->InsertNextPoint(ppPosit[i].data());
          dvels->InsertNextTuple(ppVel[i].data());
        }

        vtkSmartPointer<vtkPolyData> difpolydata =
            vtkSmartPointer<vtkPolyData>::New();
        difpolydata->SetPoints(dpoints);

        vtkSmartPointer<vtkCellArray> vertices =
            vtkSmartPointer<vtkCellArray>::New();
        for (long i = 0; i < dpoints->GetNumberOfPoints(); ++i) {
          vtkIdType pt[] = {i};
          vertices->InsertNextCell(1, pt);
        }
        difpolydata->SetVerts(vertices);

        difpolydata->GetPointData()->SetScalars(ids);
        difpolydata->GetPointData()->AddArray(ptype);
        difpolydata->GetPointData()->AddArray(dvels);
        difpolydata->GetPointData()->AddArray(density);

        std::string outFilename =
            fs::path(sp.outputPath) / (sp.outputPreffix + seqnum + "_diffuse.vtk");

        vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
				writer->SetFileTypeToBinary();
        writer->SetFileName(outFilename.c_str());
        writer->SetInputData(difpolydata);
        writer->Write();
      }

      /*
       * Write intermediary files
       */
#ifndef _MSVC
#pragma omp section
#endif
      if (sp.vtk_fluid_data) {

        vtkSmartPointer<vtkPoints> ppoints = vtkSmartPointer<vtkPoints>::New();

        // int nbucket = 0;

        for (auto &bucket : f.getBuckets()) { // Iterate over surrounding buckets
          for (auto &pi : bucket) { // Iterate over each particle in the bucket
            ppoints->InsertNextPoint(pi.pos.data());
          }
        }

        vtkSmartPointer<vtkDoubleArray> vTrappedAir =
            vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> vCrests =
            vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> vEnergy =
            vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> vDiffuse =
            vtkSmartPointer<vtkDoubleArray>::New();

        vTrappedAir->SetName("TrappedAir");
        vCrests->SetName("WaveCrests");
        vEnergy->SetName("Energy");
        vDiffuse->SetName("DiffuseParticles");

        for (long i = 0; i < npoints; i++) {
          vTrappedAir->InsertNextValue(Ita[i]);
          vCrests->InsertNextValue(waveCrest[i]);
          vEnergy->InsertNextValue(energy[i]);
          vDiffuse->InsertNextValue(ndiffuse[i]);
        }

        vtkSmartPointer<vtkPolyData> ppolydata =
            vtkSmartPointer<vtkPolyData>::New();
        ppolydata->SetPoints(ppoints);

        vtkSmartPointer<vtkCellArray> vertices =
            vtkSmartPointer<vtkCellArray>::New();
        for (long i = 0; i < ppoints->GetNumberOfPoints(); ++i) {
          vtkIdType pt[] = {i};
          vertices->InsertNextCell(1, pt);
        }
        ppolydata->SetVerts(vertices);

        ppolydata->GetPointData()->AddArray(vTrappedAir);
        ppolydata->GetPointData()->AddArray(vCrests);
        ppolydata->GetPointData()->AddArray(vEnergy);
        ppolydata->GetPointData()->AddArray(vDiffuse);

        vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
        std::string outFilename =
            fs::path(sp.outputPath) / (sp.outputPreffix + seqnum + "_fluid.vtk");
				writer->SetFileTypeToBinary();
        writer->SetFileName(outFilename.c_str());
        writer->SetInputData(ppolydata);
        writer->Write();
      }
    }
    std::cout << "Time: " << omp_get_wtime() - wall_timer << std::endl;
  }
}
