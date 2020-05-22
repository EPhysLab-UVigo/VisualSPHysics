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

#ifndef SIMULATIONPARAMS_H
#define SIMULATIONPARAMS_H

#include <string>

/**
   \brief Parámetros para iniciar una simulación.
 */
struct SimulationParams {

  std::string dataPath,			            ///< Path of the input files.
    filePrefix, 			                  ///< Prefix of the input file names.
    outputPath, 			                  ///< Path of the output files.
    outputPreffix, 			                ///< Prefix of the output file names.
    exclusionZoneFile;                  ///< FIle with the exclusion zone geometry.
  
  int nstart, 				                  ///< Initial time of the simulation.
    nend,				                        ///< Ending simulation time.
    nzeros,				                      ///< Padding of the time step string.
    text_files,                         ///< Points if output files in text format are enabled. Includes position and type of each particle.
    vtk_files,                          ///< Points if output files in Vtk format are enabled. Includes the position and the size of the particles.
    vtk_diffuse_data,                   ///< Points if output files in Vtk format are enabled. Includes position, velocity, id, type and density of the diffuse particles.
    vtk_fluid_data;                     ///< Points if output files in Vtk format are enabled. Includes information of the fluid particles for each time step of the simulation.

  double h,				                      ///< H value in meters.
    mass,                               ///< Mass of each fluid particle in Kg.
    TIMESTEP,                           ///< Time step in seconds.

    MINX,                               ///< Domain limits. Minimum in the axis x.
    MINY,                               ///< Domain limits. Minimum in the axis y.
    MINZ,                               ///< Domain limits. Minimum in the axis z.
    MAXX,                               ///< Domain limits. Maximum in the axis x.
    MAXY,                               ///< Domain limits. Maximum in the axis y.
    MAXZ,                               ///< Domain limits. Maximum in the axis z.
    
    MINTA,                              ///< Minimum threshold for trapped air.
    MAXTA,                              ///< Maximum threshold for trapped air.

    MINWC,                              ///< Minimum threshold for wave crests.
    MAXWC,                              ///< Maximum threshold for wave crests.

    MINK,                               ///< Minimum threshold for kinetic energy.
    MAXK,                               ///< Maximum threshold for kinetic energy.

    KTA,  				                      ///< Amount of diffuse material created by trapper air.
    KWC, 				                        ///< Amount of diffuse material created by wave crests.
    SPRAY,                              ///< Maximum density of fluid for spray particles.
    BUBBLES,                            ///< Minumum density of fluid for bubble particles.
    LIFEFIME, 				                  ///< Life time of diffuse particles.
    KB, 				                        ///< Buoyancy factor for bubble particles.
    KD; 				                        ///< Drag factor for buoyancy particles.
};

#endif
