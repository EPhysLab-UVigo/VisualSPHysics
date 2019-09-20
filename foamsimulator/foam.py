#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# VisualSPHysics
# Copyright (C) 2019 Orlando Garcia-Feal orlando@uvigo.es

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


## This simple Python program runs the foam simulation. It reads the
# user parameters from a ".ini" file provided by parameter.

import configparser
import sys
import xml.etree.ElementTree as ET

#Change this for the location of foam simulator module:
sys.path.append("/path/to/VisualSPHysics/build")

import diffuseparticles

#
# Main
#

if len(sys.argv) == 1:
    print("Usage: ./foam.py [config_file.ini]")
    exit()

config = configparser.ConfigParser()
suc = config.read(sys.argv[1])

if len(suc) == 0 :
    print("Error opening "+ sys.argv[1])
    exit()

try:
    # Read PATHS
    paths = config['PATHS']
    InputDataPath = paths['InputDataPath']
    InputFilesPrefix = paths['InputFilesPrefix']
    ZeroPadding = paths.getint('ZeroPadding')
    OutputDataPath = paths['OutputDataPath']
    FilesOutputPrefix = paths['FilesOutputPrefix']
    ExclusionZoneFile = paths['ExclusionZoneFile']
    XmlFile = paths['XmlFile']

    # Read OUTPUT
    ou = config['OUTPUT']

    TextFiles = ou.getboolean('TextFiles')
    PlyFiles = ou.getboolean('PlyFiles')
    VtkFiles = ou.getboolean('VtkFiles')
    VtkDiffuseData = ou.getboolean('VtkDiffuseData')
    VtkFluidData = ou.getboolean('VtkFluidData')
    
    # Read TIMESTEPS
    ts = config['TIMESTEPS']
    
    StartingTimeStep = ts.getint('StartingTimeStep')
    EndingTimeStep = ts.getint('EndingTimeStep')

    # Read FOAMPARAMETERS
    fp = config['FOAMPARAMETERS']

    MinTrappedAirThreshold = fp.getfloat('MinTrappedAirThreshold')
    MaxTrappedAirThreshold = fp.getfloat('MaxTrappedAirThreshold')
    MinWaveCrestsThreshold = fp.getfloat('MinWaveCrestsThreshold')
    MaxWaveCrestsThreshold = fp.getfloat('MaxWaveCrestsThreshold')
    MinKineticEnergyThreshold = fp.getfloat('MinKineticEnergyThreshold')
    MaxKineticEnergyThreshold = fp.getfloat('MaxKineticEnergyThreshold')
    DiffuseTrappedAirMultiplier = fp.getfloat('DiffuseTrappedAirMultiplier')
    DiffuseWaveCrestsMultiplier = fp.getfloat('DiffuseWaveCrestsMultiplier')
    SprayDensity = fp.getfloat('SprayDensity')
    BubblesDensity = fp.getfloat('BubblesDensity')
    LifefimeMultiplier = fp.getfloat('LifefimeMultiplier')
    BuoyancyControl = fp.getfloat('BuoyancyControl')
    DragControl = fp.getfloat('DragControl')

    # Read DOMAIN
    do = config['DOMAIN']

    CustomDomain = do.getboolean('CustomDomain')

    if CustomDomain:        
        DomainMinx = do.getfloat('DomainMinx')
        DomainMiny = do.getfloat('DomainMiny')
        DomainMinz = do.getfloat('DomainMinz')
        DomainMaxx = do.getfloat('DomainMaxx')
        DomainMaxy = do.getfloat('DomainMaxy')
        DomainMaxz = do.getfloat('DomainMaxz')
        if ( DomainMinx >= DomainMaxx or
             DomainMiny >= DomainMaxy or
             DomainMinz >= DomainMaxz ) :
            print("Error: domain not valid.")
            exit()

except KeyError as ke :
    print("Error reading '", ke.args[0], "' parameter.")
    exit()
except ValueError as vale :
    print("Error parsing numerical parameter:", vale.args[0])
    exit()
except:
    print("Error reading config file ", sys.exc_info()[0])
    exit()

# Ok now we have the data

tree = ET.parse(XmlFile)
root = tree.getroot()

h = float(root.find("./execution/constants/h").get("value"))
mass = float(root.find("./execution/constants/massfluid").get("value"))

TimeStep = float(root.find("./execution/parameters/parameter[@key='TimeOut']").get("value"))

if not CustomDomain:
    #Get domain from xml file
    pmin = root.find("./casedef/geometry/definition/pointmin")
    DomainMinx = float(pmin.get("x"))
    DomainMiny = float(pmin.get("y"))
    DomainMinz = float(pmin.get("z"))
    pmax = root.find("./casedef/geometry/definition/pointmax")
    DomainMaxx = float(pmax.get("x"))
    DomainMaxy = float(pmax.get("y"))
    DomainMaxz = float(pmax.get("z"))

diffuseparticles.run(InputDataPath,
                     InputFilesPrefix,
                     OutputDataPath,
                     FilesOutputPrefix,
                     ExclusionZoneFile,
                     StartingTimeStep,
                     EndingTimeStep,
                     ZeroPadding,
                     TextFiles,
                     PlyFiles,
                     VtkFiles,
                     VtkDiffuseData,
                     VtkFluidData,
                     h, mass, TimeStep,
                     DomainMinx, DomainMiny, DomainMinz, DomainMaxx, DomainMaxy, DomainMaxz,
                     MinTrappedAirThreshold, MaxTrappedAirThreshold,
                     MinWaveCrestsThreshold, MaxWaveCrestsThreshold,
                     MinKineticEnergyThreshold, MaxKineticEnergyThreshold,
                     DiffuseTrappedAirMultiplier, DiffuseWaveCrestsMultiplier,
                     SprayDensity, BubblesDensity, LifefimeMultiplier,
                     BuoyancyControl, DragControl)

