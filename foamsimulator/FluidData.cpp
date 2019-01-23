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

#include "FluidData.h"

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
//#include <vtkSelectEnclosedPoints.h>

#include <array>

FluidData::FluidData(double xmin, double xmax, double ymin, double ymax,
                     double zmin, double zmax, double h)
    : bc(xmin, xmax, ymin, ymax, zmin, zmax, h), exclude(false) {

  std::cout << "Number of buckets: " << bc.getBuckets().size() << std::endl;
}

BucketContainer<particle> *FluidData::getBucketContainer() { return &bc; }

void FluidData::setExclusionZone(std::string const &fileName) {
  exclude = true;
  exFile = fileName;
}

bool FluidData::loadFile(std::string const &fileName) {
  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fileName.c_str());
  reader->Update();

  vtkPolyData *output = reader->GetOutput();

  vtkDataArray *points = output->GetPoints()->GetData();

  vtkPointData *pointData = output->GetPointData();

  vtkDataArray *pvel = pointData->GetArray("Vel");
  vtkDataArray *rhop = pointData->GetArray("Rhop"); // Density

  // vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints;

  // if(exclude){
  //   vtkSmartPointer<vtkPolyDataReader> reader =
  //   vtkSmartPointer<vtkPolyDataReader>::New();
  //   reader->SetFileName(exFile.c_str());
  //   reader->Update();
  //   vtkPolyData * exclusion = reader->GetOutput();

  //   //Points inside test
  //   selectEnclosedPoints =  vtkSmartPointer<vtkSelectEnclosedPoints>::New();
  //   selectEnclosedPoints->Initialize(exclusion);
  // }

  for (unsigned int i = 0; i < output->GetPoints()->GetNumberOfPoints(); i++) {
    double *p = points->GetTuple(i), *v = pvel->GetTuple(i);

    // if (exclude && selectEnclosedPoints->IsInsideSurface(p)){
    //   //nPoints--;
    //   continue;
    // }

    particle pi;
    pi.pos = {p[0], p[1], p[2]};
    pi.vel = {v[0], v[1], v[2]};
    pi.rhop = rhop->GetTuple(i)[0];
    bc.addElement(pi, p[0], p[1], p[2]);
  }

  // Remove doubles
  for (auto &bucket : bc.getBuckets()) {
    for (unsigned int i = 0; i < bucket.size(); i++) {
      for (unsigned int j = 0; j < bucket.size(); j++) {
        if (j != i) {
          if (bucket[i].pos[0] == bucket[j].pos[0] &&
              bucket[i].pos[1] == bucket[j].pos[1] &&
              bucket[i].pos[2] == bucket[j].pos[2]) {
            bucket.erase(bucket.begin() + i);
            i--;
            break;
          }
        }
      }
    }
  }

  // Assign ids
  unsigned int idp = 0;
  for (auto &bucket : bc.getBuckets()) { // Iterate over all buckets
    for (auto &pi : bucket) { // Iterate over each particle in the bucket
      pi.id = idp;
      idp++;
    }
  }

  return true;
}
