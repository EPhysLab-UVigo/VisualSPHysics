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

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include "Ops.h"
#include "VtkDWriter.h"

VtkDWriter::VtkDWriter(std::string const &name, double xmin, double xmax,
                       double ymin, double ymax, double zmin, double zmax,
                       double h)
    : bc(xmin, xmax, ymin, ymax, zmin, zmax, h), h(h), fileName(name) {

  points = vtkSmartPointer<vtkPoints>::New();

  velocity = vtkSmartPointer<vtkDoubleArray>::New();
  velocity->SetName("Velocity");
  velocity->SetNumberOfComponents(3);

  size = vtkSmartPointer<vtkDoubleArray>::New();
  size->SetName("Size");
}

void VtkDWriter::setData(std::vector<std::array<double, 3>> *d,
                         std::vector<std::array<double, 3>> *v) {
  for (long i = 0; i < d->size(); i++) {
    auto p = d->at(i);
    if (p[0] != 0 && p[1] != 0 && p[2] != 0) {
      oparticle tp;
      tp.pos = p;
      tp.vel = v->at(i);
      bc.addElement(tp, p[0], p[1], p[2]);
    }
  }

  for (auto &bucket : bc.getBuckets()) {
    for (long i = 0; i < bucket.size(); i++) {
      auto newp = bucket[i];
      int psize = 1;
      for (long j = 0; j < bucket.size(); j++) {
        if (j != i && (ops::magnitude(ops::substract(
                           bucket[i].pos, bucket[j].pos)) < (h / 5))) {
          newp.pos[0] = (newp.pos[0] + bucket[j].pos[0]) / 2.;
          newp.pos[1] = (newp.pos[1] + bucket[j].pos[1]) / 2.;
          newp.pos[2] = (newp.pos[2] + bucket[j].pos[2]) / 2.;

          newp.vel[0] = (newp.vel[0] + bucket[j].vel[0]) / 2.;
          newp.vel[1] = (newp.vel[1] + bucket[j].vel[1]) / 2.;
          newp.vel[2] = (newp.vel[2] + bucket[j].vel[2]) / 2.;

          bucket.erase(bucket.begin() + j);
          j--;
          psize++;
        }
      }
      bucket.erase(bucket.begin() + i);
      i--;
      points->InsertNextPoint(newp.pos.data());
      velocity->InsertNextTuple(newp.vel.data());
      size->InsertNextValue(h / (10 / std::pow(psize, 1 / 3)));
    }
  }
}

int VtkDWriter::write() {
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  for (long i = 0; i < points->GetNumberOfPoints(); ++i){
    vtkIdType pt[] = {i};
    vertices->InsertNextCell(1, pt);
  }
  polydata->SetVerts(vertices);
  
  polydata->GetPointData()->AddArray(velocity);
  polydata->GetPointData()->AddArray(size);

  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileTypeToBinary();
  writer->SetInputData(polydata);
  writer->SetFileName(fileName.c_str());
  writer->Write();

  return 1;
}
