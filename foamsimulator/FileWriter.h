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

#ifndef FILEWRITER_H
#define FILEWRITER_H

/**
   This abstract class defines the method that a class that write output files needs to implement.
 */
class FileWriter {
 public:
  /**
     Dumps the data to a file.
     \return Returns 0 if the writing fails. Any other value otherwise.
   */
  virtual int write() =0;
};

#endif
