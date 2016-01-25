/*
 Copyright (c) 2013 Milan Voršilák

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <map>

#include <rdkit/GraphMol/GraphMol.h>

class SAScore {
public:
    static SAScore* getInstance();
    static SAScore* destroyInstance();
    virtual ~SAScore();
    double getScore(RDKit::ROMol &mol);
    static void loadData();
    static void loadData(const std::string& path);

private:

    SAScore();
    // Don't forget to declare these two. You want to make sure they
    // are unaccessable otherwise you may accidently get copies of
    // your singleton appearing.
    SAScore(const SAScore& orig);  // don't implement
    void operator=(SAScore const&);  // don't implement
    std::map<unsigned int, double> data;

    static SAScore * instance;
};
