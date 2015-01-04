/*
 Copyright (c) 2012 Vladimir Fiklik

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

#include <vector>

#include <GraphMol/ROMol.h>

class QPainter;
class QPixmap;

class MoleculePainter
{
public:
    static QPixmap DrawMolecule(const std::string &smile);

private:
    static void DrawAtom(std::vector<int>::const_iterator &pos, QPainter *painter);
    static void DrawLine(std::vector<int>::const_iterator &pos, QPainter *painter);
    static void MolToQt(const RDKit::ROMol &mol, QPainter *painter,
        int width, int height, int fontSize = 14, int maxDotsPerAngstrom = 30);
    static void SetColor(int atNum, QPainter *painter);
};
