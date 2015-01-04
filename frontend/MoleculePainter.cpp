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

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Demos/RDKit/Draw/MolDrawing.h>

#include <QtGui/QPen>
#include <QtCore/QLine>
#include <QtGui/QFontMetrics>
#include <QtGui/QPainter>
#include <QtGui/QPixmap>

#include "MoleculePainter.h"

void MoleculePainter::SetColor(int atNum, QPainter *painter)
{
    QColor color;
    switch (atNum) {
    case 7:
        color.setRgbF(0.0, 0.0, 1.0);
        break;
    case 8:
        color.setRgbF(1.0, 0.0, 0.0);
        break;
    case 9:
        color.setRgbF(0.2, 0.8, 0.8);
        break;
    case 15:
        color.setRgbF(1.0, 0.5, 0.0);
        break;
    case 16:
        color.setRgbF(0.8, 0.8, 0.0);
        break;
    case 17:
        color.setRgbF(0.0, 0.8, 0.0);
        break;
    case 35:
        color.setRgbF(0.5, 0.3, 0.1);
        break;
    case 0:
        color.setRgbF(0.5, 0.5, 0.5);
        break;
    default:
        color.setRgbF(0.5, 0.5, 0.5);
        break;
    }
    painter->setPen(color);
}

void MoleculePainter::DrawLine(
    std::vector<int>::const_iterator &pos, QPainter *painter)
{
    int width = *pos++;
    QPen pen(painter->pen());
    pen.setWidth(width * 10);
    painter->setPen(pen);

    //int dashed = *pos++; // Dashed variable never used.
    pos++; // Move the pointer without assigning a value.

    int an1 = *pos++;
    int an2 = *pos++;
    int xp1 = *pos++;
    int yp1 = *pos++;
    int xp2 = *pos++;
    int yp2 = *pos++;
    if (an1 == an2) {
        SetColor(an1, painter);
        QLine line(xp1, yp1, xp2, yp2);
        painter->drawLine(line);
    } else {
        int mx = xp1 + (xp2 - xp1) / 2;
        int my = yp1 + (yp2 - yp1) / 2;

        SetColor(an1, painter);
        QLine line(xp1, yp1, mx, my);
        painter->drawLine(line);
        SetColor(an2, painter);
        line.setLine(mx, my, xp2, yp2);
        painter->drawLine(line);
    }
}

void MoleculePainter::DrawAtom(
    std::vector<int>::const_iterator &pos, QPainter *painter)
{
    int atNum = *pos++;
    double xp = static_cast<double>(*pos++);
    double yp = static_cast<double>(*pos++);
    int slen = *pos++;
    std::string label = "";
    for (int i = 0; i < slen; ++i) {
        label += (char) *pos++;
    }
    RDKit::Drawing::OrientType orient =
        static_cast<RDKit::Drawing::OrientType>(*pos++);
    QFontMetrics fm = painter->fontMetrics();
    QString qLabel(label.c_str());

    double twidth = fm.width(qLabel);
    double theight = fm.height();

    switch (orient) {
    case RDKit::Drawing::W:
        xp -= twidth;
        yp += theight / 2;
        break;
    case RDKit::Drawing::E:
        yp += theight / 2;
        break;
    case RDKit::Drawing::S:
        xp -= twidth / 2;
        yp += theight;
        break;
    case RDKit::Drawing::N:
        xp -= twidth / 2;
        yp -= theight / 2;
        break;
    default:
        xp -= twidth / 2;
        yp += theight / 2;
    }

    QColor color;
    color.setRgbF(1.0, 1.0, 1.0);
    painter->fillRect(
        xp - 10, yp - theight - 10, twidth + 20, theight + 20, color);
    SetColor(atNum, painter);
    painter->drawText(xp, yp, qLabel);
}

void MoleculePainter::MolToQt(const RDKit::ROMol &mol, QPainter *painter,
    int width, int height, int fontSize, int maxDotsPerAngstrom)
{
    RDKit::RWMol cp(mol);
    RDKit::MolOps::Kekulize(cp);
    RDDepict::compute2DCoords(cp);
    std::vector<int> drawing = RDKit::Drawing::DrawMol(cp);

    QColor color;
    color.setRgbF(1.0, 1.0, 1.0);
    painter->setPen(color);
    painter->drawRect(0, 0, width, height);

    //cairo_select_font_face (cr, "sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);

    std::vector<int>::const_iterator pos = drawing.begin();
    int resolution = 0;
    if (*pos != RDKit::Drawing::RESOLUTION) {
        return;
    }

    resolution = *(pos + 1);

    if (!resolution > 0) {
        return;
    }

    pos += 2;
    if (*pos != RDKit::Drawing::BOUNDS) {
        return;
    }

    pos += 3;
    int dwidth;
    int dheight;
    dwidth = *pos++;
    dheight = *pos++;
    if (!(dwidth > 0 && dheight > 0)) {
        return;
    }

    // size of the image in angstroms:
    double xSize;
    double ySize;
    xSize = static_cast<double>(dwidth) / resolution;
    ySize = static_cast<double>(dheight) / resolution;
    double scale = 1.0;

    if (dwidth > width || dheight > height) {
        if ((static_cast<double>(dwidth) / width) >
                (static_cast<double>(dheight) / height)) {
            scale = static_cast<double>(width) / dwidth;
        } else {
            scale = static_cast<double>(height) / dheight;
        }
    } else {
        if (width / xSize > height / ySize) {
            if (width / xSize > maxDotsPerAngstrom) {
                scale = maxDotsPerAngstrom * xSize / width;
            }
        } else {
            if (height / ySize > maxDotsPerAngstrom) {
                scale = maxDotsPerAngstrom * ySize / height;
            }
        }
    }

    scale *= 0.80;
    painter->translate(0.5 * (width - dwidth * scale),
        0.5 * (height - dheight * scale));
    painter->scale(scale, scale);

    // scaling factors here are empirically determined
    double dFontSize = 1.5 * maxDotsPerAngstrom * fontSize / 14;
    QFont font = painter->font();
    font.setPointSizeF(dFontSize);
    painter->setFont(font);

    while (pos != drawing.end()) {
        int token = *pos++;
        switch (token) {
        case RDKit::Drawing::LINE:
            DrawLine(pos, painter);
            break;
        case RDKit::Drawing::ATOM:
            DrawAtom(pos, painter);
            break;
        default:
            break;
        }
    }
}

QPixmap MoleculePainter::DrawMolecule(const std::string &smile)
{
    QPixmap pixmap(200, 200);
    RDKit::RWMol *mol = RDKit::SmilesToMol(smile);
    if (mol) {
        pixmap.fill(Qt::white);
        QPainter painter(&pixmap);
        MolToQt(*mol, &painter, 200, 200);
        delete mol;
    }
    return pixmap;
}
