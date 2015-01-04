/*
 Copyright (c) 2012 Marek Mikes

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

#include <QtCore/QPoint>
#include <QtCore/QPointer>
#include <QtGui/QGraphicsItem>
#include <QtGui/QGraphicsSimpleTextItem>

#include "chemoper_selectors.h"

#define EDGE_DEFAULT_WIDTH 0.5
#define EDGE_HIGHLIGHT_WIDTH 5

class VisualizedMolecule;

class Edge: public QObject, public QGraphicsItem
{
public:
    Edge(VisualizedMolecule *head, VisualizedMolecule *tail,
        const ChemOperSelector chemOper, qreal scaleFactor);
    ~Edge();

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
    void SetHeadCoordinates(const qreal x, const qreal y);
    void SetTailCoordinates(const qreal x, const qreal y);
    void SetHighlighting(const bool isHighlighted);
    bool IsHighlighted() const;
    void Scale(qreal factor);
    QPointer<VisualizedMolecule> GetHead();
    QGraphicsSimpleTextItem *GetText() const;
    void RotateText(qreal angle);

protected:
    void NotifyHeadOfDeleting();
    void NotifyTailOfDeleting();
    void SetTextPosition();
    void ScaleText(qreal factor) const;

private:
    QPointF mHeadPos;
    QPointF mTailPos;
    QPointer<VisualizedMolecule> mHead;
    QPointer<VisualizedMolecule> mTail;
    ChemOperSelector mChemOper;
    bool mIsHighlighted;
    qreal mScaleFactor;

    QGraphicsSimpleTextItem *mText;
    qreal mLastRotateAngle;
};
