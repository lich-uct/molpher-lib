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

#include <cassert>

#include <QtCore/QPointer>
#include <QtGui/QFont>
#include <QtGui/QPainter>
#include <QtGui/QStyleOptionGraphicsItem>
#include <QtGui/QWidget>

#include "VisualizedMolecule.h"
#include "Edge.h"

Edge::Edge(VisualizedMolecule *head, VisualizedMolecule *tail,
        const ChemOperSelector chemOper, qreal scaleFactor) :
    mHeadPos(QPointF(head->x(), head->y())),
    mTailPos(QPointF(tail->x(), tail->y())),
    mHead(head),
    mTail(tail),
    mChemOper(chemOper),
    mIsHighlighted(false),
    mScaleFactor(scaleFactor),
    mLastRotateAngle(0.0)
{
    setCacheMode(QGraphicsItem::NoCache);
    setZValue(FLT_MIN);

    assert(mHead);
    assert(mTail);
    mHead->AddEdgeToDescendant(QPointer<Edge>(this));
    mTail->AddEdgeToParent(QPointer<Edge>(this));

    mText = new QGraphicsSimpleTextItem(ChemOperShortDesc(mChemOper));
    QFont font;
    font.setPixelSize(15);
    font.setBold(true);
    mText->setFont(font);
    mText->setZValue(FLT_MAX);
    mText->setVisible(false);
    SetTextPosition();
    ScaleText(scaleFactor);
}

Edge::~Edge()
{
    NotifyHeadOfDeleting();
    NotifyTailOfDeleting();

    delete mText;
}

QRectF Edge::boundingRect() const
{
    qreal left = 0;
    qreal top = 0;
    qreal width = 0;
    qreal height = 0;
    if (mHeadPos.x() <= mTailPos.x()) {
        width = mTailPos.x() - mHeadPos.x();
        left = mHeadPos.x();
    } else {
        width = mHeadPos.x() - mTailPos.x();
        left = mTailPos.x();
    }
    if (mHeadPos.y() <= mTailPos.y()) {
        height = mTailPos.y() - mHeadPos.y();
        top = mHeadPos.y();
    } else {
        height = mHeadPos.y() - mTailPos.y();
        top = mTailPos.y();
    }

    return QRectF(left, top, width, height);
}

void Edge::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
    QWidget *widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    QLineF edge(mHeadPos, mTailPos);

    if (mIsHighlighted) {
        painter->setPen(QPen(QColor(255, 206, 9), EDGE_HIGHLIGHT_WIDTH * mScaleFactor,
            Qt::SolidLine, Qt::FlatCap, Qt::BevelJoin));
    } else {
        painter->setPen(QPen(Qt::lightGray, EDGE_DEFAULT_WIDTH * mScaleFactor,
            Qt::SolidLine, Qt::FlatCap, Qt::BevelJoin));
    }
    painter->drawLine(edge);
}

void Edge::SetHeadCoordinates(const qreal x, const qreal y)
{
    mHeadPos.setX(x);
    mHeadPos.setY(y);

    SetTextPosition();
    update(boundingRect());
}

void Edge::SetTailCoordinates(const qreal x, const qreal y)
{
    mTailPos.setX(x);
    mTailPos.setY(y);

    SetTextPosition();
    update(boundingRect());
}

void Edge::SetTextPosition()
{
    mText->setPos((mHeadPos + mTailPos) / 2);
}

void Edge::SetHighlighting(const bool isHighlighted)
{
    if (isHighlighted == mIsHighlighted) {
        return;
    }

    mIsHighlighted = isHighlighted;
    mText->setVisible(mIsHighlighted);
    update(boundingRect());
}

bool Edge::IsHighlighted() const
{
    return mIsHighlighted;
}

void Edge::Scale(qreal factor)
{
    mScaleFactor *= factor;
    ScaleText(factor);
}

QPointer<VisualizedMolecule> Edge::GetHead()
{
    return mHead;
}

QGraphicsSimpleTextItem *Edge::GetText() const
{
    return mText;
}

void Edge::RotateText(qreal angle)
{
    mText->rotate(-mLastRotateAngle + angle);
    mLastRotateAngle = angle;
}

void Edge::NotifyHeadOfDeleting()
{
    mHead->TakeAwayEdgeToDescendant(QPointer<Edge>(this));
}

void Edge::NotifyTailOfDeleting()
{
    mTail->TakeAwayEdgeToDescendant(QPointer<Edge>(this));
}

void Edge::ScaleText(qreal factor) const
{
    mText->scale(factor, factor);
}
