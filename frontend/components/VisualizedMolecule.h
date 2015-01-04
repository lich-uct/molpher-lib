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

#include "auxiliary/QtMocHack.h"

#include <QtCore/QList>
#include <QtCore/QPointer>
#include <QtCore/QProcess>
#include <QtGui/QGraphicsItem>
#include <QtGui/QPainter>

#include "Edge.h"
#include "MolpherMolecule.h"

#define MOLECULE_DIAMETER 15.0

class Tab;

class VisualizedMolecule : public QObject, public QGraphicsItem
{
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

public:
    enum ShapeDefinition {
        SD_DEFAULT,
        SD_NEW_CANDIDATE,
        SD_NEIGHBORHOOD_ORIGIN,
        SD_IDENTITY_SEARCH
    };

    enum ColorDefinition {
        CD_SOURCE,
        CD_TARGET,
        CD_NODE, // inner node in candidate tree
        CD_LEAF, // leaf in candidate tree
        CD_NEIGHBOR,
        CD_DECOY,
        CD_PRUNED
    };

    enum Priority {
        P_PRUNED = 1, // the lowest molecule priority
        P_NEIGHBOR,
        P_NODE,
        P_LEAF,
        P_DECOY,
        P_SOURCE,
        P_TARGET,
        P_PATH_MOLECULE // the highest molecule priority
    };

    VisualizedMolecule(const MolpherMolecule &molpherMolecule,
        const ColorDefinition &color, const ShapeDefinition &shape,
        const bool pathHighlightingByMouseHover);
    ~VisualizedMolecule();

    void SetDefaultZValue(QPointer<VisualizedMolecule> molecule);
    void SetPathMoleculeZValue(QPointer<VisualizedMolecule> molecule,
        const bool pathIsHighlighted);
    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
    void SetPosition(const qreal x, const qreal y);
    void SetColor(const ColorDefinition color);
    VisualizedMolecule::ColorDefinition GetColor() const;
    void SetShape(const ShapeDefinition shape);
    VisualizedMolecule::ShapeDefinition GetShape() const;
    std::string GetSmiles() const;
    std::string GetFormula() const;
    void DeselectNeighborhoodOrigin();
    void AddEdgeToDescendant(const QPointer<Edge> edge);
    void AddEdgeToParent(const QPointer<Edge> edge);
    void TakeAwayEdgeToDescendant(const QPointer<Edge> edge);
    void TakeAwayEdgeToParent(const QPointer<Edge> edge);
    void ScaleEdgeToParent(qreal factor);
    QPointer<Edge> GetEdgeToParent() const;
    void GetEdgesToDescendants(QList<QPointer<Edge> > &edges) const;
    void DeleteEdges();
    bool IsConnected() const;
    void SetPathHighlighting(const bool byMouseHover);
    void SetHighlightingOfPathToSource(const bool isHighlighted);
    bool IsEndOfHighlightedPath();

protected:
    void contextMenuEvent(QGraphicsSceneContextMenuEvent *event);
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void StartMarvinProcess(QString &program);

    void FillColor(QColor &color) const;
    void SetTriangleWithBorder(QPolygonF &polygon) const;
    void SetTriangleWithoutBorder(QPolygonF &polygon) const;
    void SetRectangleWithBorder(QRectF &rectangle) const;
    void SetRectangleWithoutBorder(QRectF &rectangle) const;
    void SetHalfRectangleWithoutBorder(QRectF &rectangle) const;
    void SetRhombusWithoutBorder(QPolygonF &polygon) const;
    void RevertShape();
    bool PathToSourceExists() const;

protected slots:
    void OnMarvinProcessError(QProcess::ProcessError err);

    void OnActionBookmark();
    void OnActionToggleNeighborhoodOrigin();
    void OnSelectMoleculesToRoot();
    void OnActionCopySmiles();
    void OnActionCopyFormula();
    void OnActionMarvinSketch();
    void OnActionMarvinSpace();
    void OnActionMarvinView();
    void OnActionIdentityPubchem();
    void OnActionIdentityZinc();
    void OnActionIdentityChembl();
    void OnActionSimilarityPubchem();
    void OnActionSimilarityZinc();
    void OnActionSimilarityChembl();

signals:
    void ChangeNeighborhoodOrigin(const std::string &smiles);
    void ChangeHighlightingRule(const std::string &smiles);
    void RevertHighlightedPath();

private:
    // variables, which were read from MolpherMolecule class
    std::string mSmiles;
    std::string mFormula;
    std::string mParentSmile;

    ColorDefinition mColor;
    ShapeDefinition mShape;
    ShapeDefinition mLastShape; // when deselect neighborhood origin
    int mPenWidth;  // then remove and check whether QRectF is still OK

    bool mPathHighlightingByMouseHover;  // otherwise by middle button

    QPointer<Edge> mEdgeToParent;
    QList<QPointer<Edge> > mEdgesToDescendants;
};
