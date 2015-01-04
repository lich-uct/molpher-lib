/*
 Copyright (c) 2012 Marek Mikes
 Copyright (c) 2012 Petr Koupy

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

#include <boost/date_time/posix_time/posix_time.hpp>

#include <cassert>

#include <QtCore/QString>
#include <QtCore/QByteArray>
#include <QtCore/QBuffer>
#include <QtCore/QUrl>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtGui/QMenu>
#include <QtGui/QImage>
#include <QtGui/QPixmap>
#include <QtGui/QGraphicsSceneContextMenuEvent>
#include <QtGui/QDesktopServices>
#include <QtGui/QApplication>
#include <QtGui/QClipboard>
#include <QtGui/QPen>

#include "dialog/dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "auxiliary/GlobalObjectsHolder.h"
#include "inout.h"
#include "MoleculePainter.h"
#include "Tab.h"
#include "components/VisualizedMolecule.h"

VisualizedMolecule::VisualizedMolecule(const MolpherMolecule &molpherMolecule,
        const ColorDefinition &color, const ShapeDefinition &shape,
        const bool pathHighlightingByMouseHover) :
    mSmiles(molpherMolecule.smile),
    mFormula(molpherMolecule.formula),
    mParentSmile(molpherMolecule.parentSmile),
    mColor(color),
    mShape(shape),
    mLastShape(shape),
    mPenWidth(0),
    mPathHighlightingByMouseHover(pathHighlightingByMouseHover),
    mEdgeToParent(NULL)
{
    this->setCacheMode(QGraphicsItem::NoCache);
    this->setPos(molpherMolecule.posX, molpherMolecule.posY);
    this->setAcceptsHoverEvents(true);
    this->setCursor(Qt::CrossCursor);
    this->setFlag(QGraphicsItem::ItemIsSelectable);

    SetDefaultZValue(this);
}

VisualizedMolecule::~VisualizedMolecule()
{
    DeleteEdges();
}

void VisualizedMolecule::SetDefaultZValue(QPointer<VisualizedMolecule> molecule)
{
    assert(!molecule.isNull());
    switch (molecule->GetColor()) {
    case CD_SOURCE:
        molecule->setZValue(P_SOURCE);
        break;
    case CD_TARGET:
        molecule->setZValue(P_TARGET);
        break;
    case CD_NODE:
        molecule->setZValue(P_NODE);
        break;
    case CD_LEAF:
        molecule->setZValue(P_LEAF);
        break;
    case CD_NEIGHBOR:
        molecule->setZValue(P_NEIGHBOR);
        break;
    case CD_DECOY:
        molecule->setZValue(P_DECOY);
        break;
    case CD_PRUNED:
        molecule->setZValue(P_PRUNED);
        break;
    default:
        assert(false);
        break;
    }
}

void VisualizedMolecule::SetPathMoleculeZValue(
        QPointer<VisualizedMolecule> molecule, const bool pathIsHighlighted)
{
    assert(!molecule.isNull());
    if (pathIsHighlighted) {
        molecule->setZValue(P_PATH_MOLECULE);
    } else {
        SetDefaultZValue(molecule);
    }
}

QRectF VisualizedMolecule::boundingRect() const
{
    qreal radius = MOLECULE_DIAMETER / 2;
    return QRectF(-radius - mPenWidth / 2, -radius - mPenWidth / 2,
        MOLECULE_DIAMETER + mPenWidth, MOLECULE_DIAMETER + mPenWidth);
}

QPainterPath VisualizedMolecule::shape() const
{
    QPainterPath path;
    QRectF rec;
    QPolygonF polygon;
    switch (mShape) {
    case (SD_DEFAULT):
    case (SD_NEIGHBORHOOD_ORIGIN):
        SetRectangleWithBorder(rec);
        path.addEllipse(rec);
        break;
    case (SD_NEW_CANDIDATE):
        SetTriangleWithBorder(polygon);
        path.addPolygon(polygon);
        break;
    default:
        assert(false);
        break;
    }
    return path;
}

void VisualizedMolecule::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
        QWidget *widget)
{
    QColor color;
    FillColor(color);
    if (this->isSelected()) {
        painter->setBrush(color.darker());
    } else {
        painter->setBrush(color);
    }

    QPen pen;
    pen.setWidth(mPenWidth);
    painter->setPen(pen);

    QRectF rec;
    QPolygonF polygon;
    switch (mShape) {
    case (SD_DEFAULT):
        SetRectangleWithoutBorder(rec);
        painter->drawEllipse(rec);
        break;
    case (SD_NEIGHBORHOOD_ORIGIN):
        SetRectangleWithoutBorder(rec);
        painter->drawEllipse(rec);

        pen.setWidth(0);
        painter->setPen(pen);
        painter->setBrush(Qt::white);
        SetHalfRectangleWithoutBorder(rec);
        painter->drawEllipse(rec);
        break;
    case (SD_NEW_CANDIDATE):
        SetTriangleWithoutBorder(polygon);
        painter->drawPolygon(polygon);
        break;
    case (SD_IDENTITY_SEARCH):
        SetRhombusWithoutBorder(polygon);
        painter->drawPolygon(polygon);
        break;
    default:
        assert(false);
        break;
    }
}

void VisualizedMolecule::SetPosition(const qreal x, const qreal y)
{
    this->setPos(x, y);

    if ((CD_NODE == mColor) || (CD_LEAF == mColor) ||
            ((CD_DECOY == mColor) && !mEdgeToParent.isNull())) {
        assert(!mEdgeToParent.isNull());
        mEdgeToParent->SetTailCoordinates(x, y);
    }
    foreach (QPointer<Edge> edge, mEdgesToDescendants) {
        edge->SetHeadCoordinates(x, y);
    }
}

void VisualizedMolecule::SetColor(const ColorDefinition color)
{
    mColor = color;
    SetDefaultZValue(this);
}

VisualizedMolecule::ColorDefinition VisualizedMolecule::GetColor() const
{
    return mColor;
}

void VisualizedMolecule::SetShape(const ShapeDefinition shape)
{
    mLastShape = mShape;
    mShape = shape;
}

VisualizedMolecule::ShapeDefinition VisualizedMolecule::GetShape() const
{
    return mShape;
}

std::string VisualizedMolecule::GetSmiles() const
{
    return mSmiles;
}

std::string VisualizedMolecule::GetFormula() const
{
    return mFormula;
}

void VisualizedMolecule::DeselectNeighborhoodOrigin()
{
    if (SD_NEIGHBORHOOD_ORIGIN != mShape) {
        return;
    }

    RevertShape();
    update(boundingRect());
}

void VisualizedMolecule::AddEdgeToDescendant(const QPointer<Edge> edge)
{
    mEdgesToDescendants.append(edge);
}

void VisualizedMolecule::AddEdgeToParent(const QPointer<Edge> edge)
{
    assert(!mEdgeToParent);
    mEdgeToParent = edge;
}

void VisualizedMolecule::TakeAwayEdgeToDescendant(const QPointer<Edge> edge)
{
    mEdgesToDescendants.removeOne(edge);
}

void VisualizedMolecule::TakeAwayEdgeToParent(const QPointer<Edge> edge)
{
    assert(edge == mEdgeToParent);
    mEdgeToParent = NULL;
}

void VisualizedMolecule::ScaleEdgeToParent(qreal factor)
{
    if ((CD_NODE == mColor) || (CD_LEAF == mColor)) {
        assert(!mEdgeToParent.isNull());
        mEdgeToParent->Scale(factor);
    } else if (!mEdgeToParent.isNull() &&
            ((CD_DECOY == mColor) || (CD_TARGET == mColor))) {
        mEdgeToParent->Scale(factor);
    }
}

QPointer<Edge> VisualizedMolecule::GetEdgeToParent() const
{
    return mEdgeToParent;
}

void VisualizedMolecule::GetEdgesToDescendants(QList<QPointer<Edge> > &edges) const
{
    edges = mEdgesToDescendants;
}

void VisualizedMolecule::DeleteEdges()
{
    if (!mEdgeToParent.isNull()) {
        delete mEdgeToParent;
    }

    while (!mEdgesToDescendants.isEmpty()) {
        delete mEdgesToDescendants.takeFirst();
    }
}

bool VisualizedMolecule::IsConnected() const
{
    return !mEdgeToParent.isNull() || !mEdgesToDescendants.isEmpty();
}

void VisualizedMolecule::SetPathHighlighting(const bool byMouseHover)
{
    mPathHighlightingByMouseHover = byMouseHover;
}

void VisualizedMolecule::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    if (mPathHighlightingByMouseHover && PathToSourceExists()) {
        SetHighlightingOfPathToSource(true);
    }

    std::string filename;
    filename += QDir::currentPath().toStdString();
    filename += "/Temp/";
    filename += gFrontendId + "/";
    filename += boost::posix_time::to_iso_string(
        boost::posix_time::second_clock::local_time());
    filename += ".png";

    QString program("indigo-depict");

    QStringList arguments;
    arguments << "-";
    arguments << QString::fromStdString(mSmiles);
    arguments << QString::fromStdString(filename);
    arguments << "-bond" << "25";
    arguments << "-comment" << QString::fromStdString(mFormula);
    arguments << "-commentoffset" << "8";
    arguments << "-commentsize" << "14";
    arguments << "-commentalign" << "0";

    QImage image;
    if (QProcess::execute(program, arguments) == 0) {
        image = QImage(QString::fromStdString(filename));
        QFile::remove(QString::fromStdString(filename));
    } else {
        image = MoleculePainter::DrawMolecule(mSmiles).toImage();
    }

    QByteArray data;
    QBuffer buffer(&data);
    buffer.open(QIODevice::WriteOnly);
    image.save(&buffer, "PNG");
    QString tooltip = QString("<img src=\"data:image/png;base64,%1\">").arg(
        QString(buffer.data().toBase64()));
    setToolTip("<html>" + tooltip + "</html>");
}

void VisualizedMolecule::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    if (mPathHighlightingByMouseHover && PathToSourceExists()) {
        SetHighlightingOfPathToSource(false);
    }

    setToolTip("");
}

void VisualizedMolecule::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if (event->button() == Qt::MiddleButton) {
        if (mPathHighlightingByMouseHover) {
            if (PathToSourceExists()) {
                emit ChangeHighlightingRule(mSmiles);
            }
        } else {
            if (IsEndOfHighlightedPath()) {
                emit ChangeHighlightingRule("");
            } else {
                if (PathToSourceExists()) {
                    emit RevertHighlightedPath();
                    emit ChangeHighlightingRule(mSmiles);
                    SetHighlightingOfPathToSource(true);
                }
            }
        }
    } else {
        QGraphicsItem::mousePressEvent(event);
    }
}

void VisualizedMolecule::FillColor(QColor &color) const
{
    // TODO load color (depending on molecule type) from global color store
    switch (mColor) {
    case CD_SOURCE:
        color = Qt::magenta;
        break;
    case CD_TARGET:
        color = Qt::red;
        break;
    case CD_NODE:
        color = Qt::cyan;
        break;
    case CD_LEAF:
        color = Qt::green;
        break;
    case CD_NEIGHBOR:
        color = Qt::yellow;
        break;
    case CD_DECOY:
        color = QColor(255, 165, 0); // orange
        break;
    case CD_PRUNED:
        color = Qt::lightGray;
        break;
    default:
        assert(false);
        break;
    }
    assert(color.isValid());
}

void VisualizedMolecule::SetTriangleWithBorder(QPolygonF &polygon) const
{
    polygon.clear();

    qreal radius = MOLECULE_DIAMETER / 2;
    polygon.push_back(QPointF(-radius - mPenWidth / 2, radius + mPenWidth / 2));
    polygon.push_back(QPointF(radius + mPenWidth / 2, radius + mPenWidth / 2));
    polygon.push_back(QPointF(0, -radius - mPenWidth / 2));
}

void VisualizedMolecule::SetTriangleWithoutBorder(QPolygonF &polygon) const
{
    polygon.clear();

    qreal radius = MOLECULE_DIAMETER / 2;
    polygon.push_back(QPointF(-radius, radius));
    polygon.push_back(QPointF(radius, radius));
    polygon.push_back(QPointF(0, -radius));
}

void VisualizedMolecule::SetRectangleWithBorder(QRectF &rectangle) const
{
    qreal radius = MOLECULE_DIAMETER / 2;
    rectangle.setRect(-radius - mPenWidth / 2, -radius - mPenWidth / 2,
        MOLECULE_DIAMETER + mPenWidth, MOLECULE_DIAMETER + mPenWidth);
}

void VisualizedMolecule::SetRectangleWithoutBorder(QRectF &rectangle) const
{
    qreal radius = MOLECULE_DIAMETER / 2;
    rectangle.setRect(-radius, -radius, MOLECULE_DIAMETER, MOLECULE_DIAMETER);
}

void VisualizedMolecule::SetHalfRectangleWithoutBorder(QRectF &rectangle) const
{
    qreal radius = MOLECULE_DIAMETER / 2;
    qreal oneFourthOfDiameter = MOLECULE_DIAMETER / 4;
    rectangle.setRect(-oneFourthOfDiameter, -oneFourthOfDiameter, radius, radius);
}

void VisualizedMolecule::SetRhombusWithoutBorder(QPolygonF &polygon) const
{
    polygon.clear();

    qreal radius = MOLECULE_DIAMETER / 2;
    polygon.push_back(QPointF(-radius, 0));
    polygon.push_back(QPointF(0, radius));
    polygon.push_back(QPointF(radius, 0));
    polygon.push_back(QPointF(0, -radius));
}

void VisualizedMolecule::RevertShape()
{
    SetShape(mLastShape);
}

void VisualizedMolecule::SetHighlightingOfPathToSource(const bool isHighlighted)
{
    QPointer<Edge> edge = mEdgeToParent;
    QPointer<VisualizedMolecule> mol = this;
    SetPathMoleculeZValue(mol, isHighlighted);
    while (!edge.isNull()) {
        edge->SetHighlighting(isHighlighted);
        mol = edge->GetHead();
        assert(mol);
        SetPathMoleculeZValue(mol, isHighlighted);
        edge = mol->GetEdgeToParent();
    }
}

bool VisualizedMolecule::IsEndOfHighlightedPath()
{
    bool parentEdgeIsHighlighted = false;
    if (!mEdgeToParent.isNull()) {
        parentEdgeIsHighlighted = mEdgeToParent->IsHighlighted();
    }
    bool childEdgeIsHighlighted = false;
    foreach (QPointer<Edge> edge, mEdgesToDescendants) {
        if (edge->IsHighlighted()) {
            childEdgeIsHighlighted = true;
            break;
        }
    }

    return parentEdgeIsHighlighted && !childEdgeIsHighlighted;
}

bool VisualizedMolecule::PathToSourceExists() const
{
    return (CD_NODE == mColor) || (CD_LEAF == mColor) || (((CD_TARGET == mColor) ||
            (CD_DECOY == mColor)) && !mEdgeToParent.isNull());
}

void VisualizedMolecule::OnActionBookmark()
{
    Bookmarks *b = gGlobalObjectsHolder.GetBookmarks();
    MolpherMolecule self(mSmiles, mFormula);
    if (!b->Add(self)) {
        ShowWarning("There is no bookmark group.");
    }
}

void VisualizedMolecule::OnActionToggleNeighborhoodOrigin()
{
    if (SD_NEIGHBORHOOD_ORIGIN == mShape) {
        std::string empty;
        emit ChangeNeighborhoodOrigin(empty);
        return;
    }

    SetShape(SD_NEIGHBORHOOD_ORIGIN);
    setSelected(false);
    update(boundingRect());

    emit ChangeNeighborhoodOrigin(mSmiles);
}

void VisualizedMolecule::OnSelectMoleculesToRoot()
{
    QPointer<Edge> edge = mEdgeToParent;
    QPointer<VisualizedMolecule> mol = this;
    while (!edge.isNull()) {
        mol->setSelected(true);
        mol = edge->GetHead();
        assert(mol);
        edge = mol->GetEdgeToParent();
    }
    mol->setSelected(true);
}

void VisualizedMolecule::OnActionCopySmiles()
{
    QApplication::clipboard()->setText(QString::fromStdString(mSmiles));
}

void VisualizedMolecule::OnActionCopyFormula()
{
    QApplication::clipboard()->setText(QString::fromStdString(mFormula));
}

void VisualizedMolecule::StartMarvinProcess(QString &program)
{
    std::string filename;
    filename += QDir::currentPath().toStdString();
    filename += "/Temp/";
    filename += gFrontendId + "/";
    filename += boost::posix_time::to_iso_string(
        boost::posix_time::second_clock::local_time());
    filename += ".sdf";

    MolpherMolecule mol(mSmiles);
    std::vector<MolpherMolecule> mols;
    mols.push_back(mol);
    WriteMolphMolsToSDF(filename, mols);

    QStringList arguments;
    arguments << QString::fromStdString(filename);
    QProcess *process = new QProcess(this);
    QObject::connect(process, SIGNAL(error(QProcess::ProcessError)),
        this, SLOT(OnMarvinProcessError(QProcess::ProcessError)));
    process->start(program, arguments);
}

void VisualizedMolecule::OnMarvinProcessError(QProcess::ProcessError err)
{
    ShowInformation(tr(
        "To use this feature, please install MarvinBeans package from"
        " http://www.chemaxon.com/download/marvin/ and add its"
        " ChemAxon/MarvinBeans directory to the system Path variable."
        ));
}

void VisualizedMolecule::OnActionMarvinSketch()
{
    QString program("MarvinSketch");
    StartMarvinProcess(program);
}

void VisualizedMolecule::OnActionMarvinSpace()
{
    QString program("MarvinSpace");
    StartMarvinProcess(program);
}

void VisualizedMolecule::OnActionMarvinView()
{
    QString program("MarvinView");
    StartMarvinProcess(program);
}

void VisualizedMolecule::OnActionIdentityPubchem()
{
    // http://pubchem.ncbi.nlm.nih.gov/search/help_search.html#UrlSch
    std::string url;
    url.append("http://pubchem.ncbi.nlm.nih.gov/search/search.cgi?cmd=search&q_type=dt&q_data=");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    url.append("&simp_schtp=fs");
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::OnActionIdentityZinc()
{
    // http://wiki.bkslab.org/index.php/ZINCCL
    std::string url;
    url.append("http://zinc.docking.org/results?structure.smiles=");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    url.append("&structure.similarity=1.0");
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::OnActionIdentityChembl()
{
    // https://www.ebi.ac.uk/chembldb/ws
    std::string url;
    url.append("https://www.ebi.ac.uk/chemblws/compounds/smiles/");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::OnActionSimilarityPubchem()
{
    // http://pubchem.ncbi.nlm.nih.gov/search/help_search.html#UrlSch
    std::string url;
    url.append("http://pubchem.ncbi.nlm.nih.gov/search/search.cgi?cmd=search&q_type=dt&q_data=");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    url.append("&simp_schtp=80");
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::OnActionSimilarityZinc()
{
    // http://wiki.bkslab.org/index.php/ZINCCL
    std::string url;
    url.append("http://zinc.docking.org/results?structure.smiles=");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    url.append("&structure.similarity=0.8");
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::OnActionSimilarityChembl()
{
    // https://www.ebi.ac.uk/chembldb/ws
    std::string url;
    url.append("https://www.ebi.ac.uk/chemblws/compounds/similarity/");
    url.append(QByteArray(mSmiles.c_str()).toPercentEncoding().constData());
    url.append("/80");
    QDesktopServices::openUrl(QUrl::fromEncoded(QByteArray(url.c_str()), QUrl::StrictMode));
}

void VisualizedMolecule::contextMenuEvent(QGraphicsSceneContextMenuEvent *event)
{
    QMenu mainMenu;

    QAction *action = new QAction(tr("Create molecule bookmark"), &mainMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionBookmark()));
    mainMenu.addAction(action);

    action = new QAction(tr("Toggle neighborhood origin"), &mainMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionToggleNeighborhoodOrigin()));
    mainMenu.addAction(action);

    action = new QAction(tr("Select molecules to source"), &mainMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnSelectMoleculesToRoot()));
    mainMenu.addAction(action);
    if ((CD_NODE != mColor) && (CD_LEAF != mColor) &&
            ((CD_TARGET != mColor) || mEdgeToParent.isNull())) {
        action->setEnabled(false);
    }

    QMenu *copyMenu = mainMenu.addMenu(tr("Copy molecule as"));

    action = new QAction(tr("SMILES"), copyMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionCopySmiles()));
    copyMenu->addAction(action);

    action = new QAction(tr("Formula"), copyMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionCopyFormula()));
    copyMenu->addAction(action);

    QMenu *openMenu = mainMenu.addMenu(tr("Open externally in"));

    action = new QAction(tr("Marvin Sketch"), openMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionMarvinSketch()));
    openMenu->addAction(action);

    action = new QAction(tr("Marvin Space"), openMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionMarvinSpace()));
    openMenu->addAction(action);

    action = new QAction(tr("Marvin View"), openMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionMarvinView()));
    openMenu->addAction(action);

    QMenu *identityMenu = mainMenu.addMenu(tr("Find identical molecule in"));

    action = new QAction(tr("Pubchem"), identityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionIdentityPubchem()));
    identityMenu->addAction(action);

    action = new QAction(tr("ZINC"), identityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionIdentityZinc()));
    identityMenu->addAction(action);

    action = new QAction(tr("ChEMBL"), identityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionIdentityChembl()));
    identityMenu->addAction(action);

    QMenu *similarityMenu = mainMenu.addMenu(tr("Find similar molecules in"));

    action = new QAction(tr("Pubchem"), similarityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionSimilarityPubchem()));
    similarityMenu->addAction(action);

    action = new QAction(tr("ZINC"), similarityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionSimilarityZinc()));
    similarityMenu->addAction(action);

    action = new QAction(tr("ChEMBL"), similarityMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnActionSimilarityChembl()));
    similarityMenu->addAction(action);

    mainMenu.exec(event->screenPos());
}
