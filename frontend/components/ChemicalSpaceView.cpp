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
#include <cfloat>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>

#include <QtCore/QList>
#include <QtCore/QPointer>
#include <QtCore/QPointF>

#include "ChemicalSpaceView.h"
#include "tests/DbTestingManager.h"
#include "dialog/dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "dialog/NeighborhoodDialog.h"
#include "dialog/RevisualizeDialog.h"
#include "auxiliary/GlobalObjectsHolder.h"

ChemicalSpaceView::ChemicalSpaceView(IterationSnapshot &snapshot, bool pruneMeansDelete) :
    mNeighborhoodOriginSmiles(""),
    // set this option by input parameter or from some global object in future
    mPruneMeansDelete(pruneMeansDelete),
    mLastRotateAngleWithoutCorrection(0),
    mCurrentScale(1),
    mEndMoleculeOfHighlightedPath("")
{
    this->setDragMode(QGraphicsView::RubberBandDrag);
    this->setCacheMode(QGraphicsView::CacheNone);
    this->setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
    this->setRenderHint(QPainter::Antialiasing);
    this->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    this->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
    // scene's extent is defined because of correct view centring
    this->setSceneRect(-1500.0, -1500.0, 3000.0, 3000.0);

    mChemicalSpaceScene = new QGraphicsScene(this);
    // if NoIndex is not set, sometimes program crashes after removing item
    //  because of segmentation fault
    mChemicalSpaceScene->setItemIndexMethod(QGraphicsScene::NoIndex);
    this->setScene(mChemicalSpaceScene);

    SetNewSnapshot(snapshot, VisualizedMolecule::SD_DEFAULT);

    qRegisterMetaType<NeighborhoodTaskResult>();
    connect(&gCommunicator,
        SIGNAL(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)),
        this,
        SLOT(OnVisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)),
        Qt::QueuedConnection);
    connect(this,
        SIGNAL(AddPruned(JobId, std::vector<MolpherMolecule> &, std::string &)),
        &gCommunicator,
        SLOT(AddPruned(JobId, std::vector<MolpherMolecule> &, std::string &)));
    connect(this,
        SIGNAL(SkipNeighborhoodTask(boost::posix_time::ptime)),
        &gCommunicator,
        SLOT(SkipNeighborhoodTask(boost::posix_time::ptime)));
}

ChemicalSpaceView::~ChemicalSpaceView()
{
    FlushTimestamps();
    CleanSreen();
}

void ChemicalSpaceView::wheelEvent(QWheelEvent *event)
{
    int degrees = event->delta() / 8;
    qreal steps = degrees / 60.0;
    ScaleView(pow(2.0, steps));
}

void ChemicalSpaceView::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        setDragMode(QGraphicsView::ScrollHandDrag);
        QMouseEvent modifiedEvent(QEvent::MouseButtonPress, event->pos(),
            Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
        QGraphicsView::mousePressEvent(&modifiedEvent);
    } else if (event->button() == Qt::MiddleButton) {
        bool mouseClickOverMolecule = false;
        QList<QGraphicsItem *> items = this->items(event->x(), event->y());
        foreach (QGraphicsItem *item, items) {
            if (dynamic_cast<VisualizedMolecule *>(item)) {
                mouseClickOverMolecule = true;
                break;
            }
        }
        bool pathHighlightingByMouseHover = mEndMoleculeOfHighlightedPath.empty();
        if (!pathHighlightingByMouseHover && !mouseClickOverMolecule) {
            OnRevertHighlightedPath();
            OnChangeHighlightingRule("");
        }
        QGraphicsView::mousePressEvent(event);
    } else {
        QGraphicsView::mousePressEvent(event);
    }
}

void ChemicalSpaceView::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons().testFlag(Qt::RightButton)) {
        QMouseEvent modifiedEvent(QEvent::MouseMove, event->pos(), Qt::NoButton,
            Qt::LeftButton, Qt::NoModifier);
        QGraphicsView::mouseMoveEvent(&modifiedEvent);
    } else {
        QGraphicsView::mouseMoveEvent(event);
    }
}

void ChemicalSpaceView::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        QMouseEvent modifiedEvent(QEvent::MouseButtonPress, event->pos(),
            Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
        QGraphicsView::mouseReleaseEvent(&modifiedEvent);
        setDragMode(QGraphicsView::RubberBandDrag);
    } else {
        QGraphicsView::mouseReleaseEvent(event);
    }
}

void ChemicalSpaceView::keyPressEvent(QKeyEvent *e)
{
    switch (e->key()) {
    case Qt::Key_Delete:
        DeleteSelectedMolecules();
        break;
    default:
        break;
    }
    QGraphicsView::keyPressEvent(e);
}

void ChemicalSpaceView::ScaleView(qreal factor)
{
    qreal desiredScale =
        matrix().scale(factor, factor).mapRect(QRectF(0, 0, 1, 1)).width();
    if (desiredScale < 0.07 || desiredScale > 5000) {
        return;
    }
    scale(factor, factor);
    mCurrentScale *= factor;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        it->second->scale(1.0 / factor, 1.0 / factor);
        it->second->ScaleEdgeToParent(1.0 / factor);
    }
}

void ChemicalSpaceView::ScaleViewTo(qreal desiredScale)
{
    qreal actualScale = matrix().mapRect(QRectF(0, 0, 1, 1)).width();
    qreal factor = desiredScale / actualScale;
    scale(factor, factor);
    mCurrentScale *= factor;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        it->second->scale(1.0 / factor, 1.0 / factor);
        it->second->ScaleEdgeToParent(1.0 / factor);
    }
}

void ChemicalSpaceView::GoTo(const QPointF &pos)
{
    centerOn(pos);
    ScaleViewTo(CHEMSPACE_GOTO_SCALE);
}

void ChemicalSpaceView::VisualizeSourceAndTarget(const MolpherMolecule &source,
    const MolpherMolecule &target)
{
    if (mMolecules.empty()) {
        CreateMolecule(source, VisualizedMolecule::CD_SOURCE,
            VisualizedMolecule::SD_DEFAULT);
        CreateMolecule(target, VisualizedMolecule::CD_TARGET,
            VisualizedMolecule::SD_DEFAULT);
    } else {
        std::map<std::string, VisualizedMolecule *>::iterator itSource =
            mMolecules.find(source.smile);
        std::map<std::string, VisualizedMolecule *>::iterator itTarget =
            mMolecules.find(target.smile);

        if ((mMolecules.end() != itSource) && (mMolecules.end() != itTarget)) {
            assert(VisualizedMolecule::CD_SOURCE == itSource->second->GetColor());
            assert(VisualizedMolecule::CD_TARGET == itTarget->second->GetColor());
            itSource->second->SetPosition(source.posX, source.posY);
            itTarget->second->SetPosition(target.posX, target.posY);
        }
    }
}

void ChemicalSpaceView::VisualizeDecoys(const std::vector<MolpherMolecule> &decoys)
{
    DeleteOldDecoys(decoys);
    VisualizeActualDecoys(decoys);
}

void ChemicalSpaceView::VisualizeCandidates(const CandidateMap &candidates,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape)
{
    // useful when offline tab is going to open and molecule have to be created,
    //  but its parent is not created too
    std::vector<MolpherMolecule> moleculesToCreateLater;
    std::vector<MolpherMolecule> edgesToParentsToCreateLater;

    CandidateMap::const_iterator it;
    VisualizedMolecule *m;
    for (it = candidates.begin(); it != candidates.end(); ++it) {
        m = GetMolecule(it->first);
        if (NULL != m) {
            switch (m->GetColor()) {
            case (VisualizedMolecule::CD_SOURCE):
                break;
            case (VisualizedMolecule::CD_DECOY):
            case (VisualizedMolecule::CD_TARGET):
                if (!GetMolecule(it->second.parentSmile)) {
                    edgesToParentsToCreateLater.push_back(it->second);
                } else if (m->GetEdgeToParent().isNull()) {
                    CreateEdge(it->second.parentSmile, m,
                        (ChemOperSelector) it->second.parentChemOper);
                }
                break;
            case (VisualizedMolecule::CD_NODE):
            case (VisualizedMolecule::CD_LEAF):
                m->SetPosition(it->second.posX, it->second.posY);
                m->SetColor(GetCandidateColorDefinition(it->second));
                break;
            case (VisualizedMolecule::CD_NEIGHBOR):
            case (VisualizedMolecule::CD_PRUNED):
                CheckPruningCondition(m);
                m->SetPosition(it->second.posX, it->second.posY);
                m->SetColor(GetCandidateColorDefinition(it->second));
                assert(!m->IsConnected());
                if (!GetMolecule(it->second.parentSmile)) {
                    edgesToParentsToCreateLater.push_back(it->second);
                } else {
                    CreateEdge(it->second.parentSmile, m,
                        (ChemOperSelector) it->second.parentChemOper);
                }
                break;
            default:
                assert(false);
                break;
            }

            switch (m->GetShape()) {
            case (VisualizedMolecule::SD_DEFAULT):
            case (VisualizedMolecule::SD_NEIGHBORHOOD_ORIGIN):
                break;
            case (VisualizedMolecule::SD_NEW_CANDIDATE):
                m->SetShape(VisualizedMolecule::SD_DEFAULT);
                break;
            default:
                assert(false);
                break;
            }
        } else {  // add new candidate
            if (!GetMolecule(it->second.parentSmile)) {
                moleculesToCreateLater.push_back(it->second);
                continue;
            }

            VisualizedMolecule::ColorDefinition color =
                GetCandidateColorDefinition(it->second);
            CreateMolecule(it->second, color, newCandidatesShape);
        }
    }

    CreateRestOfCandidates(moleculesToCreateLater, newCandidatesShape);
    CreateRestOfEdges(edgesToParentsToCreateLater);
}

void ChemicalSpaceView::VisualizeOriginAndContext(const MolpherMolecule &origin,
    const std::vector<MolpherMolecule> &context)
{
    // we want origin and context only visualize. Some molecules can be already
    //  pruned, so will not create them

    TryToVisualize(origin);

    std::vector<MolpherMolecule>::const_iterator it;
    for (it = context.begin(); it != context.end(); ++it) {
        TryToVisualize(*it);
    }
}

void ChemicalSpaceView::VisualizeNeighborhood(const std::vector<MolpherMolecule> &neighborhood)
{
    VisualizedMolecule *m;
    std::vector<MolpherMolecule>::const_iterator it;
    for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
        m = GetMolecule(it->smile);
        if (NULL != m) {
            switch (m->GetColor()) {
            case (VisualizedMolecule::CD_SOURCE):
            case (VisualizedMolecule::CD_TARGET):
            case (VisualizedMolecule::CD_DECOY):
            case (VisualizedMolecule::CD_NODE):
            case (VisualizedMolecule::CD_LEAF):
                break;
            case (VisualizedMolecule::CD_NEIGHBOR):
                m->SetPosition(it->posX, it->posY);
                break;
            case (VisualizedMolecule::CD_PRUNED):
                CheckPruningCondition(m);
                m->SetPosition(it->posX, it->posY);
                m->SetColor(VisualizedMolecule::CD_NEIGHBOR);
                break;
            default:
                assert(false);
                break;
            }
        } else {  // add new neighbor
            CreateMolecule(*it, VisualizedMolecule::CD_NEIGHBOR,
                VisualizedMolecule::SD_DEFAULT);
        }
    }
}

void ChemicalSpaceView::PruneMolecules(const PrunedMoleculeVector &prunedMolecules)
{
    PrunedMoleculeVector::const_iterator it;
    std::map<std::string, VisualizedMolecule *>::iterator itToPrune;
    for (it = prunedMolecules.begin(); it != prunedMolecules.end(); ++it) {
        itToPrune = mMolecules.find(*it);
        if (mMolecules.end() != itToPrune) {
            switch (itToPrune->second->GetColor()) {
            case (VisualizedMolecule::CD_NODE):
            case (VisualizedMolecule::CD_LEAF):
                if (mPruneMeansDelete) {
                    Delete(itToPrune->second);
                } else {
                    itToPrune->second->SetColor(VisualizedMolecule::CD_PRUNED);
                    itToPrune->second->DeleteEdges();
                }
                break;
            case (VisualizedMolecule::CD_PRUNED):
                assert(!mPruneMeansDelete);
                assert(!itToPrune->second->IsConnected());
                break;
            case (VisualizedMolecule::CD_SOURCE):
            case (VisualizedMolecule::CD_TARGET):
                break;
            case (VisualizedMolecule::CD_DECOY):
                itToPrune->second->DeleteEdges();
                break;
            default:
                assert(false);
                break;
            }
        }
    }
}

VisualizedMolecule *ChemicalSpaceView::GetMolecule(const std::string &smile)
{
    std::map<std::string, VisualizedMolecule *>::iterator it = mMolecules.find(smile);
    if (mMolecules.end() != it) {
        return (*it).second;
    } else {
        return NULL;
    }
}

void ChemicalSpaceView::SetNewSnapshot(const IterationSnapshot &snp,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape)
{
    VisualizeSourceAndTarget(snp.source, snp.target);
    VisualizeDecoys(snp.decoys);
    VisualizeCandidates(snp.candidates, newCandidatesShape);
    PruneMolecules(snp.prunedDuringThisIter);

    PathHighlightIfHaveToBe();

    // space view rotation, which cause that source and target will be on the
    //  same y-level and target on the right side from the source
    qreal angle = GetRotateAngle(snp.source, snp.target);
    this->rotate(angle);
    QPointer<Edge> edge = NULL;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        it->second->setRotation(-mLastRotateAngleWithoutCorrection);

        edge = it->second->GetEdgeToParent();
        if (!edge.isNull()) {
            edge->RotateText(-mLastRotateAngleWithoutCorrection);
        }
    }

    this->centerOn(snp.source.posX, snp.source.posY);
}

void ChemicalSpaceView::Redraw(const IterationSnapshot &snp)
{
    FlushTimestamps();
    CleanSreen();
    SetNewSnapshot(snp, VisualizedMolecule::SD_DEFAULT);
}

void ChemicalSpaceView::CreateMolecule(const MolpherMolecule &molecule,
    const VisualizedMolecule::ColorDefinition &color,
    const VisualizedMolecule::ShapeDefinition &shape)
{
    VisualizedMolecule *m = new VisualizedMolecule(molecule, color, shape,
        mEndMoleculeOfHighlightedPath.empty());
    connect(m, SIGNAL(ChangeNeighborhoodOrigin(const std::string &)),
        this, SLOT(OnChangeNeighborhoodOrigin(const std::string &)));
    connect(m, SIGNAL(ChangeHighlightingRule(const std::string &)),
        this, SLOT(OnChangeHighlightingRule(const std::string &)));
    connect(m, SIGNAL(RevertHighlightedPath()),
        this, SLOT(OnRevertHighlightedPath()));
    mMolecules.insert(std::make_pair(molecule.smile, m));
    mChemicalSpaceScene->addItem(m);
    m->setScale(1.0 / mCurrentScale);

    if ((VisualizedMolecule::CD_NODE == color) || (VisualizedMolecule::CD_LEAF == color)) {
        CreateEdge(molecule.parentSmile, m, (ChemOperSelector) molecule.parentChemOper);
    }
}

void ChemicalSpaceView::CreateEdge(const std::string &parentSmile,
    VisualizedMolecule *descendant, const ChemOperSelector chemOper)
{
    assert(descendant);
    VisualizedMolecule *parent = GetMolecule(parentSmile);
    assert(parent);
    Edge *edge = new Edge(parent, descendant, chemOper, 1.0 / mCurrentScale);
    mChemicalSpaceScene->addItem(edge);
    mChemicalSpaceScene->addItem(edge->GetText());
}

bool ChemicalSpaceView::ApproximatelyEqual(qreal a, qreal b) const
{
    return (fabs(a - b) <= 32 * DBL_EPSILON * fmax(fabs(a), fabs(b)));
}

void ChemicalSpaceView::CreateRestOfCandidates(std::vector<MolpherMolecule> &molecules,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape)
{
    VisualizedMolecule::ColorDefinition color;
    while (!molecules.empty()) {
        for (size_t i = 0; i < molecules.size(); ++i) {
            if (GetMolecule(molecules[i].parentSmile)) {
                color = GetCandidateColorDefinition(molecules[i]);
                CreateMolecule(molecules[i], color, newCandidatesShape);
                molecules.erase(molecules.begin() + i);
                break;
            }
            assert(molecules.size() - 1 > i);
        }
    }
}

void ChemicalSpaceView::CreateRestOfEdges(std::vector<MolpherMolecule> &molecules)
{
    VisualizedMolecule *descendant;
    for (size_t i = 0; i < molecules.size(); ++i) {
        descendant = GetMolecule(molecules[i].smile);
        assert(descendant);
        CreateEdge(molecules[i].parentSmile, descendant,
            (ChemOperSelector) molecules[i].parentChemOper);
    }
}

void ChemicalSpaceView::GetSelectedMolecules(std::list<std::string> &smiles) const
{
    smiles.clear();

    std::map<std::string, VisualizedMolecule *>::const_iterator it;
    VisualizedMolecule *visualizedMolecule;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        visualizedMolecule = it->second;
        if (visualizedMolecule->isSelected()) {
            smiles.push_back(visualizedMolecule->GetSmiles());
        }
    }
}

void ChemicalSpaceView::PathHighlightIfHaveToBe()
{
    if (mEndMoleculeOfHighlightedPath.empty()) {
        return;
    }

    VisualizedMolecule *pathEnd = GetMolecule(mEndMoleculeOfHighlightedPath);
    if (pathEnd) {
        pathEnd->SetHighlightingOfPathToSource(true);
    }
}

// return positive angle
qreal ChemicalSpaceView::GetRotateAngle(const MolpherMolecule &source,
    const MolpherMolecule &target)
{
    qreal balanceAngle = 0;
    if (ApproximatelyEqual(source.posX, target.posX) ||
           ApproximatelyEqual(source.posY, target.posY)) {
        balanceAngle = GetMultipleOf90DegreeWithoutCorrection(source, target);
    } else {
        balanceAngle = GetNotMultipleOf90DegreeWithoutCorrection(source, target);
    }

    // correction - calculating the final angle to rotate with respect to last iteration
    qreal actualFinalAngleWithoutCorrection = balanceAngle;
    qreal finalAngle = balanceAngle - mLastRotateAngleWithoutCorrection;
    if (finalAngle < 0) {
        finalAngle += 360;
    }
    assert(finalAngle >= 0);
    mLastRotateAngleWithoutCorrection = actualFinalAngleWithoutCorrection;

    return finalAngle;
}

qreal ChemicalSpaceView::GetLength(qreal point1X, qreal point1Y,
    qreal point2X, qreal point2Y) const
{
    qreal length = sqrt(pow(point1X - point2X, 2) + pow(point1Y - point2Y, 2));
    return length;
}

qreal ChemicalSpaceView::GetMultipleOf90DegreeWithoutCorrection(
    const MolpherMolecule& source, const MolpherMolecule& target) const
{
    // x or y of source and target is the same
    if (ApproximatelyEqual(source.posX, target.posX)) {
        if (ApproximatelyEqual(source.posY, target.posY)) {
            // source and target are at the same position
            return 0;
        } else if (source.posY > target.posY) {
            // target is upper on the same x-level as source
            return 90;
        } else {
            // target is lower on the same x-level as source
            return 270;
        }
    }
    if (ApproximatelyEqual(source.posY, target.posY)) {
        if (ApproximatelyEqual(source.posX, target.posX)) {
            // it is not possible - it had to be caught above
            assert(false);
        } else if (source.posX < target.posX) {
            // target is on the right side on the same y-level as source
            return 0;
        } else {
            // target is on the left side on the same y-level as source
            return 180;
        }
    }

    assert(false);
    return -1;
}

qreal ChemicalSpaceView::GetNotMultipleOf90DegreeWithoutCorrection(
    const MolpherMolecule& source, const MolpherMolecule& target) const
{
    // x or y of source and target are not the same

    // check in which quadrant is target if source is center of coordinate system
    int quadrant = 0;
    if (source.posX < target.posX) {
        if (source.posY > target.posY) {
            quadrant = 1;
        } else {
            quadrant = 4;
        }
    } else {
        if (source.posY > target.posY) {
            quadrant = 2;
        } else {
            quadrant = 3;
        }
    }

    // counting angle in source point
    QPointF thirdTrianglePoint(target.posX, source.posY);
    qreal oppositeLegLength = GetLength(target.posX, target.posY,
        thirdTrianglePoint.x(), thirdTrianglePoint.y());
    qreal hypotenuseLength = GetLength(target.posX, target.posY, source.posX,
        source.posY);
    qreal sinSourceAngle = oppositeLegLength / hypotenuseLength;
    qreal sourceAngle = asin(sinSourceAngle) * 180.0 / boost::math::constants::pi<qreal>();

    // counting angle which will be used to rotate the chemical space

    qreal angle = 0;
    switch (quadrant) {
    case 1:
        angle = sourceAngle;
        break;
    case 2:
        angle = 180 - sourceAngle;
        break;
    case 3:
        angle = sourceAngle + 180;
        break;
    case 4:
        angle = 360 - sourceAngle;
        break;
    default:
        assert(false);
        break;
    }
    assert(angle > 0);

    return angle;
}

VisualizedMolecule::ColorDefinition ChemicalSpaceView::GetCandidateColorDefinition(
    const MolpherMolecule &molpherMolecule)
{
    if (molpherMolecule.descendants.empty()) {
        return VisualizedMolecule::CD_LEAF;
    } else {
        return VisualizedMolecule::CD_NODE;
    }
}

void ChemicalSpaceView::OnVisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    std::vector<boost::posix_time::ptime>::iterator itFoundTimestamp;
    itFoundTimestamp = std::find(mNeighborhoodRequestTimestamps.begin(),
        mNeighborhoodRequestTimestamps.end(), res.taskTimestamp);

    // if input timestamp does not belong to this chemical space, exit
    if (mNeighborhoodRequestTimestamps.end() == itFoundTimestamp) {
        return;
    }

    mNeighborhoodRequestTimestamps.erase(itFoundTimestamp);

    VisualizeOriginAndContext(res.origin, res.reducedContext);
    VisualizeNeighborhood(res.reducedNeighborhood);
}

void ChemicalSpaceView::OnPrune(JobId jobId)
{
    std::vector<MolpherMolecule> pruned;
    GetSelectedMolecules(pruned);

    if (pruned.empty()) {
        return;
    }

    std::string password;
    if (PasswordCache::ResolvePassword(jobId, password)) {
        emit AddPruned(jobId, pruned, password);
    } else {
        ShowPermissionError();
    }
}

void ChemicalSpaceView::OnGenerateNeighborhood()
{
    std::vector<MolpherMolecule> context;
    GetSelectedMolecules(context);

    VisualizedMolecule *neighborhoodOrigin = GetMolecule(mNeighborhoodOriginSmiles);
    bool existNeighborhoodOrigin = neighborhoodOrigin &&
        (VisualizedMolecule::SD_NEIGHBORHOOD_ORIGIN == neighborhoodOrigin->GetShape());

    if (!existNeighborhoodOrigin && context.empty()) {
        ShowWarning(tr("Origin of neighborhood or some molecules to recalculate"
            " coordinates for have to be selected."));
        return;
    }

    MolpherMolecule origin;
    if (existNeighborhoodOrigin) {
        std::string smiles = neighborhoodOrigin->GetSmiles();
        std::string formula = neighborhoodOrigin->GetFormula();
        origin = MolpherMolecule(smiles, formula);
    }

    NeighborhoodDialog *neighborhoodDialog = new NeighborhoodDialog(origin, context);
    connect(neighborhoodDialog, SIGNAL(SaveTimestamp(boost::posix_time::ptime)),
        this, SLOT(OnAcceptNeighborhoodDialog(boost::posix_time::ptime)));
    neighborhoodDialog->show();
}

void ChemicalSpaceView::OnRevisualize()
{
    std::vector<MolpherMolecule> context;
    GetSelectedMolecules(context);

    if (context.empty()) {
        ShowWarning(tr("Some molecules to recalculate coordinates for have to be"
            " selected."));
        return;
    }

    RevisualizeDialog *revisualizeDialog = new RevisualizeDialog(context);
    connect(revisualizeDialog, SIGNAL(SaveTimestamp(boost::posix_time::ptime)),
        this, SLOT(OnAcceptNeighborhoodDialog(boost::posix_time::ptime)));
    revisualizeDialog->show();
}

void ChemicalSpaceView::OnBookmark()
{
    VisualizedMolecule *mol;
    std::map<std::string, VisualizedMolecule *>::const_iterator itConst;
    bool warningWasCalled = false;
    for (itConst = mMolecules.begin(); itConst != mMolecules.end(); ++itConst) {
        mol = itConst->second;
        if (mol->isSelected()) {
            std::string smile = mol->GetSmiles();
            std::string formula = mol->GetFormula();
            Bookmarks *b = gGlobalObjectsHolder.GetBookmarks();
            MolpherMolecule self(smile, formula);
            if (!b->Add(self) && !warningWasCalled)  {
                ShowWarning("There is no bookmark group.");
                warningWasCalled = true;
            }
        }
    }
}

void ChemicalSpaceView::OnActionIdentityPubchem()
{
    std::list<std::string> selectedSmiles;
    GetSelectedMolecules(selectedSmiles);

    DbTestingManager::PerformIdentitySearch(selectedSmiles,
        DbTestingManager::DB_PUBCHEM, this);
}

void ChemicalSpaceView::OnActionSimilarityPubchem()
{
    std::list<std::string> selectedSmiles;
    GetSelectedMolecules(selectedSmiles);

    DbTestingManager::PerformSimilaritySearch(selectedSmiles,
        DbTestingManager::DB_PUBCHEM, this);
}

void ChemicalSpaceView::OnChangeNeighborhoodOrigin(const std::string &smiles)
{
    // deselecting of old molecule if it is selected as neighborhood origin
    VisualizedMolecule *neighborhoodOrigin = GetMolecule(mNeighborhoodOriginSmiles);
    if (neighborhoodOrigin) {
        neighborhoodOrigin->DeselectNeighborhoodOrigin();
    }

    mNeighborhoodOriginSmiles = smiles;
}

void ChemicalSpaceView::OnChangeHighlightingRule(const std::string &smiles)
{
    if ((mEndMoleculeOfHighlightedPath.empty() && !smiles.empty()) ||
        (!mEndMoleculeOfHighlightedPath.empty() && smiles.empty())) {
        std::map<std::string, VisualizedMolecule *>::iterator it;
        for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
            it->second->SetPathHighlighting(smiles.empty());
        }
    }

    mEndMoleculeOfHighlightedPath = smiles;
}

void ChemicalSpaceView::OnRevertHighlightedPath()
{
    std::map<std::string, VisualizedMolecule *>::iterator it;
    VisualizedMolecule *visualizedMolecule;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        visualizedMolecule = it->second;
        if (visualizedMolecule->IsEndOfHighlightedPath()) {
            visualizedMolecule->SetHighlightingOfPathToSource(false);
        }
    }
}

void ChemicalSpaceView::OnAcceptNeighborhoodDialog(boost::posix_time::ptime timestamp)
{
    mNeighborhoodRequestTimestamps.push_back(timestamp);
}

void ChemicalSpaceView::FlushTimestamps()
{
    std::vector<boost::posix_time::ptime>::iterator it;
    for (it = mNeighborhoodRequestTimestamps.begin();
            it != mNeighborhoodRequestTimestamps.end(); ++it) {
        emit SkipNeighborhoodTask(*it);
    }
    mNeighborhoodRequestTimestamps.clear();
}

void ChemicalSpaceView::CleanSreen()
{
    // delete of all visualized molecules
    unsigned int moleculesSize = mMolecules.size();
    for (unsigned int i = 0; i < moleculesSize; ++i) {
        Delete(mMolecules.begin()->second);
    }
    assert(mMolecules.empty());
    assert(mChemicalSpaceScene->items().empty());
}

void ChemicalSpaceView::GetSelectedMolecules(
        std::vector<MolpherMolecule> &molecules) const
{
    std::map<std::string, VisualizedMolecule *>::const_iterator it;
    VisualizedMolecule *visualizedMolecule;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        visualizedMolecule = it->second;
        if (visualizedMolecule->isSelected()) {
            std::string smiles = visualizedMolecule->GetSmiles();
            std::string formula = visualizedMolecule->GetFormula();
            molecules.push_back(MolpherMolecule(smiles, formula));
        }
    }
}

void ChemicalSpaceView::DeleteOldDecoys(const std::vector<MolpherMolecule> &actualDecoys)
{
    // old decoys are these decoys which were in last, but are not in actual, snapshot

    // loading old decoys
    std::vector<std::string> smilesListToErase;
    VisualizedMolecule *mol;
    std::map<std::string, VisualizedMolecule *>::const_iterator itConst;
    for (itConst = mMolecules.begin(); itConst != mMolecules.end(); ++itConst) {
        mol = itConst->second;
        if ((mol->GetColor() == VisualizedMolecule::CD_DECOY) &&
                !IsIn(itConst->first, actualDecoys)) {
            if (mol->GetEdgeToParent().isNull()) {
                smilesListToErase.push_back(itConst->first);
            } else {
                // the molecule is both decoy and candidate
                QList<QPointer<Edge> > edgesToDescendants;
                mol->GetEdgesToDescendants(edgesToDescendants);
                if (edgesToDescendants.empty()) {
                    mol->SetColor(VisualizedMolecule::CD_LEAF);
                } else {
                    mol->SetColor(VisualizedMolecule::CD_NODE);
                }
            }
        }
    }

    // erasing old decoys
    Delete(smilesListToErase);
}

void ChemicalSpaceView::DeleteSelectedMolecules()
{
    // manual delete is possible only for neighborhood nodes and pruned molecules

    // loading molecules to delete
    std::vector<std::string> smilesListToErase;
    VisualizedMolecule::ColorDefinition color;
    bool isRightColor = false;
    std::map<std::string, VisualizedMolecule *>::const_iterator itConst;
    for (itConst = mMolecules.begin(); itConst != mMolecules.end(); ++itConst) {
        color = itConst->second->GetColor();
        isRightColor = (VisualizedMolecule::CD_NEIGHBOR == color) ||
            (VisualizedMolecule::CD_PRUNED == color);
        if (isRightColor && itConst->second->isSelected()) {
            smilesListToErase.push_back(itConst->first);
        }
    }

    // erasing molecules to delete
    Delete(smilesListToErase);
}

void ChemicalSpaceView::Delete(const std::vector<std::string> &smilesList)
{
    std::string smiles;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (size_t i = 0; i < smilesList.size(); ++i) {
        smiles = smilesList[i];
        it = mMolecules.find(smiles);
        if (mMolecules.end() != it) {
            Delete(it->second);
        }
    }
}

void ChemicalSpaceView::Delete(VisualizedMolecule *molecule)
{
    assert(molecule);

    QPointer<Edge> edgeToParent = molecule->GetEdgeToParent();
    if (!edgeToParent.isNull()) {
        mChemicalSpaceScene->removeItem(edgeToParent.data());
        mChemicalSpaceScene->removeItem(edgeToParent->GetText());
    }
    QList<QPointer<Edge> > edgesToDescendants;
    molecule->GetEdgesToDescendants(edgesToDescendants);
    foreach (QPointer<Edge> edgeToDescendant, edgesToDescendants) {
        mChemicalSpaceScene->removeItem(edgeToDescendant.data());
        mChemicalSpaceScene->removeItem(edgeToDescendant->GetText());
    }

    mChemicalSpaceScene->removeItem(molecule);
    mMolecules.erase(molecule->GetSmiles());
    delete molecule;
}

void ChemicalSpaceView::VisualizeActualDecoys(
    const std::vector<MolpherMolecule> &actualDecoys)
{
    std::vector<MolpherMolecule>::const_iterator it;
    VisualizedMolecule *m;
    for (it = actualDecoys.begin(); it != actualDecoys.end(); ++it) {
        m = GetMolecule(it->smile);
        if (NULL != m) {
            switch (m->GetColor()) {
            case (VisualizedMolecule::CD_SOURCE):
            case (VisualizedMolecule::CD_TARGET):
                break;
            case (VisualizedMolecule::CD_NODE):
            case (VisualizedMolecule::CD_LEAF):
            case (VisualizedMolecule::CD_NEIGHBOR):
            case (VisualizedMolecule::CD_PRUNED):
                CheckPruningCondition(m);
                m->SetPosition(it->posX, it->posY);
                m->SetColor(VisualizedMolecule::CD_DECOY);
                break;
            case (VisualizedMolecule::CD_DECOY):
                m->SetPosition(it->posX, it->posY);
                break;
            default:
                assert(false);
                break;
            }
        } else {  // add new decoy
            CreateMolecule(*it, VisualizedMolecule::CD_DECOY,
                VisualizedMolecule::SD_DEFAULT);
        }
    }
}

bool ChemicalSpaceView::IsIn(const std::string &smiles,
    const std::vector<MolpherMolecule> &molecules) const
{
    std::vector<MolpherMolecule>::const_iterator it;
    for (it = molecules.begin(); it != molecules.end(); ++it) {
        if (it->smile == smiles) {
            return true;
        }
    }

    return false;
}

void ChemicalSpaceView::TryToVisualize(const MolpherMolecule &molecule)
{
    VisualizedMolecule *visMol = GetMolecule(molecule.smile);
    if (visMol) {
        if (visMol->GetColor() != VisualizedMolecule::CD_PRUNED) {
            visMol->SetPosition(molecule.posX, molecule.posY);
        }

        // row below is useful when similar molecule's coordinates are
        //  calculated for the first time
        visMol->setVisible(true);
    }
}

void ChemicalSpaceView::CheckPruningCondition(const VisualizedMolecule *molecule) const
{
    bool isAlwaysFalse = molecule &&
        (molecule->GetColor() == VisualizedMolecule::CD_PRUNED) &&
        mPruneMeansDelete;
    assert(!isAlwaysFalse);
}

void ChemicalSpaceView::OnAcceptIdentitySearchResult(std::string sourceMol,
    std::list<std::string> &result)
{
    if (result.empty() ||
            (std::find(result.begin(), result.end(), sourceMol) == result.end())) {
        return;
    }

    VisualizedMolecule *mol = GetMolecule(sourceMol);
    if (mol) {
        mol->SetShape(VisualizedMolecule::SD_IDENTITY_SEARCH);
    }
}

void ChemicalSpaceView::OnAcceptSimilaritySearchResult(std::string sourceMol,
    std::list<std::string> &result)
{
    if (result.empty()) {
        return;
    }

    std::vector<MolpherMolecule> simMolpherMolecules;
    MolpherMolecule simMolpherMolecule;
    VisualizedMolecule *simVisualizedMolecule;
    std::list<std::string>::const_iterator it;
    for (it = result.begin(); it != result.end(); ++it) {
        RDKit::RWMol *mol = NULL;
        try {
            mol = RDKit::SmilesToMol(*it);
            if (mol) {
                RDKit::MolOps::Kekulize(*mol);
                simMolpherMolecule.formula = RDKit::Descriptors::calcMolFormula(*mol);
                simMolpherMolecule.smile = RDKit::MolToSmiles(*mol);
            } else {
                throw ValueErrorException("");
            }
        } catch (const ValueErrorException &exc) {
            // ignore molecule from result
            continue;
        }
        delete mol;

        if (!GetMolecule(simMolpherMolecule.smile)) {
            simMolpherMolecules.push_back(simMolpherMolecule);
            CreateMolecule(simMolpherMolecule, VisualizedMolecule::CD_NEIGHBOR,
                VisualizedMolecule::SD_DEFAULT);

            simVisualizedMolecule = GetMolecule(simMolpherMolecule.smile);
            assert(simVisualizedMolecule);
            simVisualizedMolecule->setVisible(false);
        }
    }

    if (!simMolpherMolecules.empty()) {
        NeighborhoodTask task;
        task.taskTimestamp = boost::posix_time::microsec_clock::universal_time();
        task.context = simMolpherMolecules;

        mNeighborhoodRequestTimestamps.push_back(task.taskTimestamp);
        emit RevisualizeSimilarMolecules(task);
    }
}

void ChemicalSpaceView::OnSelectAll()
{
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        it->second->setSelected(true);
    }
    setFocus();
}

void ChemicalSpaceView::OnSelect(VisualizedMolecule::ColorDefinition color)
{
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        if (it->second->GetColor() == color) {
            it->second->setSelected(true);
        }
    }
    setFocus();
}

void ChemicalSpaceView::OnSelect(VisualizedMolecule::ShapeDefinition shape)
{
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        if (it->second->GetShape() == shape) {
            it->second->setSelected(true);
        }
    }
    setFocus();
}

void ChemicalSpaceView::OnGoToSource()
{
    VisualizedMolecule *visualizedMolecule = NULL;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        if (it->second->GetColor() == VisualizedMolecule::CD_SOURCE) {
            visualizedMolecule = it->second;
            break;
        }
    }

    assert(visualizedMolecule);

    GoTo(visualizedMolecule->pos());
}

void ChemicalSpaceView::OnGoToTarget()
{
    VisualizedMolecule *visualizedMolecule = NULL;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        if (it->second->GetColor() == VisualizedMolecule::CD_TARGET) {
            visualizedMolecule = it->second;
            break;
        }
    }

    assert(visualizedMolecule);

    GoTo(visualizedMolecule->pos());
}

void ChemicalSpaceView::OnGoToNeighborhoodOrigin()
{
    VisualizedMolecule *visualizedMolecule = NULL;
    std::map<std::string, VisualizedMolecule *>::iterator it;
    for (it = mMolecules.begin(); it != mMolecules.end(); ++it) {
        if (it->second->GetShape() == VisualizedMolecule::SD_NEIGHBORHOOD_ORIGIN) {
            visualizedMolecule = it->second;
            break;
        }
    }

    if (visualizedMolecule) {
        GoTo(visualizedMolecule->pos());
    }
}
