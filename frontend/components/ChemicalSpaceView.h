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

#include <map>

#include <QtCore/QObject>
#include <QtGui/QGraphicsView>
#include <QtGui/QWheelEvent>
#include <QtGui/QMouseEvent>

#include "global_types.h"
#include "IterationSnapshot.h"
#include "NeighborhoodTask.h"
#include "VisualizedMolecule.h"

#define CHEMSPACE_GOTO_SCALE 10.0

class Tab;

class ChemicalSpaceView : public QGraphicsView
{
    Q_OBJECT

public:
    ChemicalSpaceView(IterationSnapshot &snapshot, bool pruneMeansDelete);
    ~ChemicalSpaceView();

    void SetNewSnapshot(const IterationSnapshot &snp,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape);
    void Redraw(const IterationSnapshot &snp);

public slots:
    // slots used by communicator
    void OnVisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);

    void OnAcceptNeighborhoodDialog(boost::posix_time::ptime timestamp);
    void OnAcceptIdentitySearchResult(std::string sourceMol, std::list<std::string> &result);
    void OnAcceptSimilaritySearchResult(std::string sourceMol, std::list<std::string> &result);
    void GetSelectedMolecules(std::vector<MolpherMolecule> &molecules) const;

    // slots used by visualized molecule
    void OnChangeNeighborhoodOrigin(const std::string &smiles);
    void OnChangeHighlightingRule(const std::string &smiles);
    void OnRevertHighlightedPath();

    // slots used by tab
    void OnSelectAll();
    void OnSelect(VisualizedMolecule::ColorDefinition color);
    void OnSelect(VisualizedMolecule::ShapeDefinition shape);
    void OnGoToSource();
    void OnGoToTarget();
    void OnGoToNeighborhoodOrigin();
    void OnPrune(JobId jobId);
    void OnGenerateNeighborhood();
    void OnRevisualize();
    void OnBookmark();
    void OnActionIdentityPubchem();
    void OnActionSimilarityPubchem();

protected:
    typedef IterationSnapshot::CandidateMap CandidateMap;
    typedef IterationSnapshot::PrunedMoleculeVector PrunedMoleculeVector;

    // Functions called only internally.
    void wheelEvent(QWheelEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *e);
    void ScaleView(qreal factor);
    void ScaleViewTo(qreal desiredScale);
    void GoTo(const QPointF &pos);
    void VisualizeSourceAndTarget(const MolpherMolecule &source, const MolpherMolecule &target);
    void VisualizeDecoys(const std::vector<MolpherMolecule> &decoys);
    void VisualizeCandidates(const CandidateMap &candidates,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape);
    void VisualizeOriginAndContext(const MolpherMolecule &origin,
        const std::vector<MolpherMolecule> &context);
    void VisualizeNeighborhood(const std::vector<MolpherMolecule> &neighborhood);
    void PruneMolecules(const PrunedMoleculeVector &prunedMolecules);
    VisualizedMolecule *GetMolecule(const std::string &smile);
    void FlushTimestamps();
    void CleanSreen(); // if redraw another snapshot, we want delete all molecules
    void DeleteOldDecoys(const std::vector<MolpherMolecule> &actualDecoys);
    void DeleteSelectedMolecules();
    void Delete(const std::vector<std::string> &smilesList);
    void Delete(VisualizedMolecule *molecule);
    void VisualizeActualDecoys(const std::vector<MolpherMolecule> &actualDecoys);
    bool IsIn(const std::string &smiles, const std::vector<MolpherMolecule> &molecules) const;
    void TryToVisualize(const MolpherMolecule &molecule);
    void CheckPruningCondition(const VisualizedMolecule *molecule) const;
    void CreateMolecule(const MolpherMolecule &molecule,
        const VisualizedMolecule::ColorDefinition &color,
        const VisualizedMolecule::ShapeDefinition &shape);
    void CreateEdge(const std::string &parentSmile, VisualizedMolecule *descendant,
        const ChemOperSelector chemOper);
    qreal GetRotateAngle(const MolpherMolecule &source, const MolpherMolecule &target);
    qreal GetLength(qreal point1X, qreal point1Y, qreal point2X, qreal point2Y) const;
    qreal GetMultipleOf90DegreeWithoutCorrection(const MolpherMolecule &source,
        const MolpherMolecule &target) const;
    qreal GetNotMultipleOf90DegreeWithoutCorrection(const MolpherMolecule &source,
        const MolpherMolecule &target) const;
    VisualizedMolecule::ColorDefinition GetCandidateColorDefinition(
        const MolpherMolecule &molpherMolecule);
    bool ApproximatelyEqual(qreal a, qreal b) const;
    void CreateRestOfCandidates(std::vector<MolpherMolecule> &molecules,
        const VisualizedMolecule::ShapeDefinition newCandidatesShape);
    void CreateRestOfEdges(std::vector<MolpherMolecule> &molecules);
    void GetSelectedMolecules(std::list<std::string> &smiles) const;
    void PathHighlightIfHaveToBe();

signals:
    void AddPruned(JobId jobId, std::vector<MolpherMolecule> &pruned, std::string &password);
    void SkipNeighborhoodTask(boost::posix_time::ptime timestamp);
    void RevisualizeSimilarMolecules(NeighborhoodTask &task);

private:
    QGraphicsScene *mChemicalSpaceScene;
    // find method is logarithmic, i can not copy QObject
    std::map<std::string, VisualizedMolecule *> mMolecules;
    std::string mNeighborhoodOriginSmiles;
    std::vector<boost::posix_time::ptime> mNeighborhoodRequestTimestamps;
    bool mPruneMeansDelete;
    // it is an angle from visualization of last snapshot, which should be used
    //  when transformation matrix is not modified (space is not rotated)
    qreal mLastRotateAngleWithoutCorrection;
    qreal mCurrentScale;
    // empty string means path highlighting by mouse hover
    std::string mEndMoleculeOfHighlightedPath;
};
