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

#include <QtCore/QObject>
#include <QtNetwork/QNetworkAccessManager>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QNetworkRequest>
#include <QtCore/QUrl>

#include <DataStructs/BitOps.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "MolpherMolecule.h"
#include "inout.h"

class SimilarityTester : public QObject
{
    Q_OBJECT

public:
    SimilarityTester() {}
    virtual ~SimilarityTester() {}

public slots:
    virtual void SimilaritySearch(const std::list<std::string> &smileList) = 0;
    virtual void IdentitySearch(const std::list<std::string> &smileList) = 0;

signals:
    void PublishSimilarityResult(std::string testedMol, std::list<std::string> &result);
    void PublishIdentityResult(std::string testedMol, std::list<std::string> &result);

protected:
    QNetworkAccessManager mNetworkManager;
};
