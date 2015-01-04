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

#include <QtXml/QDomDocument>
#include <QtCore/QTimer>

#include "SimilarityTester.h"

class CidTranslationResult;

class PubchemSimilarityTester : public SimilarityTester
{
    Q_OBJECT

public:
    PubchemSimilarityTester();
    ~PubchemSimilarityTester();

public slots:
    virtual void SimilaritySearch(const std::list<std::string> &smileList);
    virtual void IdentitySearch(const std::list<std::string> &smileList);
    void Reply();
    void CidReply();
    void WaitingTimeout();
    void ResultCompleted(std::string smile);

protected:
    void ParseListKey(const QDomDocument &doc, std::string &listKey);
    void AskForResult(std::string listKey);
    void SendNextQuery();
    void SendSimilarityQuery();
    void SendIdentityQuery();
    void ParseCids(const QDomDocument &doc, std::list<std::string> &cidList);
    void ParseRequestListKey(const QNetworkRequest &request, std::string &listKey);
    void ParseSmilesProperty(const QDomDocument &doc, std::string &cid, std::string &smile);
    void TranslateCids(std::string &smile, std::list<std::string> &cidList);
    void StoreTranslatedCid(const std::string &cid, const std::string &smile);
    void RemoveListKey(const std::string &listKey);

private:
    QTimer *mTimer;
    // Map <listKey, smiles>
    std::map<std::string, std::string> mListKeys;
    // Map <smiles, result>
    std::map<std::string, CidTranslationResult *> mCidTranslationMap;
    // Map <listKey, errorCount>
    std::map<std::string, int> mErrorCounter;
    std::list<std::string> mQueryList;
    std::string mCurrentSmile;
    bool mIdentitySearch;
};

class CidTranslationResult : public QObject
{
    Q_OBJECT

public:
    CidTranslationResult(const std::string &sourceSmile,
        const std::map<std::string, std::string> &cidSmileMap);
    ~CidTranslationResult();
    bool ContainsCid(const std::string &cid);
    void AddSmile(const std::string &cid, const std::string &smile);
    bool IsAllFinished();
    void GetResultSmiles(std::list<std::string> &resultList);
    void PublishResult();

public slots:
    void ResultTimeout();

signals:
    void ResultCompleted(std::string smile);

protected:
    unsigned int mResolvedCount;
    std::string mSourceSmile;
    // Map <cid, smile>
    std::map<std::string, std::string> mCidSmileMap;
    QTimer *mTimer;
};
