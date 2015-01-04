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

#include "PubchemSimilarityTester.h"

#define ERROR_THRESHOLD 3

PubchemSimilarityTester::PubchemSimilarityTester()
{
    mTimer = new QTimer(this);
    connect(mTimer, SIGNAL(timeout()), this, SLOT(WaitingTimeout()));
    mTimer->setInterval(2000);
    mTimer->start();
}

PubchemSimilarityTester::~PubchemSimilarityTester()
{

    mTimer->stop();
    std::map<std::string, CidTranslationResult *>::iterator it;
    while (!mCidTranslationMap.empty()) {
        it = mCidTranslationMap.begin();
        delete it->second;
        mCidTranslationMap.erase(it);
    }
}

void PubchemSimilarityTester::SimilaritySearch(const std::list<std::string> &smileList)
{
    mQueryList = smileList;
    mIdentitySearch = false;
    SendNextQuery();
}

void PubchemSimilarityTester::IdentitySearch(const std::list<std::string> &smileList)
{
    mQueryList = smileList;
    mIdentitySearch = true;
    SendNextQuery();
}

void PubchemSimilarityTester::SendSimilarityQuery()
{
    SynchCout("Pubchem similarity tester: Send similarity query");
    assert(!mIdentitySearch);

    mCurrentSmile = mQueryList.front();
    mQueryList.pop_front();
    std::stringstream url;
    url << "http://pubchem.ncbi.nlm.nih.gov/rest/pug";
    url << "/compound";
    url << "/similarity";
    url << "/smiles";
    url << "/" << QByteArray(mCurrentSmile.c_str()).toPercentEncoding().constData();
    url << "/XML";
    url << "?Treshold=95&MaxRecords=10";

    SynchCout(url.str());
    QNetworkRequest request(QUrl::fromEncoded(QByteArray(url.str().c_str())));
    QNetworkReply *reply = mNetworkManager.get(request);

    connect(reply, SIGNAL(finished()), this, SLOT(Reply()));
}

void PubchemSimilarityTester::SendIdentityQuery()
{
    SynchCout("Pubchem similarity tester: Send identity query");
    assert(mIdentitySearch);

    std::stringstream url;
    mCurrentSmile = mQueryList.front();
    mQueryList.pop_front();
    url << "http://pubchem.ncbi.nlm.nih.gov/rest/pug";
    url << "/compound";
    url << "/identity";
    url << "/smiles";
    url << "/" << QByteArray(mCurrentSmile.c_str()).toPercentEncoding().constData();
    url << "/XML";

    QNetworkRequest request(QUrl::fromEncoded(QByteArray(url.str().c_str())));
    QNetworkReply *reply = mNetworkManager.get(request);

    connect(reply, SIGNAL(finished()), this, SLOT(Reply()));
}

void PubchemSimilarityTester::CidReply()
{
    SynchCout("Pubchem similarity tester: Cid reply received");
    QNetworkReply *reply = qobject_cast<QNetworkReply *>(sender());

    if (reply) {
        if (reply->error() == QNetworkReply::NoError) {
            QString content(reply->readAll());
            QDomDocument doc;
            doc.setContent(content);
            std::string smile;
            std::string cid;
            ParseSmilesProperty(doc, cid, smile);
            StoreTranslatedCid(cid, smile);
        } else {
            //get http status code
            int httpStatus = reply->attribute(
                QNetworkRequest::HttpStatusCodeAttribute).toInt();

            std::stringstream s;
            s << "HTTP error: " << httpStatus;
            SynchCout(s.str());
            SynchCout(reply->errorString().toStdString());
            //do some error management
        }

        reply->deleteLater();
    }
}

void PubchemSimilarityTester::Reply()
{
    SynchCout("Pubchem similarity tester: Reply received");
    QNetworkReply *reply = qobject_cast<QNetworkReply *>(sender());

    if (reply) {
        if (reply->error() == QNetworkReply::NoError) {

            QString content = (QString)reply->readAll();
            QDomDocument doc;
            doc.setContent(content);

            std::string rootElement;
            QDomElement root = doc.documentElement();
            rootElement = root.tagName().toStdString();

            std::string listKey;
            if (rootElement == "Waiting") {
                ParseListKey(doc, listKey);

                // If the received listkey is not registered yet, register it
                // along with the current smile
                if (mListKeys.find(listKey) == mListKeys.end()) {
                    std::stringstream s;
                    s << "Pubchem similarity tester: Registering list key " << listKey;
                    SynchCout(s.str());
                    mListKeys.insert(std::make_pair(listKey, mCurrentSmile));
                    mErrorCounter.insert(std::make_pair(listKey, 0));
                    SendNextQuery();
                }
            } else {
                std::list<std::string> cidList;
                ParseCids(doc, cidList);

                std::string reqListKey;
                ParseRequestListKey(reply->request(), reqListKey);
                if (reqListKey != "") {

                    std::map<std::string, std::string>::iterator it =
                        mListKeys.find(reqListKey);
                    if (it != mListKeys.end()) {
                        std::string smile = it->second;
                        mListKeys.erase(it);
                        std::stringstream s;
                        s << "Pubchem similarity tester: Reply received for " << "listkey: " << reqListKey << ", molecule: " << smile;
                        SynchCout(s.str());
                        TranslateCids(smile, cidList);
                    } else {
                        SynchCout("Pubchem similarity tester: Reply received for non-registered list key");
                    }

                } else {
                    SynchCout("Pubchem similarity tester: Request list key could not be parsed");
                }

                std::stringstream s;
                s << "Pubchem similarity tester: Cids parsed, count: " << cidList.size();
                SynchCout(s.str());
            }

        } else {
            //get http status code
            int httpStatus = reply->attribute(
                QNetworkRequest::HttpStatusCodeAttribute).toInt();

            std::stringstream s;
            s << "HTTP error: " << httpStatus;
            SynchCout(s.str());
            SynchCout(reply->errorString().toStdString());
            //do some error management

            std::string reqListKey;
            ParseRequestListKey(reply->request(), reqListKey);
            if (reqListKey != "") {
                std::map<std::string, int>::iterator it = mErrorCounter.find(reqListKey);
                if (it != mErrorCounter.end()) {
                    it->second++;
                    SynchCout("Pubchem similarity tester: Reply error, increasing errcount");
                    if (it->second > ERROR_THRESHOLD) {
                        std::stringstream s;
                        s << "Pubchem similarity tester: Reply threshold reached, removing list key " << reqListKey;
                        SynchCout(s.str());
                        RemoveListKey(reqListKey);
                    }
                }
            }
        }
        reply->deleteLater();
    }
}

void PubchemSimilarityTester::ParseCids(const QDomDocument &doc,
    std::list<std::string> &cidList)
{
    QDomElement root = doc.documentElement();
    if (root.tagName() != "IdentifierList") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid root element. " <<
            "Expected: IdentifierList, Found: " << root.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    QDomElement cidElement = root.firstChild().toElement();
    while (!cidElement.isNull()) {
        if (cidElement.tagName() != "CID") {
            std::stringstream s;
            s << "Pubchem similarity tester: Returned document has invalid element. " <<
                "Expected: CID Found: " << cidElement.tagName().toStdString();
            SynchCout(s.str());

        }
        cidList.push_back(cidElement.text().toStdString());
        cidElement = cidElement.nextSibling().toElement();
    }
}

void PubchemSimilarityTester::ParseListKey(const QDomDocument &doc, std::string &listKey)
{
    QDomElement root = doc.documentElement();
    if (root.tagName() != "Waiting") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid root element. " <<
            "Expected: Waiting, Found: " << root.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    QDomElement listKeyElem = root.firstChild().toElement();
    if (listKeyElem.tagName() != "ListKey") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid element. " <<
            "Expected: ListKey Found: " << listKeyElem.tagName().toStdString();
        SynchCout(s.str());
    }

    listKey = listKeyElem.text().toStdString();
}

void PubchemSimilarityTester::AskForResult(std::string listKey)
{
    std::stringstream url;
    url << "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/";
    url << listKey << "/cids";
    url << "/XML";

    std::stringstream s;
    s << "Pubchem similarity tester: Asking for result of listkey " << listKey;
    SynchCout(s.str());

    QUrl qUrl(url.str().c_str());
    QNetworkRequest request(qUrl);
    QNetworkReply *reply = mNetworkManager.get(request);

    connect(reply, SIGNAL(finished()), this, SLOT(Reply()));
}

void PubchemSimilarityTester::WaitingTimeout()
{
    std::map<std::string, std::string>::iterator it = mListKeys.begin();
    while (it != mListKeys.end()) {
        AskForResult(it->first);
        it++;
    }
}

void PubchemSimilarityTester::ParseRequestListKey(
    const QNetworkRequest &request, std::string &listKey)
{
    QUrl qUrl = request.url();
    std::string url = qUrl.path().toStdString();

    std::string listKeyPattern("listkey/");
    unsigned int start = url.find(listKeyPattern);
    if (start == url.npos) {
        listKey = "";
        return;
    }
    start = start + listKeyPattern.size();

    std::string tmp = url.substr(start);

    unsigned int end = tmp.find("/");
    if (end == tmp.npos) {
        listKey = "";
        return;
    }

    listKey = tmp.substr(0, end);
}

void PubchemSimilarityTester::SendNextQuery()
{
    if (mQueryList.empty()) {
        return;
    }

    if (mIdentitySearch) {
        SendIdentityQuery();
    } else {
        SendSimilarityQuery();
    }
}

void PubchemSimilarityTester::TranslateCids(
    std::string &smile, std::list<std::string> &cidList)
{
    std::map<std::string, std::string> cidMap;

    std::list<std::string>::iterator it = cidList.begin();
    while (it != cidList.end()) {
        std::stringstream url;
        url << "http://pubchem.ncbi.nlm.nih.gov/rest/pug/";
        url << "compound/cid/";
        url << *it << "/";
        url << "property/CanonicalSMILES/XML";

        cidMap.insert(std::make_pair(*it, ""));

        QUrl qUrl(url.str().c_str());
        QNetworkRequest request(qUrl);
        QNetworkReply *reply = mNetworkManager.get(request);

        connect(reply, SIGNAL(finished()), this, SLOT(CidReply()));

        it++;
    }

    CidTranslationResult *res = new CidTranslationResult(smile, cidMap);
    connect(res, SIGNAL(ResultCompleted(std::string)),
        this, SLOT(ResultCompleted(std::string)));
    mCidTranslationMap.insert(std::make_pair(smile, res));
}

void PubchemSimilarityTester::ParseSmilesProperty(
    const QDomDocument &doc, std::string &cid, std::string &smile)
{
    QDomElement root = doc.documentElement();
    if (root.tagName() != "PropertyTable") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid element. " <<
            "Expected: PropertyTable Found: " << root.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    QDomElement propertiesElem = root.firstChildElement("Properties");
    if (propertiesElem.isNull() || propertiesElem.tagName() != "Properties") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid element. " <<
            "Expected: Properties Found: " << propertiesElem.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    QDomElement cidElem = propertiesElem.firstChildElement("CID");
    if (cidElem.isNull() || cidElem.tagName() != "CID") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid element. " <<
            "Expected: CID Found: " << cidElem.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    cid = cidElem.text().toStdString();

    QDomElement smilesElem = propertiesElem.firstChildElement("CanonicalSMILES");
    if (smilesElem.isNull() || smilesElem.tagName() != "CanonicalSMILES") {
        std::stringstream s;
        s << "Pubchem similarity tester: Returned document has invalid element. " <<
            "Expected: CanonicalSMILES Found: " << smilesElem.tagName().toStdString();
        SynchCout(s.str());
        return;
    }

    smile = smilesElem.text().toStdString();
}

void PubchemSimilarityTester::StoreTranslatedCid(
    const std::string &cid, const std::string &smile)
{
    std::map<std::string, CidTranslationResult *>::iterator it =
        mCidTranslationMap.begin();
    CidTranslationResult *res;
    while (it != mCidTranslationMap.end()) {
        res = it->second;
        it++;

        if (res->ContainsCid(cid)) {
            res->AddSmile(cid, smile);
        }
    }
}

void PubchemSimilarityTester::ResultCompleted(std::string smile)
{
    SynchCout("Pubchem similarity tester: Result completed for molecule " + smile );
    std::map<std::string, CidTranslationResult *>::iterator it =
        mCidTranslationMap.find(smile);
    if (it == mCidTranslationMap.end()) {
        SynchCout("Pubchem similarity tester: Translation result not found in the translation map");
        return;
    }

    std::list<std::string> smileList;
    CidTranslationResult *result = it->second;
    result->GetResultSmiles(smileList);

    if (mIdentitySearch) {
        emit PublishIdentityResult(smile, smileList);
    } else {
        emit PublishSimilarityResult(smile, smileList);
    }

    mCidTranslationMap.erase(it);
    result->deleteLater();
}

void PubchemSimilarityTester::RemoveListKey(const std::string& listKey)
{
    std::map<std::string, std::string>::iterator it = mListKeys.find(listKey);
    if (it != mListKeys.end()) {
        mListKeys.erase(it);
    }

    std::map<std::string, int>::iterator errIt = mErrorCounter.find(listKey);
    if (errIt != mErrorCounter.end()) {
        mErrorCounter.erase(errIt);
    }
}

CidTranslationResult::CidTranslationResult(const std::string &sourceSmile,
    const std::map<std::string, std::string> &cidSmileMap
    ) :
    mResolvedCount(0),
    mSourceSmile(sourceSmile)
{
    mCidSmileMap = cidSmileMap;
    mTimer = new QTimer(this);
    unsigned int timeout = (floor(mCidSmileMap.size() / 10) + 1) * 1000;
    mTimer->setInterval(timeout);
    connect(mTimer, SIGNAL(timeout()), this, SLOT(ResultTimeout()));
    mTimer->start();
}

CidTranslationResult::~CidTranslationResult()
{
}

bool CidTranslationResult::ContainsCid(const std::string &cid)
{
    std::map<std::string, std::string>::iterator it = mCidSmileMap.find(cid);
    return it != mCidSmileMap.end();
}

void CidTranslationResult::AddSmile(const std::string &cid, const std::string &smile)
{
    std::map<std::string, std::string>::iterator it = mCidSmileMap.find(cid);

    if (it == mCidSmileMap.end() || it->second != "") {
        return;
    }

    mCidSmileMap.erase(it);
    mCidSmileMap.insert(std::make_pair(cid, smile));
    mResolvedCount++;

    if (IsAllFinished()) {
        SynchCout("Pubchem similarity tester: Testing finished, emiting result");
        PublishResult();
    }
}

bool CidTranslationResult::IsAllFinished()
{
    return mResolvedCount == mCidSmileMap.size();
}

void CidTranslationResult::ResultTimeout()
{
    std::stringstream s;
    s << "Pubchem similarity tester: Translation timeout - result ";
    if (IsAllFinished()) {
        s << " IS FINISHED";
    } else {
        s << " NOT FINISHED";
    }
    SynchCout(s.str());

    PublishResult();
}

void CidTranslationResult::GetResultSmiles(std::list<std::string> &resultList)
{
    std::map<std::string, std::string>::iterator it;
    for (it = mCidSmileMap.begin(); it != mCidSmileMap.end(); ++it) {
        resultList.push_back(it->second);
    }
}

void CidTranslationResult::PublishResult()
{
    mTimer->stop();
    emit ResultCompleted(mSourceSmile);
}
