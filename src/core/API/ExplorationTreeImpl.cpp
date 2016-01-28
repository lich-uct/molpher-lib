
//void ExplorationTree::ExplorationTreeImpl::updateFromData(ExplorationData &data)
//{
//    if (!data.isValid()) {
//        std::runtime_error("Supplied exploration data is invalid.");
//    }
//    
//    candidates.clear();
//    for (auto& mol_data : data.candidates) {
//        candidates.push_back(std::make_shared<MolpherMol::MolpherMolImpl>(mol_data));
//    }
//    
//    candidatesMask.swap(data.candidatesMask)
//    chemOpers.swap(data.chemOpers)
//    fingerprint = data.fingerprint;
//    generationCnt = data.generationCnt;
//    
//    if (morphDerivations.empty()) {
//        for (auto& mol_data : data.morphDerivations) {
//            morphDerivations.insert(mol_data);
//        }
//    }
//    
//    params = data.params
//    simCoeff = data.simCoeff;
//    source = std::make_shared<MolpherMol::MolpherMolImpl>(data.source);
//    target = std::make_shared<MolpherMol::MolpherMolImpl>(data.target);
//    
//    if (treeMap.empty()) {
//        for (auto& mol_data : data.treeMap) {
//            treeMap.insert(mol_data);
//        }
//    }
//    
//}
//
//std::shared_ptr<ExplorationData> ExplorationTree::ExplorationTreeImpl::asData() {
//    auto data = std::make_shared<ExplorationData>();
//    
//    for (auto& mol : candidates) {
//        data->candidates.push_back(*(mol->asData()));
//    }
//    
//    data->candidatesMask.swap(candidatesMask)
//    data->chemOpers.swap(chemOpers)
//    data->fingerprint = fingerprint;
//    data->generationCnt = generationCnt;
//    
//    for (auto& item : morphDerivations) {
//        data->morphDerivations.insert(item);
//    }
//    
//    data->params = params;
//    data->simCoeff = simCoeff;
//    data->source = *(source->asData());
//    data->target = *(target->asData());
//    
//    for (auto& item : treeMap) {
//        data->treeMap.insert(std::make_pair<std::string, MolpherMolData>(item.first, *(item.second->asData())));
//    }
//}



//void PathFinderContext::clear()
//{
//    chemOperSelectors.clear();
//    decoys.clear();
//    candidates.clear();
//    morphDerivations.clear();
//    prunedDuringThisIter.clear();
//    prunedDuringThisIter.shrink_to_fit();
//}
