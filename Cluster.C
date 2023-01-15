#include "Cluster.h"

int RightBX = -380;

Cluster::Cluster(TriggerPrimitive Tps[], vector<int> ClusterIndex){
    //Takes a list of TP and the indeces of the ones to put in the cluster
    for (int i : ClusterIndex) {
        if (Tps[i].hasHigherQuality) {
            _BestQuality = Tps[i];
        }
        else {
            if (Tps[i].BX != RightBX) _OutofTimeGhosts.push_back(Tps[i]);
            else _InTimeGhosts.push_back(Tps[i]);
        }
    }
};

Cluster::Cluster(TriggerPrimitive Tps[], int size, double xcut, int st, int wh, int sec){
    // in a wheel-station-sector finds the higher quality TP and looks for LQ ones in a xcut intervall
    wheel = wh;
    station = st;
    sector = sec;
    //find best quality in wh st
    int maxQ = 0;
    for (int i = 0; i < size; ++i ) {
        TriggerPrimitive element = Tps[i];
        if (element.wheel != wh) continue;
        if (element.station != st) continue;
        if (element.sector != sec) continue;
        if (abs(element.BX - RightBX) > 1 ) continue;
        if (element.quality > maxQ) {
            maxQ = element.quality;
            _BestQuality = element;
        }
    }
    // find cluster around best quality
    double DeltaxLoc;
    for (int i = 0; i < size; ++i ) {
        TriggerPrimitive element = Tps[i];
        if (element.index == _BestQuality.index) continue;
        if (element.wheel != wh) continue;
        if (element.station != st) continue;
        if (element.sector != sector) continue;
        DeltaxLoc = abs(element.xLoc - _BestQuality.xLoc);
    // cout << DeltaxLoc << "      " << element.wheel << "     " << element.station << "       "<< wh << "     " <<st  << endl;
        if (DeltaxLoc < xcut) {
            if ( abs(element.BX - RightBX) > 1) _OutofTimeGhosts.push_back(element);
            else _InTimeGhosts.push_back(element);
        }
    }
    
    //cout << _OutofTimeGhosts.size() << "            ||          " << _InTimeGhosts.size() << endl;
};


int Cluster::GetOoTSize(){
    return _OutofTimeGhosts.size();
};

int Cluster::GetITSize(){
    return _InTimeGhosts.size();
};

int Cluster::GetBestQualityIndex(){
    return _BestQuality.index;
};

vector<double> Cluster::GetOoTBX(){
    vector<double> ootBX;
    for (auto i : _OutofTimeGhosts){
        ootBX.push_back(i.BX);
    }
    return ootBX;
};

vector<double> Cluster::GetITBX(){
    vector<double> itBX;
    for (auto i : _InTimeGhosts){
        itBX.push_back(i.BX);
    }
    return itBX;
};

vector<double> Cluster::GetOoTResidual(){
    vector<double> ootResidual;
    for (auto i : _OutofTimeGhosts){
        ootResidual.push_back(i.xLoc - _BestQuality.xLoc);
    }
    return ootResidual;
};

vector<double> Cluster::GetITResidual(){
    vector<double> itResidual;
    for (auto i : _InTimeGhosts){
        itResidual.push_back(i.xLoc- _BestQuality.xLoc);
    }
    return itResidual;
};

vector<int> Cluster::GetOoTQualities(){
    vector<int> ootQualities;
    for (auto i : _OutofTimeGhosts){
        ootQualities.push_back(i.quality);
    }
    return ootQualities;
};

vector<int> Cluster::GetITQualities(){
    vector<int> itQualities;
    for (auto i : _InTimeGhosts){
        itQualities.push_back(i.quality);
    }
    return itQualities;
};

int Cluster::GetBestQuality(){
    return _BestQuality.quality;
};


