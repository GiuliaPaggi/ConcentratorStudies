#ifndef Cluster_h
#define Cluster_h

#include "TriggerPrimitive.h"

class Cluster{
    private:
    vector <TriggerPrimitive> _OutofTimeGhosts;
    vector <TriggerPrimitive> _InTimeGhosts;
    TriggerPrimitive _BestQuality;
    int wheel = -5;
    int station = -1;
    int sector = -1;

    public: 
    Cluster(){};
    Cluster(TriggerPrimitive Tps[], vector<int> ClusterIndex);
    Cluster(TriggerPrimitive Tps[], int size, double xcut, int st, int wh, int sec);
    int GetOoTSize();
    int GetITSize();
    int GetBestQualityIndex();
    vector<double> GetOoTBX();
    vector<double> GetITBX();
    vector<double> GetOoTResidual();
    vector<double> GetITResidual();
    vector<int> GetOoTQualities();
    vector<int> GetITQualities();
    int GetBestQuality();

    //void MakeClusters(vector <Cluster> Clusters, TriggerPrimitive TPS[], int sz, double xCut);

};

#endif