#ifndef TriggerPrimitive_h
#define TriggerPrimitive_h


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <vector>

class TriggerPrimitive{
    
    public:

    //TP information
    short int index = -1;
    short int wheel = -5;
    short int sector = -1;
    short int station = -1;
    short int quality = -1;
    double phi = -1000.;
    double phiB = -1000.;
    double psi = phi + phiB;
    int BX = -1000;
    int t0 = -1;
    float xLoc = -1000; //in cm
    double phiExpected[4];
    bool computedPhi = false;
    vector<int> Matches;
    bool hasMatched = false;

    bool hasRightBX = false;
    bool isGhostOutOfTime = false;

    // ADD LIKE ISONTIMEHHOST ISOFFTIMEGHOST
    

    //Cluster info
    bool inCluster = false;
    bool hasHigherQuality = false;

    TriggerPrimitive(){};
    TriggerPrimitive( short i, short ph2TpgPhiEmuAm_wheel, short ph2TpgPhiEmuAm_sector, short ph2TpgPhiEmuAm_station, short ph2TpgPhiEmuAm_quality, 
       int ph2TpgPhiEmuAm_phi, int ph2TpgPhiEmuAm_phiB, int ph2TpgPhiEmuAm_BX, int ph2TpgPhiEmuAm_t0, float ph2TpgPhiEmuAm_posLoc_x);
    ~TriggerPrimitive(){};

    void ComputeExpectedPhi();
    bool Match(TriggerPrimitive &TP, double PhiCut, double TimeCut);
    vector<int> MakeCluster(TriggerPrimitive listOfPrimitives[], int size,  double phicut); // taglio in cm, valore ragionevole 5 cm
    void FindHigherQuality (TriggerPrimitive listOfPrimitives[], vector<int> clusterIndeces);
    //vector<int> SelectRightBX(TriggerPrimitive listOfPrimitives[], vector<int> clusterIndeces);
    void CheckBX(); 

};


#endif

// classe cluster-> qualità massima al bx giusto, n di ghost al bx giusto, n ghost al bx sbagliato, accedere a tutti (vettore indici di trigger primitive)
// parti da prompt