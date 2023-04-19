#include "Digi.h"
#include <iostream>
#include <fstream>
#include <string>

const double cell = 4.2;
bool readfile = false;
double geometry[5][4][12][3][4][2];


void ReadGeometry(){
    std::ifstream geofile;
    geofile.open("DTGeom.txt");
    std::string line;
    int WH, ST, SEC, SL, LA, FIRSTW;
    double XCOORD;

    if (geofile.is_open() ){
        while (getline(geofile ,line)) {
            while(geofile >> WH >> ST >> SEC >> SL >> LA >> FIRSTW >> XCOORD){
                geometry[WH+2][ST-1][SEC-1][SL-1][LA-1][0] = FIRSTW;
                geometry[WH+2][ST-1][SEC-1][SL-1][LA-1][1] = XCOORD;
            }
        }
    }
    geofile.close();
}


Digi::Digi(int i, int wh, int sec, int stat, int SL, int L, int w, double t) :  wheel{wh}, sector{sec}, station{stat}, superlayer{SL}, layer{L}, wire{w}, time{t} {
    index = i;

    if (sec == 13) sector = 4;
    else if (sec == 14) sector = 10;

    if (!readfile) {
        ReadGeometry();
        readfile = true;
    }
    
    int firstWire = geometry[wheel+2][station-1][sector-1][superlayer-1][layer-1][0];
    double xFirstWire = geometry[wheel+2][station-1][sector-1][superlayer-1][layer-1][1]; 

    xLoc = xFirstWire + (wire - firstWire) *cell;
}

std::vector<Digi> Digi::FindCluster(std::vector<Digi> digisToCluster, double cut){
    std::vector<Digi> DigiCluster;

    // find digis close to the one that calls the function
    std::copy_if( digisToCluster.begin(), digisToCluster.end(), std::back_inserter(DigiCluster),
                    [=](auto& d) {return std::abs(xLoc -d.xLoc) < cut;});
    
    // if there are more digis in the chamber-> find mean value of xLoc in the cluster we have for now and use that to look for close digis
    if (digisToCluster.size() > DigiCluster.size()){
        double mean_xLoc = xLoc;
        for (auto d : digisToCluster) mean_xLoc+=d.xLoc;
        mean_xLoc = mean_xLoc/DigiCluster.size();
        //clear cluster to fill it again with new center
        DigiCluster.clear();
        std::copy_if( digisToCluster.begin(), digisToCluster.end(), std::back_inserter(DigiCluster),
                    [=](auto& d) {return std::abs(mean_xLoc -d.xLoc) < cut;});
    }

    return DigiCluster;

}