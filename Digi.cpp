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


// dovrebbero essere uguali (1, 2, 3, 5, 6, 7, 8, 12) (9, 11) (13, 4) (10, 14)
