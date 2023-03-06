#ifndef Segment_h
#define Segment_h

#include "TMath.h"

class Segment{
    public:
    Segment(){};
    Segment(int i, int st, int wh, int sec, int nHits, double x );
    //bool MatchCluster(Cluster cluster, double xCut){};
    
    //private:
    std::size_t index{9999};
    double xLoc{-1000};
    int wheel{-5};
    int sector{-1};
    int station{-1};
    int nPhiHits{-1};

};

#endif  