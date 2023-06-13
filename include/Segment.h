#ifndef Segment_h
#define Segment_h

#include "TMath.h"

class Segment{
    public:
    Segment() = default;
    Segment(int i, int st, int wh, int sec, int nHits, double x);
    //bool MatchCluster(Cluster cluster, double xCut){};
    
    //private:
    std::size_t index{9999};
    double xLoc{-5000};
    int wheel{-5};
    int sector{-1};
    int station{-1};
    int nPhiHits{-1};

    bool inCluster{false};

};

inline bool operator==(Segment const &lseg, Segment const &rseg) {
  return lseg.index == rseg.index;
}


#endif  