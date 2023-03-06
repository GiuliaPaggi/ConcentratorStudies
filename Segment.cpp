#include "Segment.h"

Segment::Segment(int i, int st, int wh, int sec, int nHits, double x) /*: station(st), wheel(wh), sector(sec), nPhiHits(nHits), xLoc(x)*/ {
    index = i;
    station = st;
    wheel = wh;
    if (sec == 13) sector = 4;
    else if (sec == 14) sector = 10;
    else  sector = sec;
    nPhiHits = nHits;
    xLoc = x;
}

/*bool Segment::MatchCluster(Cluster cluster, double xCut){
    if ( std::abs(cluster.bestTP().xLoc - xLoc) < xCut) return true;
    else return false;
}*/