#include "include/Segment.h"

#include "include/Geometry.h"

Segment::Segment(int i, int st, int wh, int sec, int nHits, double x)
    : index(i), station(st), wheel(wh), sector(GEOM.sectorInRange(sec)), nPhiHits(nHits), xLoc(x) {}
