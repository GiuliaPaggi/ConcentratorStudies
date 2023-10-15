#ifndef Segment_h
#define Segment_h

#include <vector>
struct Segment {
 public:
  Segment() = default;
  Segment(int i, int st, int wh, int sec, int nHits, double x);

  std::size_t index{9999};
  int station{-1};
  int wheel{-5};
  int sector{-1};
  int nPhiHits{-1};
  double xLoc{-5000};

  bool inCluster{false};
};

inline bool operator==(Segment const &lseg, Segment const &rseg) { return lseg.index == rseg.index; }

#endif