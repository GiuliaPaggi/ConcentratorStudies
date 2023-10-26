#ifndef Digi_h
#define Digi_h

#include <vector>

struct Digi {
 public:
  Digi() = delete;
  Digi(int i, int wh, int sec, int st, int sl, int la, int w, double t);

  std::vector<Digi> findCluster(const std::vector<Digi> & digis, double cut) const;

  std::size_t index{9999};
  int wheel{-1};
  int sector{-1};
  int station{-1};
  int superlayer{-1};
  int layer{-1};
  int wire{-1};
  double time{-1.0};
  double xLoc{-1.0};

  bool inCluster{false};
};

inline bool operator==(Digi const &ldigi, Digi const &rdigi) { return ldigi.index == rdigi.index; }

#endif
