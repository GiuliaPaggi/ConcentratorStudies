#ifndef Cluster_h
#define Cluster_h

#include "TriggerPrimitive.h"

class Cluster{
public:
  int wheel{-5};
  int station{-1};
  int sector{-1};

 private:
  TriggerPrimitive _bestTP{};

  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;

 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive> const& tps, double x_cut, int st, int wh, int sec);


  int itSize() const;
  int ootSize() const;

  bool isolated() const{ return bestTPQuality() < 0 && !hasGhosts(); };
  bool hasGhosts() const { return itSize() || ootSize(); };

  const std::vector<TriggerPrimitive>& itGhosts() const;
  const std::vector<TriggerPrimitive>& ootGhosts() const;

  int itCountIf(std::function<bool(TriggerPrimitive const&)> f) const;
  int ootCountIf(std::function<bool(TriggerPrimitive const&)>  f) const;

  int bestTPIndex() const;
  int bestTPQuality() const;
  const TriggerPrimitive& bestTP() const;

    //void MakeClusters(vector <Cluster> Clusters, TriggerPrimitive TPS[], int sz, double xCut);

};

#endif