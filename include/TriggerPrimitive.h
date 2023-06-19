#ifndef TriggerPrimitive_h
#define TriggerPrimitive_h

#include <array>
#include <vector>

class TriggerPrimitive {
 public:
  // TP information
  std::size_t index{9999};
  int wheel{-5};
  int sector{-1};
  int station{-1};
  int quality{-1};

  int BX{-1000};
  int t0{-100000};  // in ns

  float xLoc{-1000.0};  // in cm

  double phi{-1000.0};   // in rad
  double phiB{-1000.0};  // in rad
  double psi{-1000.0};   // in rad

  std::array<double, 4> phiExpected;
  std::vector<int> matches;

  bool computedPhi{false};
  bool hasMatched{false};

  bool hasRightBX{false};
  bool isGhostOutOfTime{false};

  // Cluster info
  bool inCluster{false};
  bool hasHighestQuality{false};

  TriggerPrimitive() = default;
  TriggerPrimitive(std::size_t i, int tpg_wheel, int tpg_sector, int tpg_station, int tpg_quality, int tpg_phi,
                   int tpg_phiB, int tpg_BX, int tpg_t0, float tpg_posLoc_x);

  ~TriggerPrimitive(){};

  void computeExpectedPhi();
  bool match(TriggerPrimitive &tp, double phiCut, double timeCut);
  // std::vector<int> makeCluster(TriggerPrimitive listOfPrimitives[], int size,
  //                         double phiCut);  // taglio in cm, valore ragionevole 5 cm
  void findHigherQuality(TriggerPrimitive listOfPrimitives[], const std::vector<int> &clusterIndices);
  // std::vector<int> selectRightBX(TriggerPrimitive listOfPrimitives[], const std::vector<int>& clusterIndices);
  void checkBX();

  // Filter functions
  bool MatchFromLQ(TriggerPrimitive &TP, double PhiCut, double TimeCut);
  bool MatchFromHQ(TriggerPrimitive &TP, double PhiCut, double TimeCut);
};

inline bool operator==(TriggerPrimitive const &ltp, TriggerPrimitive const &rtp) { return ltp.index == rtp.index; }

#endif
