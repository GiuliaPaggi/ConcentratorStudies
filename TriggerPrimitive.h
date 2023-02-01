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
  int t0{-100000};

  float xLoc{-1000.0};  // in cm

  double phi{-1000.0};   // in rad
  double phiB{-1000.0};  // in rad
  double psi{-1000.0};   // in rad

  std::array<double, 4> phiExpected;
  vector<int> Matches;

  bool computedPhi{false};
  bool hasMatched{false};

  bool hasRIGHT_BX{false};
  bool isGhostOutOfTime{false};

  // ADD LIKE ISONTIMEHHOST ISOFFTIMEGHOST

  // Cluster info
  bool inCluster{false};
  bool hasHigherQuality{false};

  TriggerPrimitive(){};
  TriggerPrimitive(std::size_t i, int tpg_wheel, int tpg_sector, int tpg_station,
                                   int tpg_quality, int tpg_phi, int tpg_phiB, int tpg_BX,
                                   int tpg_t0, float tpg_posLoc_x);

  ~TriggerPrimitive(){};

  void ComputeExpectedPhi();
  bool Match(TriggerPrimitive &TP, double PhiCut, double TimeCut);
  vector<int> MakeCluster(TriggerPrimitive listOfPrimitives[], int size,
                          double phicut);  // taglio in cm, valore ragionevole 5 cm
  void FindHigherQuality(TriggerPrimitive listOfPrimitives[], vector<int> clusterIndeces);
  vector<int> SelectRIGHT_BX(TriggerPrimitive listOfPrimitives[], vector<int> clusterIndeces);
  void CheckBX();
};

bool operator==(TriggerPrimitive const &ltp, TriggerPrimitive const &rtp) {
  return ltp.index == rtp.index;
}

#endif

// classe cluster-> qualità massima al bx giusto, n di ghost al bx giusto, n
// ghost al bx sbagliato, accedere a tutti (vettore indici di trigger primitive)
// parti da prompt