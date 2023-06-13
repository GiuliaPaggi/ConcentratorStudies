#ifndef Geometry_h
#define Geometry_h

#include <fstream>
#include <string>
#include <vector>

class Geometry {

  double geometry[5][4][14][3][4][2]; // CB sectors are 14

public:
  // CB not optimal, but readable
  const std::vector<int> WHEELS{-2, -1, 0, 1, 2};
  const std::vector<int> SECTORS{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  const std::vector<int> STATIONS{1, 2, 3, 4};
  static constexpr double CELL_WIDTH{4.2};

  Geometry(std::string geoFileName = "geometry/DTGeom.txt") {
    std::ifstream geoFile;
    geoFile.open(geoFileName);

    std::string line;
    int wh, st, sec, sl, la;
    double firstW, xCoord;

    if (geoFile.is_open()) {
      while (getline(geoFile, line)) {
        while (geoFile >> wh >> st >> sec >> sl >> la >> firstW >> xCoord) {
          geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][0] = firstW;
          geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][1] = xCoord;
        }
      }
    }
    geoFile.close();
  }

  double xFirstWire(int wh, int st, int sec, int sl, int la) const {
    return geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][1];
  };

  double firstWire(int wh, int st, int sec, int sl, int la) const {
    return geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][0];
  };
};

#endif
