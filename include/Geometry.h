#ifndef Geometry_h
#define Geometry_h

#include <array>
#include <fstream>
#include <string>

class Geometry {
 public:
  // CB not optimal, but readable
  static constexpr std::array<int, 5> WHEELS{-2, -1, 0, 1, 2};
  static constexpr std::array<int, 12> SECTORS{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  static constexpr std::array<int, 4> STATIONS{1, 2, 3, 4};
  static constexpr double CELL_WIDTH{4.2};

  double xFirstWire(int wh, int st, int sec, int sl, int la) const {
    return geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][1];
  }

  double firstWire(int wh, int st, int sec, int sl, int la) const {
    return geometry[wh + 2][st - 1][sec - 1][sl - 1][la - 1][0];
  }

  static int sectorInRange(int sec) {
    if (sec == 13)
      sec = 4;
    else if (sec == 14)
      sec = 10;
    return sec;
  }

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

 private:
  double geometry[5][4][14][3][4][2];  // CB sectors are 14
};

const Geometry GEOM;

#endif
