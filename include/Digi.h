#ifndef Digi_h
#define Digi_h

class Digi{
    public:
    Digi(){};
    Digi(int i, int wh, int sec, int stat, int SL, int L, int w, double t);

    std::vector<Digi> FindCluster(std::vector<Digi> digisToCluster, double cut);


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

#endif
    