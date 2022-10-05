#pragma once

#include <armadillo>

// corner point on fractal (x, y)
struct Point {
    // constructor
    Point(double x, double y);
    
    // position
    double x;
    double y;
};

class Fractal {
public:
    // initialize square around (0, 0) with side lengths L
    Fractal(double L, int l);
    
    // length between two corners
    double L;
    
    // fractal level
    int l;
    
    // generate fractal level l for current fractal
    void generator(int l);
    
    // positions of corners in the fractal
    std::vector<Point> r;
    
    // helper functions
    void write_to_file(std::string filename);
    friend std::ostream& operator<<(std::ostream& os, Fractal& f);
};
