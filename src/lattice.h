#pragma once

#include <queue>
#include "fractal.h"

const std::map<std::string, unsigned char> point_type {
    {"unknown", 0},
    {"outside", 1},
    {"boundary", 2},
    {"inside", 3}
};

class Lattice {
private:
    // matrix holding lattice points
    std::vector<std::vector<unsigned char>> r;

    // lattice spacing
    double d; 

    // original L
    double L;

    // number of grid points between two corners in fractal
    int scale; 

    // ray casting helper function
    unsigned char ray_casting(int x, int y, Fractal& fractal);

    // BFS helper functions
    void fill_boundary(Fractal& fractal);
    void fill_inside();

    // helmholtz helper functions
    void assert_lattice_points_type();
    std::vector<std::pair<int, int>> get_inside_indices();
    void write_arma_matrix_to_file(
            arma::sp_mat& A, std::string filename);
    void write_arma_complex_matrix_to_file(
            arma::cx_mat& A, std::string filename);
    void write_arma_complex_vector_to_file(
            arma::cx_vec& V, std::string filename);

    // helmholtz five-point helper functions
    arma::sp_mat helmholtz_matrix_five_point(
            std::vector<std::pair<int, int>>& inside_indices);

    // helmholtz thirteen-point helper functions
    arma::sp_mat helmholtz_matrix_thirteen_point(
            std::vector<std::pair<int, int>>& inside_indices);

    // biharmonic nine-point helper functions
    arma::sp_mat biharmonic_matrix_nine_point(
            std::vector<std::pair<int, int>>& inside_indices);

public:
    // constructor
    Lattice(Fractal& fractal, int scale);

    // point inside/outside boundary algorithms
    void ray_casting_all_points(Fractal& fractal);
    unsigned char BFS(int x, int y, Fractal& fractal);

    // helmholtz equation solver with five-point stencil
    void helmholtz_five_point();
    
    // helmholtz equation solver with thirteen-point stencil
    void helmholtz_thirteen_point();

    // biharmonic equation solver with nine-point stencil
    void biharmonic_nine_point();
    
    // helper functions
    void write_to_file(std::string filename);
    void reset_lattice_points();
    void write_data_to_file(std::string filename);
};
