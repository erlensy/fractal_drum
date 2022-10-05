#include <thread>

#include "lattice.h"

void generate_fractals() {
    double L = 5.0; 

    Fractal fractal0{L, 0};
    fractal0.write_to_file("../data/fractal_level_0.txt");
        
    Fractal fractal1{L, 1};
    fractal1.write_to_file("../data/fractal_level_1.txt");

    Fractal fractal2{L, 2};
    fractal2.write_to_file("../data/fractal_level_2.txt");

    Fractal fractal3{L, 3};
    fractal3.write_to_file("../data/fractal_level_3.txt");

    Fractal fractal4{L, 4};
    fractal4.write_to_file("../data/fractal_level_4.txt");

    Fractal fractal5{L, 5};
    fractal5.write_to_file("../data/fractal_level_5.txt");
}

void test_point_inside_outside() {
    double L = 5.0;
    Fractal fractal{L, 3};
    Lattice lattice{fractal, 0};

    lattice.ray_casting_all_points(fractal);
    lattice.write_to_file("../data/lattice_ray_casting.txt");
    
    lattice.reset_lattice_points();

    lattice.BFS(0, 0, fractal);
    lattice.write_to_file("../data/lattice_BFS.txt");
}

void helmholtz_five_point() {
    double L = 5.0; int l = 3;

    Fractal fractal(L, l);
    Lattice lattice{fractal, 0};

    lattice.BFS(0, 0, fractal);
    
    lattice.helmholtz_five_point();
    lattice.write_to_file("../data/helmholtz_five_point_lattice.txt");
    lattice.write_data_to_file("../data/helmholtz_five_point_lattice_data.txt");
}

void helmholtz_thirteen_point() {
    double L = 5.0; int l = 3;

    Fractal fractal(L, l);
    Lattice lattice{fractal, 0};

    lattice.BFS(0, 0, fractal);
    
    lattice.helmholtz_thirteen_point();
    lattice.write_to_file("../data/helmholtz_thirteen_point_lattice.txt");
    lattice.write_data_to_file("../data/helmholtz_thirteen_point_lattice_data.txt");
}

void biharmonic_nine_point() {
    double L = 5.0; int l = 3;

    Fractal fractal(L, l);
    Lattice lattice{fractal, 0};

    lattice.BFS(0, 0, fractal);
    
    lattice.biharmonic_nine_point();
    lattice.write_to_file("../data/biharmonic_nine_point_lattice.txt");
    lattice.write_data_to_file("../data/biharmonic_nine_point_lattice_data.txt");
}

int main() {

    std::vector<std::thread> threads;
    
    threads.push_back(std::thread(
            generate_fractals));
    
    threads.push_back(std::thread(
            test_point_inside_outside));

    threads.push_back(std::thread(
            helmholtz_five_point));

    threads.push_back(std::thread(
            helmholtz_thirteen_point));
    
    threads.push_back(std::thread(
            biharmonic_nine_point));
    
    for (auto& thread : threads) {
        thread.join();
    }

    return 0;
}
