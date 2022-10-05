import numpy as np
from multiprocessing import Process
from matplotlib import pyplot as plt

def plot_fractal(filename, lw):
    f = open(filename) # open file
    lines = f.readlines() # read lines
    f.close()
    r = np.zeros((len(lines), 2))
    for i, line in enumerate(lines):
        r[i, :] = list(map(float, line.strip().split(" ")))

    fig = plt.figure();
    plt.plot(r[:, 0], r[:, 1], color = "blue", linewidth = lw)
    plt.plot([r[-1, 0], r[0, 0]], [r[-1, 1], r[0, 1]], color = "blue", linewidth = lw)
    plt.savefig("../figures" + filename[7:-3] + "pdf", dpi = 600)
    plt.close(fig)

def generate_fractals():
    plot_fractal("../data/fractal_level_0.txt", 1.0)
    plot_fractal("../data/fractal_level_1.txt", 1.0)
    plot_fractal("../data/fractal_level_2.txt", 1.0)
    plot_fractal("../data/fractal_level_3.txt", 1.0)
    plot_fractal("../data/fractal_level_4.txt", 0.5)
    plot_fractal("../data/fractal_level_5.txt", 0.1)

def plot_lattice_type(filename):
    f = open(filename) # open file
    lines = f.readlines() # read lines
    f.close()
    r = np.zeros((len(lines), len(lines)))

    for i, line in enumerate(lines):
        r[i] = list(map(int, line.strip().split(" ")))

    fig = plt.figure()
    plt.imshow(r)
    plt.savefig("../figures" + filename[7:-3] + "pdf", dpi = 600)
    plt.close(fig)

def test_point_inside_outside_algorithms():
    plot_lattice_type("../data/lattice_BFS.txt")
    plot_lattice_type("../data/lattice_ray_casting.txt")

def plot_matrix(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    n = len(lines)
    r = np.zeros((n, n))
    for i, line in enumerate(lines):
        r[i] = list(map(float, line.strip().split(" ")))

    fig = plt.figure()
    plt.imshow(r)
    plt.savefig("../figures" + filename[7:-3] + "pdf", dpi = 600)

def plot_eigenmode(filename_eigvec, filename_lattice, filename_lattice_data, eigenmode):
    f = open(filename_lattice)
    lines = f.readlines()
    f.close()
    n = len(lines)
    r = np.zeros((n, n))
    lattice_values = np.zeros((n, n))
    for i, line in enumerate(lines):
        r[i] = list(map(int , line.strip().split(" ")))

    f = open(filename_eigvec)
    lines = f.readlines()
    f.close()
    eigvecs = np.zeros((len(lines), len(list(map(float, lines[0].strip().split(" "))))))
    for i, line in enumerate(lines):
        eigvecs[i] = list(map(float, line.strip().split(" ")))
    
    k = 0
    boundary = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if r[j, i] == 3:
                lattice_values[j, i] = eigvecs[k, eigenmode]
                k += 1        
            if r[j, i] == 2:
                boundary[j, i] = 1

    f = open(filename_lattice_data)
    L, d, scale = list(map(float, f.readlines()[0].strip().split(" ")))
    f.close()

    a = int(n/2) * d / L
    extent = (-a, a, -a, a)

    lattice_values = np.ma.masked_where(abs(lattice_values) < 1e-16, lattice_values)
                
    fig = plt.figure()
    plt.contourf(lattice_values, origin = "upper", cmap = "seismic", alpha = 0.9, extent = extent)
    plt.colorbar()
    plt.imshow(boundary, cmap = "Greys", origin = "upper", extent = extent)
    plt.savefig("../figures" + filename_eigvec[7:-5] + str(eigenmode) + ".pdf", dpi = 600)
    plt.close(fig)

def helmholtz_thirteen_point():
    for i in range(10):
        plot_eigenmode("../data/helmholtz_thirteen_point_eigvecs.txt", 
                       "../data/helmholtz_thirteen_point_lattice.txt", 
                       "../data/helmholtz_thirteen_point_lattice_data.txt", i)
    #plot_matrix("../data/helmholtz_thirteen_point_matrix.txt")

def helmholtz_five_point():
    for i in range(10):
        plot_eigenmode("../data/helmholtz_five_point_eigvecs.txt", 
                       "../data/helmholtz_five_point_lattice.txt", 
                       "../data/helmholtz_five_point_lattice_data.txt", i)
    #plot_matrix("../data/helmholtz_five_point_matrix.txt")

def biharmonic_nine_point():
    for i in range(10):
        plot_eigenmode("../data/biharmonic_nine_point_eigvecs.txt", 
                       "../data/biharmonic_nine_point_lattice.txt", 
                       "../data/biharmonic_nine_point_lattice_data.txt", i)

    #plot_matrix("../data/biharmonic_nine_point_matrix.txt")

if __name__ == "__main__":
    generate_fractals()

    #test_point_inside_outside_algorithms()

    helmholtz_five_point()

    helmholtz_thirteen_point()

    biharmonic_nine_point()

