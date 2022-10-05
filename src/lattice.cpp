#include "lattice.h"

Lattice::Lattice(Fractal& fractal, int scale) : d{fractal.L / (scale+1)}, scale{scale} {
    // calculate original length
    L = fractal.L * pow(4.0, fractal.l);
    
    // find number of points needed in lattice
    double extent = L / 2.0;
    for (int i = 1; i <= fractal.l; i++) {
        extent += L / pow(4.0, i);
    }
    int n = extent * 2.0 / d + 3;
    r.resize(n);

    // initialize "outside" around fractal and unknown inside
    for (int i = 0; i < n; i++) {
        r[0].push_back(point_type.at("outside"));
        r[n-1].push_back(point_type.at("outside"));
        if (i > 0 && i < n-1) {
            r[i].push_back(point_type.at("outside"));
            for (int j = 0; j < n - 2; j++) {
                r[i].push_back(point_type.at("unknown"));
            }
            r[i].push_back(point_type.at("outside"));
        }
    }
}

// ray casting functions
void Lattice::ray_casting_all_points(Fractal& fractal) {
    // ray cast all points in lattice
    // the border is already set to outside in the constructor
    for (int i = 1; i < r.size() - 1; i++) {
        for (int j = 1; j < r.size() - 1; j++) {
            r[i][j] = ray_casting(i, j, fractal);
        }
    }
}

unsigned char Lattice::ray_casting(int x, int y, Fractal& fractal) {
    // cast ray to the right and see how many times it crosses the boundary
    double n2 = r.size() / 2.0;
    int crossings = 0;

    // go through all fractal points to calculate number of crossings
    for (int i = 0; i < fractal.r.size(); i++) {
        Point p1 = fractal.r[i]; Point p2 = fractal.r[(i+1)%fractal.r.size()];
        
        
        int x_min = std::min(p1.x, p2.x)/d + n2; int x_max = std::max(p1.x, p2.x)/d + n2;
        int y_min = std::min(p1.y, p2.y)/d + n2; int y_max = std::max(p1.y, p2.y)/d + n2;
        
        if ((x >= x_min && x <= x_max) && (y >= y_min && y <= y_max)) {
            return point_type.at("boundary");
        }
        
        // special case if boundary is horizontal
        if ((y >= y_min && y < y_max) && (x <= x_min)) {
            if (y_min == y_max) {
                int index = i - 1;
                if (index == -1) {
                    index = fractal.r.size() - 1;
                }
                Point p0 = fractal.r[index];
                Point p3 = fractal.r[(i+2)%fractal.r.size()];
                if (p0.y < p1.y && p3.y > p2.y) {
                    crossings++;
                }
                else if (p0.y > p1.y && p3.y < p2.y) {
                    crossings++;
                }
            }
            else {
                crossings++;
            }
        }
    }
    
    if (crossings % 2 == 0) {
        return point_type.at("outside");
    }
    return point_type.at("inside");
}
// end of ray casting functions

// BFS functions
unsigned char Lattice::BFS(int x, int y, Fractal& fractal) {
    // initialize all outside, fill boundary, do BFS from known point inside
    for (int i = 1; i < r.size() - 1; i++) {
        for (int j = 1; j < r.size() - 1; j++) {
            r[i][j] = point_type.at("outside");
        }
    }
    fill_boundary(fractal); // fill boundary from fractal
    fill_inside(); // use BFS to fill inside
    return r[x][y];
}

void Lattice::fill_boundary(Fractal& fractal) {
    double n = r.size(); // need double division
    for (int i = 0; i < fractal.r.size(); i++) {
        Point p1 = fractal.r[i]; 
        Point p2 = fractal.r[(i+1)%fractal.r.size()];
        int l_x1 = p1.x/d + n/2.0; int l_x2 = p2.x/d + n/2.0;
        int l_y1 = p1.y/d + n/2.0; int l_y2 = p2.y/d + n/2.0;
        int diffx, diffy;
        if (l_x1 != l_x2) {
            diffx = (l_x2 - l_x1) / abs(l_x1 - l_x2);
        } else {diffx = 0;}
        if (l_y1 != l_y2) {
            diffy = (l_y2 - l_y1) / abs(l_y1 - l_y2);
        } else {diffy = 0;}
        for (int j = 0; j < scale + 1; j++) {
            r[l_x1 + j*diffx][l_y1 + j*diffy]= point_type.at("boundary");
        }
    }
}

void Lattice::fill_inside() {
    int n = r.size();
    std::queue<std::pair<int, int>> queue;
    queue.push({n/2, n/2});
    std::vector<std::pair<int, int>> indices{
        {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
    std::vector<std::pair<int, int>> visited;
    while (!queue.empty()) {
        int x = queue.front().first; int y = queue.front().second;
        for (int i = 0; i < indices.size(); i++) {
            int a = x + indices[i].first;
            int b = y + indices[i].second;
            bool found = false;
            for (int j = 0; j < visited.size(); j++) {
                if (std::pair<int, int>{a, b} == visited[j]) {
                    found = true;
                }
            }
            if (!found && r[a][b] == point_type.at("outside")) {
                queue.push({a, b});
            }
            visited.push_back({a, b});
        }
        if (r[x][y] != point_type.at("boundary")) {
            r[x][y] = point_type.at("inside");
        }
        queue.pop();
    }
}
// end of BFS functions

// helmholtz five-point functions
void Lattice::helmholtz_five_point() {
    
    // check if all lattice points have a type != "unknown"
    assert_lattice_points_type();
    
    // indices of points inside boundary
    auto inside_indices = get_inside_indices();
    
    // five-point finite difference matrix
    auto A = helmholtz_matrix_five_point(inside_indices);
    
    // get 10 smallest eigenvalues and eigenvectors
    arma::cx_vec eigval; arma::cx_mat eigvec;
    arma::eigs_gen(eigval, eigvec, A, 10, "sm");
    
    // correct eigenvalues
    eigval *= pow(d, 2);

    //write_arma_matrix_to_file(A, "../data/five_point_matrix.txt");
    write_arma_complex_matrix_to_file(eigvec, "../data/helmholtz_five_point_eigvecs.txt");
    write_arma_complex_vector_to_file(eigval, "../data/helmholtz_five_point_eigvals.txt");
}

void Lattice::assert_lattice_points_type() {
    for (int i = 0; i < r.size(); i++) {
        for (int j = 0; j < r.size(); j++) {
            assert(r[i][j] != point_type.at("unknown"));
        }
    }
}

std::vector<std::pair<int, int>> Lattice::get_inside_indices() {
    std::vector<std::pair<int, int>> inside_indices;
    for (int col = 1; col < r.size() - 1; col++) {
        for (int row = 1; row < r.size() - 1; row++) {
            if (r[row][col] == point_type.at("inside")) {
                inside_indices.push_back({row, col});
            }
        }
    }
    return inside_indices;
}

arma::sp_mat Lattice::helmholtz_matrix_five_point(
        std::vector<std::pair<int, int>>& inside_indices) {
    
    // create matrix
    auto A = arma::sp_mat(inside_indices.size(), inside_indices.size());
    
    // relative indices to be iterated over ((x+1, y+0), (x-1, y+0), (x+0, y+1), (x+0, y-1))
    const std::vector<std::pair<int, int>> indices{
        {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    const std::vector<double> constants{
        -1.0, -1.0, -1.0, -1.0};
    const double diag = 4.0;
    
    // go through all points inside and set matrix elements according to five-point stencil
    for (int i = 0; i < inside_indices.size(); i++) {
        // set diagonal
        A(i, i) = diag;
        int x = inside_indices[i].first; int y = inside_indices[i].second;
        // go through points around (x, y)
        for (int j = 0; j < indices.size(); j++) {
            int a = indices[j].first; int b = indices[j].second;

            // if the point next to (x, y) is inside we need to change a matrix element
            if (r[x+a][y+b] == point_type.at("inside")) {
                // go through inside_indices and find which column (x+a, y+b) represent
                int col = -1;
                for (int k = 0; k < inside_indices.size(); k++) {
                    if (std::pair{x+a, y+b} == inside_indices[k]) {
                        col = k;
                        break;
                    }
                }
                A(i, col) = constants[j];
            }
        }
    }
    return A;
}
// end of helmholtz five-point functions

// helmholtz thirteen-point functions
void Lattice::helmholtz_thirteen_point() {
    assert_lattice_points_type();

    auto inside_indices = get_inside_indices();
    auto A = helmholtz_matrix_thirteen_point(inside_indices);
    arma::cx_vec eigval; arma::cx_mat eigvec;
    arma::eigs_gen(eigval, eigvec, A, 10, "sm");
    eigval *= pow(d, 2);

    // write_arma_matrix_to_file(A, "../data/thirteen_point_matrix.txt");
    write_arma_complex_matrix_to_file(eigvec, "../data/helmholtz_thirteen_point_eigvecs.txt");
    write_arma_complex_vector_to_file(eigval, "../data/helmholtz_thirteen_point_eigvals.txt");
}


arma::sp_mat Lattice::helmholtz_matrix_thirteen_point (
        std::vector<std::pair<int, int>>& inside_indices) {

    auto A = arma::sp_mat(inside_indices.size(), inside_indices.size());
    
    const std::vector<std::pair<int, int>> thirteen_point_indices {
        {3, 0}, {2, 0}, {1, 0}, {-1, 0}, {-2, 0}, {-3, 0}, 
        {0, 3}, {0, 2}, {0, 1}, {0, -1}, {0, -2}, {0, -3}};
    const std::vector<double> thirteen_point_constants {
        -1.0/90.0, 3.0/20.0, -3.0/2.0, -3.0/2.0, 3.0/20.0, -1.0/90.0,
        -1.0/90.0, 3.0/20.0, -3.0/2.0, -3.0/2.0, 3.0/20.0, -1.0/90.0};
    double thirteen_point_diag = 49.0/18.0 * 2.0;
    
    const std::vector<std::pair<int, int>> five_point_indices{
        {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    const std::vector<double> five_point_constants {
        -1.0, -1.0, -1.0, -1.0};
    double five_point_diag = 4.0;

    double diag;
    std::vector<std::pair<int, int>> indices;
    std::vector<double> constants;
    for (int i = 0; i < inside_indices.size(); i++) {
        int x = inside_indices[i].first; int y = inside_indices[i].second;
        
        // if point is in the third, second or nearest layer to boundary
        // we need to perform a five-point finite difference scheme instead of
        // thirteen-point because boundary conditions is not defined outside boundary
        if (x >= r.size() - 3 || x <= 2 || y >= r.size() - 3 || y <= 2) {
            diag = five_point_diag;
            indices = five_point_indices;
            constants = five_point_constants;
        }
        else {
            diag = thirteen_point_diag;
            indices = thirteen_point_indices;
            constants = thirteen_point_constants;
        }
        
        A(i, i) = diag;
        for (int j = 0; j < indices.size(); j++) {
            int a = indices[j].first; int b = indices[j].second;
            if (r[x+a][y+b] == point_type.at("inside")) {
                int col = -1;
                for (int k = 0; k < inside_indices.size(); k++) {
                    if (std::pair<int, int>{x+a, y+b} == inside_indices[k]) {
                        col = k;
                        break;
                    }
                }
                A(i, col) = constants[j];
            }
        }
    }
    return A;
}
// end of helmholtz thirteen-point functions

// biharmonic nine-point functions
arma::sp_mat Lattice::biharmonic_matrix_nine_point(
        std::vector<std::pair<int, int>>& inside_indices) {

    auto A = arma::sp_mat(inside_indices.size(), inside_indices.size());
    
    const std::vector<std::pair<int, int>> indices {
        {2, 0}, {1, 0}, {-1, 0}, {-2, 0}, 
        {0, 2}, {0, 1}, {0, -1}, {0, -2}};
    const std::vector<double> constants {
        1.0, -4.0, -4.0, 1.0,
        1.0, -4.0, -4.0, 1.0};
    double diag = 12.0;
    
    for (int i = 0; i < inside_indices.size(); i++) {
        A(i, i) = diag;
        int x = inside_indices[i].first; int y = inside_indices[i].second;
        for (int j = 0; j < indices.size(); j++) {
            int a = indices[j].first; int b = indices[j].second;
            
            if (r[x+a][y+b] == point_type.at("inside")) {
                int col = -1;
                for (int k = 0; k < inside_indices.size(); k++) {
                    if (std::pair{x+a, y+b} == inside_indices[k]) {
                        col = k;
                        break;
                    }
                }
                A(i, col) = constants[j];
            }
            if (r[x+a][y+b] == point_type.at("outside")) {
                A(i, i) += 1;
            }
        }
    }
    return A;
}

void Lattice::biharmonic_nine_point() {
    assert_lattice_points_type();

    auto inside_indices = get_inside_indices();
    auto A = biharmonic_matrix_nine_point(inside_indices);
    
    arma::cx_vec eigval; arma::cx_mat eigvec;
    arma::eigs_gen(eigval, eigvec, A, 10, "sm");
    eigval *= pow(d, 2);

    // write_arma_matrix_to_file(A, "../data/nine_point_matrix.txt");
    write_arma_complex_matrix_to_file(eigvec, "../data/biharmonic_nine_point_eigvecs.txt");
    write_arma_complex_vector_to_file(eigval, "../data/biharmonic_nine_point_eigvals.txt");
}
// end of biharmonic nine-points functions

// helper functions
void Lattice::write_to_file(std::string filename) {
    std::ofstream out;
    out.open(filename);
    for (int i = 0; i < r.size(); i++) {
        for (int j = 0; j < r.size(); j++) {
            out << static_cast<int>(r[i][j]) << " ";
        }
        out << "\n";
    }
    out.close();
}

void Lattice::write_data_to_file(std::string filename) {
    std::ofstream out; out.open(filename);
    out << L << " " << d << " " << scale;
    out.close();
}

void Lattice::reset_lattice_points() {
   for (int i = 1; i < r.size() - 1; i++) {
      for (int j = 1; j < r.size() - 1; j++) {
         r[i][j] = point_type.at("unknown");
      }
   }
}

void Lattice::write_arma_matrix_to_file(arma::sp_mat& A, std::string filename) {
    std::ofstream out; out.open(filename);
    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
            out << A(i, j) << " ";
        }
        out << "\n";
    }
    out.close();
}

void Lattice::write_arma_complex_matrix_to_file(arma::cx_mat& A, std::string filename) {
    std::ofstream out; out.open(filename);
    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
           out << A(i, j).real() << " "; 
        }
        out << "\n";
    }
    out.close();
}

void Lattice::write_arma_complex_vector_to_file(arma::cx_vec& V, std::string filename) {
    std::ofstream out; out.open(filename);
    for (int i = 0; i < V.n_elem; i++) {
        out << V(i).real() << "\n";
    }
    out.close();
}
// end of helper functions
