#include "fractal.h"

Point::Point(double x, double y) : x{x}, y{y} {}

Fractal::Fractal(double L, int l) : L{L}, l{l} {
    double L2 = L / 2.0;
    r.push_back(Point(-L2, -L2));
    r.push_back(Point(-L2, L2));
    r.push_back(Point(L2, L2));
    r.push_back(Point(L2, -L2));

    generator(l);
}

void Fractal::generator(int l) {
    for (int i = 0; i < l; i++) {
        std::vector<Point> temp;
        double L4;
        for (int j = 0; j < r.size(); j++) {
            Point p1 = r[j]; Point p9 = r[(j+1) % r.size()];
            L4 = L / 4.0;
            double dir = atan2(p9.y - p1.y, p9.x - p1.x);
            double a = L4 * cos(dir);
            double b = L4 * sin(dir);

            Point p2{p1.x + a, p1.y + b};
            Point p3{p2.x - b, p2.y + a};
            Point p4{p3.x + a, p3.y + b};
            Point p5{p4.x + b, p4.y - a};
            Point p6{p5.x + b, p5.y - a};
            Point p7{p6.x + a, p6.y + b};
            Point p8{p7.x - b, p7.y + a};
            temp.insert(temp.end(), {p1, p2, p3, p4, p5, p6, p7, p8});
        }
        r = temp;
        L = L4;
    }
}

void Fractal::write_to_file(std::string filename) {
    std::ofstream out; out.open(filename);
    for (auto p : r) {
        out << p.x << " " << p.y << "\n";
    }
    out.close();
}

std::ostream& operator<<(std::ostream& os, Fractal& f) {
    for (Point p : f.r) {
        os << "\n" << p.x << " " << p.y;
    }
    return os;
}
