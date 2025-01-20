#include <iostream>
#include <vector>
#include <ostream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <fstream>
#include <list>
#include <iomanip>

double area(const std::vector<std::pair<double, double>>& poly)
{
    double area = 0;
    size_t s = poly.size();
    for (size_t i = 0; i < s; ++i)
    {
        area += poly[i].first * poly[(i + 1) % s].second - poly[(i + 1) % s].first * poly[i].second;
    }
    area = 0.5 * abs(area);
    return area;
}

double area3(const std::pair<double, double>& p1,const std::pair<double, double>& p2,const std::pair<double, double>& p3)
{
    return abs(p1.first * p2.second + p2.first * p3.second + p3.first * p1.second\
        - p2.first * p1.second - p3.first * p2.second - p1.first * p3.second) / 2;
}

std::pair<double, double> operator-(const std::pair<double, double>& a, const std::pair<double, double>& b)
{
    return { a.first - b.first, a.second - b.second };
}

int sign_of_det(const std::pair<double, double>& p1, const std::pair<double, double>& p2, const std::pair<double, double>& q) {
    std::pair<double, double> a = { p2.first - p1.first, p2.second - p1.second };
    std::pair<double, double> b = { q.first - p1.first, q.second - p1.second };
    double determinant = a.first * b.second - b.first * a.second;
    return determinant < -std::numeric_limits<double>::epsilon();
}

double epsilon = 3.6e-15;

bool sumz_of_area(const std::vector<std::pair<double, double>>& poly,const std::pair<double, double>& point) {
    long double full_area = area(poly);
    long double sum_of_areas = 0;
    size_t s = poly.size();
    for (size_t i = 0; i < s; ++i)
    {
        sum_of_areas += area3(point, poly[i], poly[(i + 1) % s]);
    }
    return abs(full_area - sum_of_areas) < epsilon;
}

bool orientation(const std::vector<std::pair<double, double>>& poly,const std::pair<double, double>& q) {
    size_t s = poly.size();
    for (int i = 0; i < s; i++)
    {
        if (sign_of_det(poly[i], poly[(i + 1) % s], q))
        {
            return false;
        }
    }
    return true;
}

struct Grid {
    double x0, y0;
    double h;

    Grid(double x0_, double y0_, double h_) : x0(x0_), y0(y0_), h(h_) {}

    std::pair<int, int> get_cell(double x, double y) const {
        int i = std::floor((x - x0) / h);
        int j = std::floor((y - y0) / h);
        return { i, j };
    }
    std::vector <std::pair<int, int>> get_cells(double x, double y) const {
        std::vector <std::pair<int, int>> cells; 
        int i = std::floor((x - x0) / h);
        int j = std::floor((y - y0) / h);
        if ((x - x0) / h == std::floor((x - x0) / h) && (y - y0) / h == std::floor((y - y0) / h)) {
            cells.push_back({ i - 1, j });
            cells.push_back({ i - 1, j - 1 });
            cells.push_back({ i, j - 1 });
        }
        else if ((x - x0) / h == std::floor((x - x0) / h)) {
            cells.push_back({ i - 1, j });
        }
        else if ((y - y0) / h == std::floor((y - y0) / h)) {
            cells.push_back({ i, j - 1 });
        }

        cells.push_back({ i, j });
        return cells;
    }

    std::pair<double, double> get_cell_center(int i, int j) const {
        double x = x0 + (i + 0.5) * h;
        double y = y0 + (j + 0.5) * h;
        return { x, y };
    }
};

bool is_point_inside(const std::vector<std::pair<double, double>>& polygon, const std::pair<double, double>& point) {
    int count = 0;
    int n = polygon.size();
    double x0 = point.first;
    double y0 = point.second;

    for (int i = 0; i < n; ++i) {
        std::pair<double, double> p1 = polygon[i];
        std::pair<double, double> p2 = polygon[(i + 1) % n]; 

        if ((p1.first < x0 && p2.first >= x0) || (p2.first < x0 && p1.first >= x0)) {
            double y_intersect;
            if (p1.first == p2.first) {
                continue; 
            }

            y_intersect = p1.second + (x0 - p1.first) * (p2.second - p1.second) / (p2.first - p1.first);
            
            if (y_intersect < y0) {
                count++;
            }
        }
    }
    return count % 2 != 0;
}


std::vector<std::pair<int, int>> get_polygon_cells(const Grid& grid, const std::vector<std::pair<double, double>>& vertices) {
    if (vertices.empty()) return {};

    //double min_x = std::numeric_limits<double>::max();
    //double min_y = std::numeric_limits<double>::max();
    //double max_x = std::numeric_limits<double>::lowest();
    //double max_y = std::numeric_limits<double>::lowest();   
    double min_x = vertices.back().first;
    double min_y = vertices.back().second;
    double max_x = vertices.back().first;
    double max_y = vertices.back().second;

    for (const auto& vertex : vertices) {
        min_x = std::min(min_x, vertex.first);
        min_y = std::min(min_y, vertex.second);
        max_x = std::max(max_x, vertex.first);
        max_y = std::max(max_y, vertex.second);
    }

    std::pair<int, int> min_cell = grid.get_cell(min_x, min_y);
    std::pair<int, int> max_cell = grid.get_cell(max_x, max_y);
    std::vector<std::pair<int, int>> polygon_cells;

    for (int i = min_cell.first; i <= max_cell.first; ++i) {
        for (int j = min_cell.second; j <= max_cell.second; ++j) {
            std::pair<double, double> cell_center = grid.get_cell_center(i, j);
            if (is_point_inside(vertices, cell_center)) {
                polygon_cells.push_back({ i, j });
            }
        }
    }
    std::sort(polygon_cells.begin(), polygon_cells.end());
    polygon_cells.erase(std::unique(polygon_cells.begin(), polygon_cells.end()), polygon_cells.end());
    return polygon_cells;
}

bool grid_method(const Grid& grid, const std::vector<std::pair<int, int>>& polygonCells, const std::pair<double, double>& q) {
    std::vector <std::pair<int, int>> cells = grid.get_cells(q.first, q.second);
    for (const auto& cell : cells) {
        int i = cell.first;
        int j = cell.second;
        for (const auto& polygonCell : polygonCells) {
            if (polygonCell.first == i && polygonCell.second == j) {
                return true;
            }
        }
    }
    return false;
}

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec) {
    for (const auto& v : vec) {
        out << v << " "; 
    }
    return out;
}

int main() {
    std::vector<double> times_sumz;
    std::vector<double> times_orien;
    std::vector<double> times_grid;

    std::list<std::string> files_data = { "Data100.txt", "Data1000.txt", "Data10.txt", "Data10000.txt", "Data100000.txt"};
    for (auto& file : files_data) {
        std::ifstream inputFile(file);
        std::string res_file = "Res" + file;
        std::ofstream outputFile(res_file);
        std::vector<std::pair<double, double>> points;
        std::vector<std::pair<double, double>> polygon;
        int numPoints;
        inputFile >> numPoints;
        for (int i = 0; i < numPoints; ++i) {
            double x, y;
            inputFile >> x >> y;
            points.emplace_back(x, y);
        }
        int numVertices;
        inputFile >> numVertices;
        for (int i = 0; i < numVertices; ++i) {
            double x, y;
            inputFile >> x >> y;
            polygon.emplace_back(x, y);
        }
        inputFile.close();
        int count1 = 0;
        int count2 = 0;
        std::vector<std::pair<double, double>> match_points;
        for (const auto& point : points) {
            if (sumz_of_area(polygon, point)) {
                count1 += 1;
                match_points.push_back(point);
            }
        }
        for (const auto& point : points) {
            if (orientation(polygon, point)) {
                count2 += 1;
            }
        }
        auto begin = std::chrono::high_resolution_clock::now();
        for (const auto& point : points) {
            sumz_of_area(polygon, point);
        }
        auto end = std::chrono::high_resolution_clock::now();
        double t1 = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;
        times_sumz.push_back(t1);
        
        begin = std::chrono::high_resolution_clock::now();
        for (const auto& point : points) {
            orientation(polygon, point);
        }
        end = std::chrono::high_resolution_clock::now();
        double t2 = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;
        times_orien.push_back(t2);
        outputFile << count1 << std::endl;
        for (auto point : match_points) {
            outputFile << std::fixed << std::setprecision(16)
            << point.first << " " << point.second << std::endl;
        }
        outputFile.close();
        std::cout << count1 << "   " << count2 << std::endl;
    }

    std::list<std::string> files_grid_data = { "DataGrid10.txt", "DataGrid100.txt", "DataGrid100.txt", "DataGrid1000.txt", "DataGrid100000.txt" };
        for (auto& file : files_grid_data) {
            std::ifstream inputFile(file);
            std::string res_file = "Res" + file;
            std::ofstream outputFile(res_file); 
            std::vector<std::pair<double, double>> points;
            std::vector<std::pair<double, double>> polygon;
            int numPoints;
            inputFile >> numPoints;
            for (int i = 0; i < numPoints; ++i) {
                double x, y;
                inputFile >> x >> y;
                points.emplace_back(x, y);
            }
            int numVertices;
            inputFile >> numVertices;
            for (int i = 0; i < numVertices; ++i) {
                double x, y;
                inputFile >> x >> y;
                polygon.emplace_back(x, y);
            }
            inputFile.close();
            double x0 = 0.0;
            double y0 = 0.0;
            double h = 1.0;
            int count = 0;
            std::vector<std::pair<double, double>> match_points;
            Grid grid(x0, y0, h);
            std::vector<std::pair<int, int>> polygonCells = get_polygon_cells(grid, polygon);
            for (const auto& point : points) {
                if (grid_method(grid, polygonCells, point)) {
                    count += 1;
                    match_points.push_back(point);
                }
            }
            auto begin = std::chrono::high_resolution_clock::now();
            for (const auto& point : points) {
                grid_method(grid, polygonCells, point);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double t = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;
            times_grid.push_back(t);
            outputFile << count << std::endl;
            for (const auto& point : match_points) {
                outputFile << std::fixed << std::setprecision(16)
                    << point.first << " " << point.second << std::endl;
            }
            outputFile.close();
            std::cout << count << std::endl;
        }
        std::cout << "times_sumz:  " << std::fixed << std::setprecision(7) << times_sumz << std::endl;
        std::cout << "times_orien: " << std::fixed << std::setprecision(7) << times_orien << std::endl;
        std::cout << "times_grid: " << std::fixed << std::setprecision(7) << times_grid << std::endl;
        return 0;
}




