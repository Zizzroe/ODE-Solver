#include <iostream>
#include <vector>
#include <functional>

using EquationSystem = std::function<std::vector<double>(double, const std::vector<double>&)>;

std::vector<std::pair<double, std::vector<double>>> runge_kutta_system(const EquationSystem& system, double x0, const std::vector<double>& y0, double xn, int n) {
    if (n <= 0) {
        std::cerr << "Error: Number of steps should be greater than 0.\n";
        return {};
    }

    double h = (xn - x0) / n;
    std::vector<std::pair<double, std::vector<double>>> results;
    results.reserve(n + 1); // Reserve space for efficiency

    for (int step = 0; step <= n; ++step) {
        results.emplace_back(x0, y0);

        std::vector<double> k1 = system(x0, y0);
        std::vector<double> k2 = system(x0 + h / 2, y0 + k1 / 2);
        std::vector<double> k3 = system(x0 + h / 2, y0 + k2 / 2);
        std::vector<double> k4 = system(x0 + h, y0 + k3);

        for (size_t i = 0; i < y0.size(); ++i) {
            y0[i] = y0[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }
        x0 = x0 + h;
    }

    return results;
}