#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

using EquationSystem = std::function<std::vector<double>(double, const std::vector<double>&)>;

std::vector<std::pair<double, std::vector<double>>> nsv_euler(const EquationSystem& system, double t0, const std::vector<double>& y0, double tn, int n) {
    if (n <= 0) {
        std::cerr << "Error: Number of steps should be greater than 0.\n";
        return {};
    }

    double h = (tn - t0) / n;
    std::vector<std::pair<double, std::vector<double>>> results;
    results.reserve(n + 1); // Reserve space for efficiency

    std::vector<double> y = y0;

    for (int step = 0; step <= n; ++step) {
        results.emplace_back(t0 + step * h, y);

        std::vector<double> dy = system(t0 + step * h, y);

        for (size_t i = 0; i < y.size(); ++i) {
            y[i] = y[i] + h * dy[i];
        }
    }

    return results;
}

int main() {
    // Define your NSV ODE system
    EquationSystem nsvSystem = [](double t, const std::vector<double>& y) -> std::vector<double> {
        // Example system: dy1/dt = -y1, dy2/dt = y1 - y2
        return { -y[0], y[0] - y[1] };
        };

    // Initial conditions
    double t0 = 0.0;
    std::vector<double> y0 = { 1.0, 0.0 };  // Initial values of y1 and y2

    // Time span and number of steps
    double tn = 5.0;
    int n = 100;

    // Solve the NSV ODE system using Euler's method
    std::vector<std::pair<double, std::vector<double>>> result = nsv_euler(nsvSystem, t0, y0, tn, n);

    // Display the results
    std::cout << "Results for NSV ODE system using Euler's method:\n";
    for (const auto& pair : result) {
        std::cout << "t: " << pair.first << ", y: " << pair.second[0] << ", " << pair.second[1] << '\n';
    }

    return 0;
}
