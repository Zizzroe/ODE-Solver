#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <string>
#include <queue>
#include <stack>
#include <cctype>
#include <map>
#include <functional>
#include <regex>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <fstream>

#define SIZE 10 

using EquationFunction = std::function<double(float, float)>;
using EquationSystem = std::function<std::vector<double>(double, const std::vector<double>&)>;


void gaussJordanElimination(int n, float coefficients[SIZE][SIZE + 1], float solutions[SIZE]) {
    float ratio;

    // Applying Gauss Jordan Elimination
    for (int i = 1; i <= n; i++) {
        if (coefficients[i][i] == 0.0) {
            std::cout << "Mathematical Error!";
            std::exit(0);
        }
        for (int j = 1; j <= n; j++) {
            if (i != j) {
                ratio = coefficients[j][i] / coefficients[i][i];
                for (int k = 1; k <= n + 1; k++) {
                    coefficients[j][k] = coefficients[j][k] - ratio * coefficients[i][k];
                }
            }
        }
    }

    // Obtaining Solution
    for (int i = 1; i <= n; i++) {
        solutions[i] = coefficients[i][n + 1] / coefficients[i][i];
    }
}

EquationFunction getEquation() {
    return[](float x, float y) {
        return (y * y - x * x) / (y * y + x * x); //TYPE IN YOUR EQUATION HERE
        };
}

float getStep(float x0, float xn, int n) {
    float stepSize = (xn - x0) / n;
    return stepSize;
}

float trapezoidal_integration(float lowerBound, float upperBound, float interval) { //RECOMMENDED VALUE FOR INTERVAL IS MINIMUM OF 500
#define integrand(x) 1 / (1 + pow(x, 2)) //TYPE IN YOUR INTEGRAND HERE
    float integration = 0.0;
    float stepSize, k;

    stepSize = (upperBound - lowerBound) / interval;
    integration = integrand(lowerBound) + integrand(upperBound);

    for (int i = 1; i < interval; i++) {
        k = lowerBound + i * stepSize;
        integration += 2 * integrand(k);
    }

    integration = integration * (stepSize / 2);

    return integration;
}

std::vector<std::pair<double, double>> runge_kutta(const EquationFunction& function, float x0, float y0, float xn, float n) {
    if (n <= 0) {
        std::cerr << "Error: Number of steps should be greater than 0.\n";
        return {};
    }

    double h = getStep(x0, xn, n);
    std::vector<std::pair<double, double>> results;
    results.reserve(n + 1); // Reserve space for efficiency

    for (int step = 0; step <= n; ++step) {
        results.emplace_back(x0, y0);

        double k1 = h * function(x0, y0);
        double k2 = h * function(x0 + h / 2, y0 + k1 / 2);
        double k3 = h * function(x0 + h / 2, y0 + k2 / 2);
        double k4 = h * function(x0 + h, y0 + k3);

        y0 = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x0 = x0 + h;
    }

    return results;
}

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

void cashKarpMethod(double t, std::vector<double>& y, double& h, const EquationSystem& system, double tolerance) {
    const double a2 = 1.0 / 5.0;
    const double a3 = 3.0 / 10.0;
    const double a4 = 3.0 / 5.0;
    const double a5 = 1.0;
    const double a6 = 7.0 / 8.0;

    const double b21 = 1.0 / 5.0;

    const double b31 = 3.0 / 40.0;
    const double b32 = 9.0 / 40.0;

    const double b41 = 3.0 / 10.0;
    const double b42 = -9.0 / 10.0;
    const double b43 = 6.0 / 5.0;

    const double b51 = -11.0 / 54.0;
    const double b52 = 5.0 / 2.0;
    const double b53 = -70.0 / 27.0;
    const double b54 = 35.0 / 27.0;

    const double b61 = 1631.0 / 55296.0;
    const double b62 = 175.0 / 512.0;
    const double b63 = 575.0 / 13824.0;
    const double b64 = 44275.0 / 110592.0;
    const double b65 = 253.0 / 4096.0;

    const double c1 = 37.0 / 378.0;
    const double c3 = 250.0 / 621.0;
    const double c4 = 125.0 / 594.0;
    const double c6 = 512.0 / 1771.0;

    const double dc1 = c1 - 2825.0 / 27648.0;
    const double dc3 = c3 - 18575.0 / 48384.0;
    const double dc4 = c4 - 13525.0 / 55296.0;
    const double dc5 = -277.0 / 14336.0;
    const double dc6 = c6 - 0.25;

    std::vector<double> k1 = system(t, y);

    std::vector<double> y2(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y2[i] = y[i] + b21 * h * k1[i];
    }

    std::vector<double> k2 = system(t + a2 * h, y2);

    std::vector<double> y3(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y3[i] = y[i] + h * (b31 * k1[i] + b32 * k2[i]);
    }

    std::vector<double> k3 = system(t + a3 * h, y3);

    std::vector<double> y4(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y4[i] = y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
    }

    std::vector<double> k4 = system(t + a4 * h, y4);

    std::vector<double> y5(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y5[i] = y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
    }

    std::vector<double> k5 = system(t + a5 * h, y5);

    std::vector<double> y6(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y6[i] = y[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
    }

    std::vector<double> k6 = system(t + a6 * h, y6);

    // Calculate the error
    double error = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        double delta = c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i];
        error += delta * delta;
    }
    error = std::sqrt(error / y.size());

    // Adjust the step size
    double factor = 0.9 * std::pow(tolerance / error, 0.2);
    h *= factor;

    // Update the solution
    if (error < tolerance) {
        for (size_t i = 0; i < y.size(); ++i) {
            y[i] = y[i] + dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + dc6 * k6[i];
        }
    }
}

std::vector<std::pair<double, std::vector<double>>> nsv_cash_karp(const EquationSystem& system, double t0, const std::vector<double>& y0, double tn, double tolerance) {
    const int maxSteps = 10000; // Maximum number of steps to avoid infinite loop

    double h = 100000; // Initial step size
    double minH = 1e-6; // Minimum step size
    double t = t0;
    std::vector<double> y = y0;

    std::vector<std::pair<double, std::vector<double>>> results;
    results.emplace_back(t, y);

    for (int step = 0; t < tn && step < maxSteps; ++step) {
        cashKarpMethod(t, y, h, system, tolerance);

        // Ensure that h is not too small
        h = std::max(h, minH);

        t += h;
        results.emplace_back(t, y);
    }

    return results;
}


std::vector<double> getLorenzSystem(double t, const std::vector<double>& state) {
    const double sigma = 10.0;
    const double rho = 28.0;
    const double beta = 8.0 / 3.0;

    double x = state[0];
    double y = state[1];
    double z = state[2];

    double dxdt = sigma * (y - x);
    double dydt = x * (rho - z) - y;
    double dzdt = x * y - beta * z;

    return { dxdt, dydt, dzdt };
}


int main() {

    //CASE DETERMINANT
    int EvalType;
    //SHARED VARIABLES
    int n; //SLE & STRANGE ATTRACTORS
    //RK4
    EquationFunction function;
    std::vector<std::pair<double, double>> resultRK4;
    //RKCK5
    double tolerance;
    //TRAPEZOIDAL_INTEGRATION
    float resultIntegral;
    double t;
    double tLower, tUpper, stepSize;
    std::string equation;
    //SLE
    float coefficients[SIZE][SIZE + 1], solutions[SIZE];
    //STRANGE ATTRACTORS;
    double t0, tn;
    std::vector<double> initialState;
    std::vector<std::pair<double, std::vector<double>>> results;

    std::cout << "\n**********ODE SOLVER********** \n";
    std::cout << "1. Evaluate Using 4th Order Runge-Kutta Method \n";
    std::cout << "2. Evaluate Integral of A Function \n";
    std::cout << "3. Evaluate System Of Linear Equations \n";
    std::cin >> EvalType;

    switch (EvalType) {
    
    case 1:
        float x0, y0, xn, n;
        std::cout << "Enter Initial Conditions \n";
        std::cout << "y0: \n";
        std::cout << "x0: \n";
        std::cin >> y0;
        std::cin >> x0;
        std::cout << "Enter Calculation Point: \n";
        std::cout << "xn: \n";
        std::cin >> xn;
        std::cout << "Enter Number Of Steps: \n";
        std::cin >> n;
        function = getEquation();
        resultRK4 = runge_kutta(function, x0, y0, xn, n);
        std::cout << "Results for Runge-Kutta method:\n";
        for (const auto& point : resultRK4) {
            std::cout << "x: " << point.first << ", y: " << point.second << "\n";
        }
        break;
    case 2:
        float lowerBound, upperBound, interval;
        std::cout << "Lower Integration Bound: \n";
        std::cin >> lowerBound;
        std::cout << "Upper Integration Bound: \n";
        std::cin >> upperBound;
        std::cout << "Enter Number Of Sub Intervals: \n";
        std::cin >> interval;
        resultIntegral = trapezoidal_integration(lowerBound, upperBound, interval);
        std::cout << "Integral " << equation << "From " << lowerBound << " To " << upperBound << " Is: " << resultIntegral << "\n";
        break;
    case 3:
        std::cout << "Enter The # Of Unknown Variables In Your System Of Linear Equations: \n";
        std::cin >> n;
        std::cout << "Enter Coefficients of Augmented Matrix: " << std::endl;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n + 1; j++) {
                std::cout << "a[" << i << "]" << j << "= ";
                std::cin >> coefficients[i][j];
            }
        }
        gaussJordanElimination(n, coefficients, solutions);
        std::cout << std::endl << "Solution: " << std::endl;
        for (int i = 1; i <= n; i++) {
            std::cout << "x[" << i << "] = " << solutions[i] << std::endl;
        }
        break;
    case 4:
        std::cout << "Case 4 is just that you already defined some other stuff in the code because you are code literate \n";
        std::cout << "And In This Case Its Just Some Strange Attractors \n";

        t0 = 0.0; // initial time;
        tn = 25.0; // run-to time;
        n = 5000;
        initialState = { 0.0, 0.1, 0.0 }; // initial conditions
        results = nsv_euler(getLorenzSystem, t0, initialState, tn, n);

        // Print the results
        for (const auto& point : results) {
            std::cout << "t: " << point.first << ", x: " << point.second[0] << ", y: " << point.second[1] << ", z: " << point.second[2] << "\n";
        }

        // Ask the user if they want to save the output
        char saveOutput;
        std::cout << "Do you want to save the output to a file? (y/n): ";
        std::cin >> saveOutput;

        if (saveOutput == 'y' || saveOutput == 'Y') {
            // Open a file for writing
            std::ofstream outputFile("output.txt");

            // Check if the file is successfully opened
            if (outputFile.is_open()) {
                // Write the results to the file
                for (const auto& point : results) {
                    outputFile << "x: " << point.second[0] << " y : " << point.second[1] << " z : " << point.second[2] << "\n";
                }

                // Close the file
                outputFile.close();
                std::cout << "Results written to output.txt\n";
            }
            else {
                std::cerr << "Error opening output.txt for writing\n";
            }
        }
        break;

    default:
        std::cerr << "Invalid option\n";
        break;
    }

    return 0;
};

//TODO
//CREATE A METHOD FOR SOLVING NSV ODE's
//CREATE A SYSTEM OF LINEAR EQUATIONS HANDLER FOR NSV ODE
