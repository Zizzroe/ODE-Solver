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
#include <tuple>

#define SIZE 10 

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


using EquationFunction = std::function<double(float, float)>;

EquationFunction getEquation() {
    return[](float x, float y) {
        return (y * y - x * x ) / ( y * y + x * x ); //TYPE IN YOUR EQUATION HERE
        };
}

float getStep(float x0, float xn, int n) {
    float stepSize = (xn - x0) / n;
    return stepSize;
}

float trapezoidal_integration(float lowerBound, float upperBound, float interval) { //RECOMMENDED VALUE FOR INTERVAL IS MINIMUM OF 500
#define integrand(x) pow(x, x) //TYPE IN YOUR INTEGRAND HERE
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

std::vector<std::tuple<double, double, double>> runge_kutta(const EquationFunction& function, float x0, float y0, float xn, float n) {
    double yn;
    if (n <= 0) {
        std::cerr << "Error: Number of steps should be greater than 0.\n";
        return {};
    }

    double h = getStep(x0, xn, n);
    std::vector<std::tuple<double, double, double>> results;
    results.reserve(n + 1); // Reserve space for efficiency

    for (int step = 0; step <= n; ++step) {

        double k1 = h * function(x0, y0);
        double k2 = h * function(x0 + h / 2, y0 + k1 / 2);
        double k3 = h * function(x0 + h / 2, y0 + k2 / 2);
        double k4 = h * function(x0 + h, y0 + k3);
        double k = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        yn = y0 + k;
        x0 = x0 + h;
        y0 = yn;
        
        results.emplace_back(x0, y0, yn);
    }

    return results;
}


int main() {

    //CASE DETERMINANT
    int EvalType;
    //RK4
    EquationFunction function;
    std::vector<std::tuple<double, double, double>> resultRK4;
    //TRAPEZOIDAL_INTEGRATION
    float resultIntegral;
    double t;
    double tLower, tUpper, stepSize;
    std::string equation;
    //SLE
    float coefficients[SIZE][SIZE + 1], solutions[SIZE];
    int n;

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
            std::cout << "x: " << std::get<0>(point) << ", y: " << std::get<1>(point) << " yn: " << std::get<2>(point) <<"\n";
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
        std::cout << "Integral " << equation << "From " << lowerBound << " To " << upperBound << " Is: " << std::setprecision(15) << resultIntegral << "\n";
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
    default:
        std::cerr << "Invalid option\n";
        break;
    }

    return 0;
};

//TODO
//CREATE A METHOD FOR SOLVING NSV ODE's
//CREATE A SYSTEM OF LINEAR EQUATIONS HANDLER FOR NSV ODE
