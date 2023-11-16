#include <iostream>
#include <math.h>
#include <string>
#include <queue>
#include <stack>
#include <cctype>
#include <map>
#include <functional>
#include <regex>

int getPrecedence(char op) {
    switch (op) {
    case '+':
    case '-':
        return 1;
    case '*':
    case '/':
        return 2;
    case '^':
        return 3;
    default:
        return 0;
    }
}

std::queue<std::string> shuntingYard(const std::string& infix) {
    std::queue<std::string> outputQueue;
    std::stack<char> operatorStack;

    for (char token : infix) {
        if (std::isspace(token)) {
            continue;
        }

        if (std::isalnum(token)) {
            outputQueue.push(std::string(1, token));
        }
        else if (token == '(') {
            operatorStack.push(token);
        }
        else if (token == ')') {
            while (!operatorStack.empty() && operatorStack.top() != '(') {
                outputQueue.push(std::string(1, operatorStack.top()));
                operatorStack.pop();
            }
            if (!operatorStack.empty() && operatorStack.top() == '(') {
                operatorStack.pop();
            }
            else {
                std::cerr << "Error: incorrect parenthesis \n";
                return {};
            }
        }
        else if (std::ispunct(token)) {
            while (!operatorStack.empty() && getPrecedence(operatorStack.top()) >= getPrecedence(token)) {
                outputQueue.push(std::string(1, operatorStack.top()));
                operatorStack.pop();
            }
            operatorStack.push(token);
        }
        else {
            std::cerr << "Error: unidentified token" << token << "\n";
            return {};
        }

        std::cout << "Token: " << token << std::endl;
    }
    while (!operatorStack.empty()) {
        if (operatorStack.top() == '(') {
            std::cerr << "Error: mismatched parenthesis \n";
            return {};
        }
        outputQueue.push(std::string(1, operatorStack.top()));
        operatorStack.pop();
    }

    return outputQueue;
}


std::queue<std::string> postfixEval(const std::string& infix) {
    std::queue<std::string> postfixExpression = shuntingYard(infix);
    return postfixExpression;
}

double evaluatePostfix(std::queue<std::string>& postfix, double t) {
    std::stack<double> operandStack;

    while (!postfix.empty()) {
        const std::string& token = postfix.front();

        if (std::isdigit(token[0]) || (token[0] == '-' && token.size() > 1 && std::isdigit(token[1]))) {
            operandStack.push(std::stod(token));
        }
        else if (token == "t") {
            operandStack.push(t);
        }
        else if (token == "+") {
            double b = operandStack.top();
            operandStack.pop();
            double a = operandStack.top();
            operandStack.pop();
            operandStack.push(a + b);
        }
        else if (token == "-") {
            double b = operandStack.top();
            operandStack.pop();
            double a = operandStack.top();
            operandStack.pop();
            operandStack.push(a - b);
        }
        else if (token == "*") {
            double b = operandStack.top();
            operandStack.pop();
            double a = operandStack.top();
            operandStack.pop();
            operandStack.push(a * b);
        }
        else if (token == "/") {
            double b = operandStack.top();
            operandStack.pop();
            double a = operandStack.top();
            operandStack.pop();
            operandStack.push(a / b);
        }

        std::cout << "Token: " << token << " Result: " << operandStack.top() << std::endl;

        postfix.pop();
    }

    if (!operandStack.empty()) {
        return operandStack.top();
    }

    return 0.0;
}


std::string getUserInput() {
    std::string equation;

    std::cout << "**********ODE SOLVER ********** \n";
    std::cout << "Input Your Differential Equation \n";
    std::getline(std::cin, equation);
    return equation;
}

//int gauss_seidel_SLE() {
    //TODO
//}


using EquationFunction = std::function<double(float, float)>;

EquationFunction getEquation() {
    return[](float x, float y) {
        return (y * y - x * x) / (y * y + x * x); //TYPE IN YOUR EQUATION HERE
        };
}

float getStep(float x0, float xn, int n) {
    float stepSize = (xn - x0) / n;
    return stepSize;
}

float trapezoidal_integration(float lowerBound, float upperBound, float interval) {
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

double runge_kutta(const EquationFunction& function, float x0, float y0, float xn, float n) {
    float yn, k1, k2, k3, k4, k;
    int i;

    float h = getStep(x0, xn, n);

    for (i = 0; i < n; i++) {
        k1 = h * function(x0, y0);
        k2 = h * function((x0 + h / 2), (y0 + k1 / 2));
        k3 = h * function((x0 + h / 2), (y0 + k2 / 2));
        k4 = h * function((x0 + h), (y0 + k3));
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        yn = y0 + k;
        x0 = x0 + h;
        y0 = yn;
    }
    return yn;
}

int main() {
    std::string infix = getUserInput();
    std::queue<std::string> postfix = postfixEval(infix);

    if (postfix.empty()) {
        std::cerr << "Error: Invalid Postfix Expression \n";
        return 1;
    }

    std::queue<std::string> copyPostfix = postfix;  // Make a copy for printing
    std::queue<std::string> evalPostfix = postfix;  // Make another copy for evaluation

    // Print the postfix expression
    std::cout << "\nPostfix Expression: ";
    while (!copyPostfix.empty()) {
        std::cout << copyPostfix.front() << " ";
        copyPostfix.pop();
    }
    std::cout << std::endl;

    int EvalType;
    double result;
    double t;
    double tLower, tUpper, stepSize;
    EquationFunction function;
    EquationFunction integrand;
    std::string equation;

    std::cout << "\n**********ODE SOLVER********** \n";
    std::cout << "1. Evaluate ODE At A Static Time \n";
    std::cout << "2. Evaluate ODE Over An Interval Of Time \n";
    std::cout << "3. Evaluate Using 4th Order Runge-Kutta Method \n";
    std::cout << "4. Evaluate Integral of A Function \n";
    std::cout << "5. Evaluate System Of Linear Equations \n";
    std::cin >> EvalType;

    switch (EvalType) {
    case 1:
        std::cout << "Evaluate At What Value Of T: \n";
        std::cin >> t;
        result = evaluatePostfix(evalPostfix, t);
        std::cout << "\nResult: " << result << std::endl;
        break;

    case 2:
        std::cout << "Evaluate Equation From T = : \n";
        std::cin >> tLower;
        std::cout << "Evaluate Equation To T = : \n";
        std::cin >> tUpper;
        std::cout << "Enter Step Size: \n";
        std::cin >> stepSize;

        t = tLower;
        for (double i = tLower; i <= tUpper; i += stepSize) {
            result = evaluatePostfix(evalPostfix, i);
            std::cout << "y(" << i << "): " << result << std::endl;
        }
        break;
    case 3:
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
        runge_kutta(function, x0, y0, xn, n);
        std::cout << "f(" << x0 << ") = " << y0 << "\n";
        break;
    case 4:
        float lowerBound, upperBound, interval;
        std::cout << "Lower Integration Bound: \n";
        std::cin >> lowerBound;
        std::cout << "Upper Integration Bound: \n";
        std::cin >> upperBound;
        std::cout << "Enter Number Of Sub Intervals: \n";
        std::cin >> interval;
        trapezoidal_integration(lowerBound, upperBound, interval);
        std::cout << "Integral Of " << equation << "From " << lowerBound << "To " << upperBound << "Is: \n";
        break;
    case 5:
        break;
    default:
        std::cerr << "Invalid option\n";
        break;
    }

    return 0;
};

//TODO
//FIX CALCULATIONS OVER AN INTERVAL
//CREATE A METHOD FOR SOLVING NSV ODE's
//IMPLEMENT EULER METHOD FOR SOLVING ODE, KEEP SHUNTING YARD ALGORITHM FOR PARSING INPUT INTO EULER METHOD
//CREATE A SYSTEM OF LINEAR EQUATIONS HANDLER FOR NSV ODE
