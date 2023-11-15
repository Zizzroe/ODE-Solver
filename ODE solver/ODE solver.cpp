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

int gauss_seidel_SLE() {
//TODO
}


float getStep(float x0, float xn, int n) {
    float stepSize = (xn - x0) / n;
    return stepSize;
}

int trapezoidal_integration(std::string equation, float lowerBound, float upperBound, float interval) {
    #define f(x) equation
    float integration = 0.0;
    float stepSize, k;
    int i;

    stepSize = (upperBound - lowerBound) / interval;
    integration = f(lowerBound) + f(upperBound);

    for (i = 0; i <= interval - 1; i++) {
        k = lowerBound + i * stepSize;
        integration = integration + 2 * (f(k));
    }
    integration = integration * (stepSize / 2);

    return integration;
}

int runge_kutta(std::string equation, float x0, float y0, float xn*, float n*) {
    #define f(x,y) equation;
    float yn, k1, k2, k3, k4, k;
    int i;

    float h = getStep(x0, xn, n);
    
    for (i = 0, i < n; i++) {
        k1 = h * (f(x0, y0));
        k2 = h * (f((x0 + h / 2), (y0 + k1 / 2);
        k3 = h * (f((x0 + h / 2), (y0 + k2 / 2);
        k4 = h * (f((x0 + h), (y0 + k3);
        k = (k1 + 2 * k2 + 2 * k3 * k4) / 6;
        yn = y0 + k;
        x0 = x0 + h;
        y0 = yn;
    }
    return yn, xn;
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

    std::cout << "\n**********ODE SOLVER********** \n"
        << "1. Evaluate ODE At A Static Time \n"
        << "2. Evaluate ODE Over An Interval Of Time \n";
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

    default:
        std::cerr << "Invalid option\n";
        break;
    }

    return 0;
}

//TODO
//FIX CALCULATIONS OVER AN INTERVAL
//CREATE A METHOD FOR SOLVING NSV ODE's
//IMPLEMENT EULER METHOD FOR SOLVING ODE, KEEP SHUNTING YARD ALGORITHM FOR PARSING INPUT INTO EULER METHOD
//CREATE A SYSTEM OF LINEAR EQUATIONS HANDLER FOR NSV ODE's
