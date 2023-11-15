#include <iostream>
#include <string>
#include <queue>
#include <stack>
#include <cctype>
#include <map>
#include <functional>

int getPrecedence(char op) { //checks order of operations, assigns precedence value for operator
	switch (op) {			 //variable 'operator' cannot be used as operator is a build in function
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
			continue; //ignores all whitespace in input
		}

		if (std::isalnum(token)) {
			//if token is a number, push it to the outputeQueue
			outputQueue.push(std::string(1, token));

		}
		else if (token == '(') {
			//if token is a openinng parenthesis, push it to the operatorStack
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
			// If the token is a number, push it onto the operand stack
			operandStack.push(std::stod(token));
		}
		else if (token == "t") {
			// If the token is 't', push the current time onto the operand stack
			operandStack.push(t);
		}
		else if (token == "+") {
			// Addition
			double b = operandStack.top();
			operandStack.pop();
			double a = operandStack.top();
			operandStack.pop();
			operandStack.push(a + b);
		}
		else if (token == "-") {
			// Subtraction
			double b = operandStack.top();
			operandStack.pop();
			double a = operandStack.top();
			operandStack.pop();
			operandStack.push(a - b);
		}
		else if (token == "*") {
			// Multiplication
			double b = operandStack.top();
			operandStack.pop();
			double a = operandStack.top();
			operandStack.pop();
			operandStack.push(a * b);
		}
		else if (token == "/") {
			// Division
			double b = operandStack.top();
			operandStack.pop();
			double a = operandStack.top();
			operandStack.pop();
			operandStack.push(a / b);
		}

		postfix.pop();
	}

	if (!operandStack.empty()) {
		return operandStack.top();
	}

	// Return 0.0 if the stack is empty (indicating an issue in evaluation)
	return 0.0;
}


std::string getUserInput() {
	std::string equation;
	std::cout << "Enter the differential equation, please use 't' as the independent variable" << std::endl;
	std::getline(std::cin, equation);
	return equation;
}


int main() {
	std::string infix = getUserInput();
	std::queue<std::string> postfix = postfixEval(infix);

	if (postfix.empty()) {
		std::cerr << "Error: Invalid postfix expression." << std::endl;
		return 1; // Return an error code
	}

	std::queue<std::string> copyPostfix = postfix;  // Make a copy

	while (!copyPostfix.empty()) {
		std::cout << copyPostfix.front() << " ";
		copyPostfix.pop();
	}

	// Now you can proceed with evaluating the postfix expression
	double result = evaluatePostfix(postfix, 2);
	std::cout << "\nResult: " << result << std::endl;

	return 0;
}
