#pragma once
#include <exception>
#include <string>

class CustomException : public std::runtime_error {
private:
	std::string message;
public:
	CustomException(std::string message) : message(message), std::runtime_error(message) {}
	const char* what() const override;
};
