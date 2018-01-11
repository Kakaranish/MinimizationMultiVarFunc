#include "CustomException.h"

const char* CustomException::what() const { 
	return message.c_str(); 
}