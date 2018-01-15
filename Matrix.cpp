#include "Matrix.h"

//Private functions
void Matrix::deallocate2DArray(double** _arr, unsigned rowsNum) { // Deallocating memory 
	for (unsigned i = 0; i < rowsNum; i++)
		delete[] _arr[i];
	delete[] _arr;
}

//Constructors and destructor
Matrix::Matrix(unsigned rowNum, unsigned colNum) : size(MatrixSize(rowNum,colNum)){
	arr = new double*[size.getRowsNum()];
	for (int i = 0; i < size.getRowsNum(); i++)
		arr[i] = new double[size.getColNum()];

	init();
}
Matrix::Matrix(MatrixSize size): size(size) {	
	arr = new double*[size.getRowsNum()];
	for (int i = 0; i < size.getRowsNum(); i++)
		arr[i] = new double[size.getColNum()];
	init();
}
Matrix::Matrix(const Matrix& matrix) {
	size = matrix.size;
	arr = new double*[size.getRowsNum()];
	for (int i = 0; i < size.getRowsNum(); i++) {
		arr[i] = new double[size.getColNum()];
		for (int j = 0; j < size.getColNum(); j++) {
			arr[i][j] = matrix.arr[i][j];
		}
	}
}
Matrix::Matrix(Matrix&& movedMatrix) noexcept{
	arr = std::move(movedMatrix.arr);
	size = std::move(movedMatrix.size);
	movedMatrix.arr = nullptr;
}
Matrix::~Matrix() {
	if (arr == nullptr)
	{
		//std::cout << "Nothing to deallocate" << std::endl;
		return; //nothing to deallocate
	}
	
	deallocate2DArray(this->arr, size.getRowsNum());
}


//Operations on matrices
void Matrix::multiplyByScalar(double scalar) {
	Matrix temp(*this);
	for (unsigned i = 0; i < size.getRowsNum(); i++)
		for (unsigned j = 0; j < size.getColNum(); j++)
			temp(i, j) = temp(i, j) * scalar;
	*this = std::move(temp);
}
Matrix Matrix::getMultipliedByScalar(double scalar){
	Matrix temp(*this);
	for(unsigned i=0;i<size.getRowsNum();i++)
		for(unsigned j=0;j<size.getColNum();j++)
			temp(i,j) = temp(i,j) * scalar;
	return std::move(temp);
}
void Matrix::multiplyByMatrix(Matrix const& B) {
	if (size.getColNum() != B.size.getRowsNum())
		throw CustomException("Unable to multiply matrices! Wrong dimensions!");

	double** tmpArr = arr;
	unsigned prevColNum = size.getColNum();
	size = MatrixSize(size.getRowsNum(), B.size.getColNum());

	arr = nullptr;
	arr = new double*[size.getRowsNum()];
	for (int i = 0; i < size.getRowsNum(); i++) {
		arr[i] = new double[size.getColNum()];
		for (int j = 0; j < size.getColNum(); j++) {
			double sum = 0;
			for (unsigned k = 0; k < prevColNum; k++)
				sum += (tmpArr[i][k] * B.arr[k][j]);
			arr[i][j] = sum;
		}
	}
	deallocate2DArray(tmpArr, size.getRowsNum());
}
Matrix Matrix::getMultipliedByMatrix(Matrix const &B) {
	if (size.getColNum() != B.size.getRowsNum())
		throw CustomException("Unable to get multiplied matrices! Wrong dimensions!");

	Matrix rMatrix(size.getRowsNum(), B.size.getColNum());
	for (unsigned i = 0; i < size.getRowsNum(); i++) //iterate A matrix through columns
		for (unsigned j = 0; j < B.size.getColNum(); j++) { //iterate B matrix through rows
			double sum = 0;
			for (unsigned k = 0; k < size.getColNum(); k++)
				sum += (arr[i][k] * B.arr[k][j]);
			rMatrix(i, j) = sum;
		}
	return rMatrix;
}
void Matrix::transpose() {
	if (size.getColNum() != size.getRowsNum()) {

		double** tmpPtr = arr;
		arr = nullptr;


		this->size = MatrixSize(size.getColNum(), size.getRowsNum());
		arr = new double*[size.getRowsNum()];
		for (int i = 0; i < size.getRowsNum(); i++) {
			arr[i] = new double[size.getColNum()];
			for (int j = 0; j < size.getColNum(); j++) 
				arr[i][j] = tmpPtr[j][i];
		}
		//We want to deallocate previous arr after copying values
		deallocate2DArray(tmpPtr, size.getColNum());
	}
	else {
		for (int i = 0; i < size.getRowsNum(); i++)
			for (int j = 0; j < size.getColNum(); j++)
				if(i < j) //2x swap gives the same result, so we wanna avoid this
					std::swap(arr[i][j], arr[j][i]);
	}
}
Matrix Matrix::getTransposed() {
	Matrix tMatrix(size.getColNum(), size.getRowsNum());

	for (unsigned i = 0; i < size.getColNum(); i++) 
		for (unsigned j = 0; j < size.getRowsNum(); j++)
			tMatrix(i,j) = arr[j][i];

	return tMatrix;
}
void Matrix::interchangeRows(unsigned r1, unsigned r2) {
	if (r1 >= size.getRowsNum() || r2 >= size.getRowsNum()) 
		throw CustomException("Unable to interchangeRows!");
		if (r1 == r2)
			return;
		for (unsigned i = 0; i < size.getColNum(); i++)
			std::swap(arr[r1][i], arr[r2][i]);
}
void Matrix::addRow(unsigned baseRow, unsigned addedRow, double multiplier) {
	//because baseRow and addedRow are unsigned, that means if we put values that's less than zero
	//in function it will be equals to max unsigned value 
	if (baseRow >= size.getRowsNum() || addedRow >= size.getRowsNum()) 
		throw CustomException("Out of range! Unable to add row!");
	
	for (unsigned i = 0; i < size.getColNum(); i++)
		arr[baseRow][i] += arr[addedRow][i] * multiplier;
}


//Vector functions
bool Matrix::isVector() const {
	return ((size.getColNum() > 0 && size.getRowsNum() > 0)
		&& (size.getRowsNum() == 1 || size.getColNum() == 1)) ? true : false;
}
bool Matrix::isVerticalVector()  const{
	if (!isVector())
		throw CustomException("Can't check if vector is vertical! Matrix isn't vector!");
	return (size.getRowsNum() == 1) ? false : true;
}
double Matrix::getVectorLength(){
	if (!isVector())
		throw CustomException("Unable to get vector lenght! Matrix isn't matrix");

	double sum=0;
	if(size.getColNum() > 1){
		for(unsigned i=0;i<size.getColNum();i++)
			sum += pow(arr[0][i],2);
	}else{
		for(unsigned i=0;i<size.getRowsNum();i++)
			sum += pow(arr[i][0],2);
	}
	return sqrt(sum);
}
void Matrix::normalizeVector(){
	if (!this->isVector())
		throw CustomException("Unnable to normalize vector! Matrix isn't vector");
	this->multiplyByScalar(1.f/getVectorLength());
}
Matrix Matrix::getNormalizedVector() {
	if (!this->isVector())
		throw CustomException("Unnable to get normalized vector! Matrix isn't vector");
	Matrix tmpMatrix(*this);
	tmpMatrix.multiplyByScalar(1.f / getVectorLength());
	return tmpMatrix;
}
double Matrix::getScalarProduct(const Matrix& matrix){
	if (!matrix.isVector() || !this->isVector() || matrix.size != this->size)
		throw CustomException("Unable to get scalar product!");
	
	double sum = 0;
	if(size.getColNum() > 1)
		for(int i=0;i<size.getColNum();i++)
			sum += arr[0][i] * matrix.getArrPtr()[0][i];
	else
		for(int i=0;i<size.getRowsNum();i++)
			sum += arr[i][0] * matrix.getArrPtr()[i][0];

	return sum;
}
Matrix Matrix::findPerpendicularVector(){
	if(!this->isVector())
		throw CustomException("Unable to find perpendicular vector! Matrix isn't vector!");
	Matrix temp(this->size);
	double sum = 0;
	if(size.getColNum() > 1)
		for(unsigned i=0;i<size.getColNum()-1;i++){
			temp(0,i) = 1;
			sum += temp(0,i) * arr[0][i];
		}
	else
		for(unsigned i=0;i<size.getRowsNum()-1;i++){
			temp(i,0) = 1;
			sum += temp(i,0) * arr[i][0];
		}

	double last;
	if(size.getColNum() > 1){
		last = - (sum / arr[0][size.getColNum()-1]);
		temp(0, size.getColNum() - 1) = last;
	}	else{
		last = -(sum / arr[size.getRowsNum()-1][0]);
		temp(size.getRowsNum() - 1, 0 ) = last;
	}

	return std::move(temp);
}


//Overloading operators
double& Matrix::operator()(unsigned row, unsigned col) {
	if (row >= size.getRowsNum() || col >= size.getColNum())
		throw CustomException("Out of range! Unable to get access to specific element!");
	return arr[row][col];
}
bool Matrix::operator==(const Matrix& matrix) {
	if (size != matrix.size)
		return false;
	for (int i = 0; i < size.getRowsNum(); i++)
		for (int j = 0; j < size.getColNum(); j++)
			if (arr[i][j] != matrix.arr[i][j])
				return false;
	return true;
}
bool Matrix::operator!=(const Matrix& matrix) {
	if (size != matrix.size)
		return true;
	for (int i = 0; i < size.getRowsNum(); i++)
		for (int j = 0; j < size.getColNum(); j++)
			if (arr[i][j] != matrix.arr[i][j])
				return true;
	return false;
}
Matrix& Matrix::operator=(Matrix & B) {
	if (arr != nullptr) 
		deallocate2DArray(arr, size.getRowsNum());

	size = B.size;
	arr = new double*[size.getRowsNum()];
	for (int i = 0; i < size.getRowsNum(); i++) {
		arr[i] = new double[size.getColNum()];
		for (int j = 0; j < size.getColNum(); j++) 
			arr[i][j] = B.arr[i][j];		
	}
	return *this;
}
Matrix& Matrix::operator=(Matrix && B){
	deallocate2DArray(this->arr, size.getRowsNum());
	this->size = std::move(B.size);

	arr = std::move(B.arr);

	B.arr = nullptr;
	return *this;
}
Matrix& Matrix::operator+=(const Matrix& matrix) {
	if (this->size != matrix.size)
		throw CustomException("Can't use operator+= because of different dimensions!");
	
	for (unsigned i = 0; i<size.getRowsNum(); i++)
		for (unsigned j = 0; j<size.getColNum(); j++)
			arr[i][j] += matrix.arr[i][j];
	return *this;
}
Matrix& Matrix::operator-=(const Matrix& matrix){
	if(this->size != matrix.size)
		throw CustomException("Can't use operator-= because of different dimensions!");

	for(unsigned i=0;i<size.getRowsNum();i++)
		for(unsigned j=0;j<size.getColNum();j++)
			arr[i][j] -= matrix.arr[i][j];
	return *this;
}
Matrix Matrix::operator+(const Matrix& matrix) {
	if (this->size != matrix.size)
		throw CustomException("Can't use operator+ because of different dimensions!");

	Matrix tMatrix(size);
	for (unsigned i = 0; i<size.getRowsNum(); i++)
		for (unsigned j = 0; j<size.getColNum(); j++)
			tMatrix.arr[i][j] = (arr[i][j] + matrix.arr[i][j]);
	return tMatrix;
}
Matrix Matrix::operator-(const Matrix& matrix) {
	if (this->size != matrix.size)
		throw CustomException("Can't use operator+ because of different dimensions!");

	Matrix tMatrix(getSize());
	for (unsigned i = 0; i<size.getRowsNum(); i++)
		for (unsigned j = 0; j<size.getColNum(); j++)
			tMatrix.arr[i][j] = (arr[i][j] - matrix.arr[i][j]);
	return tMatrix;
}
std::ostream& operator<<(std::ostream& output, const Matrix& _m) {
	for (unsigned i = 0; i < _m.size.getRowsNum(); i++) {
		for (unsigned j = 0; j < _m.size.getColNum(); j++)
			output << std::setw(8) << _m.arr[i][j] << " ";
		output << std::endl;
	}
	return output;
}


//Misc functions
bool Matrix::isInitialized(){
	if (size.getRowsNum() > 0 && size.getColNum() > 0)
		return true;
	return false;
}
void Matrix::init() { // Filling the array with 0
	for (unsigned i = 0; i < size.getRowsNum(); i++)
		for (unsigned j = 0; j < size.getColNum(); j++)
			arr[i][j] = 0;
}
void Matrix::initRandomValues(int intervalBeg, int intervalEnd) {	
	if (intervalBeg > intervalEnd)
		std::swap(intervalBeg, intervalEnd);

	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_int_distribution<int> distribution(intervalBeg, intervalEnd);
	
	for (int i = 0; i < size.getRowsNum(); i++)
		for (int j = 0; j < size.getColNum(); j++)
			arr[i][j] = distribution(generator);
}
