#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

ifstream fileOut;

static double *matrix;
static int numOfVariables; 
static int numOfEquations;
static int numOfIterations = 0;

void getMatrix(int numOfVar, int numOfEq, double *matrix/*out*/) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			fileOut >> matrix[i * numOfVar + j];
		}
	}
}

void getMartixSize(int *numOfVar/*out*/, int *numOfEq/*out*/) {
	fileOut >> *numOfEq >> *numOfVar;
	return;
}

void initMatrixByNum(double *matrix, int numOfVar, int numOfEq, int num) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			matrix[i * numOfVar + j] = num;
		}
	}
	return;
}

double normMatrix(double *A, int numOfVar, int numOfEq) {
	double sum = 0;
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			double tmpValue = 0;
			tmpValue = A[i * numOfVar + j];
			tmpValue = powl(tmpValue, 2);
			sum += tmpValue;
		}
	}
	return sqrt(sum);
}

double BisectionMethod(double(*Function)(double root), double error, double a, double b) {

	if (!(Function(a) * Function(b) < 0)) {
		return INFINITY;
	}

	double Fa = Function(a);
	double c = (a + b) / 2;
	double Fc = Function(c);

	while (abs(Fc) > error) {

		if (Fa * Fc > 0) {
			a = c;
		}
		else if (Fa * Fc < 0) {
			b = c;
		}

		c = (a + b) / 2;
		Fc = Function(c);
		//numOfIterations++;
	}
	//cout<<numOfIterations;

	return c;
}

double HessenbergCharacterPolyVal(double eigenValue) {
	double *masOfDeterminants = new double [numOfEquations + 1];

    masOfDeterminants[0] = 1;
	masOfDeterminants[1] = matrix[0] - eigenValue;
	
	for (int k = 0; k < numOfEquations; k++) {
		masOfDeterminants[k + 1] = (matrix[k * numOfVariables + k] - eigenValue) * masOfDeterminants[k];

		for (int j = 1; j <= k; j++) {
			double tmpVal = 1;

			for (int i = 1; i <= j; i++) {
				tmpVal *= matrix[(k + 1 - i) * numOfVariables + (k - i)];
			}

			masOfDeterminants[k + 1] += ((j % 2) == 0 ? 1 : -1) * matrix[(k - j) * numOfVariables + k] * tmpVal * masOfDeterminants[k - j];
		}
	}
	
	double result = masOfDeterminants[numOfEquations];
	delete [] masOfDeterminants;

    return result;
}

double *GetRoots(double *matrix, int numOfVar, int numOfEq) {
	double error = 0.001;//affect accuracy
	double delta = 0.1;//does not affect accuracy, only affects the gap size for the half division method
	double *masOfRoots = new double[numOfVar];
	double norma = abs(normMatrix(matrix, numOfVar, numOfEq));
	double xForPositiveVal = INFINITY;
	double xForNegativeVal = -INFINITY;

	int rootCounter = 0;
	bool negativeValSearched = false;
	bool positiveValSearched = false;
	for (double i = -norma; i <= norma; i += delta) {
		double value = HessenbergCharacterPolyVal(i);
		if (value > 0) {
			xForPositiveVal = i;
			positiveValSearched = true;
		}
		if (value < 0) {
			xForNegativeVal = i;
			negativeValSearched = true;
		}
		if (positiveValSearched && negativeValSearched) {
			double root = BisectionMethod(HessenbergCharacterPolyVal, error, xForNegativeVal, xForPositiveVal);
			masOfRoots[rootCounter] = root;
			rootCounter++;
			if (value > 0) {
			   negativeValSearched = false;
			}
			else if(value < 0) {
			   positiveValSearched = false;
			}
		}
	}
	return masOfRoots;
}


void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	ofstream file;
	file.open("result.txt");

	file.precision(3);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file << matrix[i * numOfVariables + j] << "     ";
		}
		file << endl;
	}
	file << "--------------------------------------";
	file.close();
}


int main(void) {
	fileOut.open("matr1.txt");

	getMartixSize(&numOfVariables, &numOfEquations);

	//calculated the matrix of coefficients
	matrix = new double[numOfVariables * numOfEquations];
	getMatrix(numOfVariables, numOfEquations, matrix);

	fileOut.close();

	//--------------------------------------------------
	
	double *X = GetRoots(matrix, numOfVariables, numOfEquations);
	PrintMatrix(X, 1, numOfEquations);
}