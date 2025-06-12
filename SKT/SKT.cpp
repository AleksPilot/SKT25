#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
//#include "PTS_GEN.cpp"

using namespace std;

struct cell {
	double x = 0;
	double y = 0;
};


string line;
int ROWS = 6;
int COLS = 9;
int NUM = COLS; //количество приемников 


cell start;

//в данной задаче будем считать расположение приемников идентично клеткам, размер ячейки 1х1, MesOmegaK = 1;
vector<cell> receiver_positions; // = (i, 0) for 0 <= i <= NUM (см. ниже)
vector<vector<cell>> P; //P0
vector<vector<cell>> P_out; //вычисленные P

double coeff;
void Gen();
void inverse();
void direct(bool flag, vector<vector<cell>>P);
double Functional();
vector<cell> initialS; //сигналы в начальный момент времени
vector<cell> clS;
vector<cell> CalculatedS; //после вычисления p

vector<double> b; //b из уравнения

void writing()
{
	Gen();

	direct(false, P);
	inverse();
	direct(true, P_out);
	// Создаем или открываем файл для записи
	ofstream outFile("points.txt");

	if (!outFile.is_open()) {
		cerr << "Error in opening the file points.txt!" << endl;
	}

	// Записываем точки для первого графика (из P)
	for (size_t i = 0; i < P.size(); i++) {
		for (size_t j = 0; j < P[i].size(); j++) {
			outFile << P[i][j].x << " " << P[i][j].y << " ";
		}
		outFile << endl;
	}

	// Записываем разделитель
	outFile << "---" << endl;

	// Записываем точки для второго графика (из P_out)
	for (size_t i = 0; i < P_out.size(); i++) {
		for (size_t j = 0; j < P_out[i].size(); j++) {
			outFile << P_out[i][j].x << " " << P_out[i][j].y << " ";
		}
		outFile << endl;
	}

	// Закрываем файл
	outFile.close();

	cout << "Data in points.txt renewed!" << endl;

	//cout << Functional();


}


void direct(bool flag, vector<vector<cell>>P) {
	for (size_t i = 0; i < NUM; i++) 
	{
		clS[i].x = 0;
		clS[i].y = 0;
		for (size_t y = 0; y < ROWS; y++)
		{
			for (size_t x = 0; x < COLS; x++)
			{
				cell cell_center;
				cell_center.x = start.x + double((x + 0.5));
				cell_center.y = start.y - double((y + 0.5));
				
				//P ->-V
				double Local[2][2];
				//r = abs(vec(receiver, cell_center)), для формул расчёта L и S.
				double r = sqrt(pow((cell_center.x - receiver_positions[i].x), 2) + pow((cell_center.y - receiver_positions[i].y), 2));
				
				Local[0][0] = P[y][x].x * (((3 * pow((receiver_positions[i].x - cell_center.x), 2)) / (r * r)) - 1);
				Local[0][1] = P[y][x].y * ((3 * (receiver_positions[i].x - cell_center.x) * (receiver_positions[i].y - cell_center.y)) / (r * r));
				Local[1][0] = P[y][x].x * ((3 * (receiver_positions[i].x - cell_center.x) * (receiver_positions[i].y - cell_center.y)) / (r * r));
				Local[1][1] = P[y][x].y * (((3 * pow((receiver_positions[i].y - cell_center.y), 2)) / (r * r)) - 1);

				// так я делал все с допущением, что у нас все ячейки 1*1, а следовательно MesOmegaK = 1 => то формулы выше нужно будет переписывать
				clS[i].x += 1 * coeff / (r * r * r) * (Local[0][0] + Local[0][1]);
				clS[i].y += 1 * coeff / (r * r * r) * (Local[1][0] + Local[1][1]);
				CalculatedS[i].x = clS[i].x;
				CalculatedS[i].y = clS[i].y;
			}
		}
	}

	if (!flag) {
		swap(initialS, CalculatedS);
	}

}

void SLAE(vector<vector<double>>& A, vector<double>& b) { //GAUSS, не было времени использовать Pardiso/MFE
	int n = A.size();

	// Прямой ход (приведение к верхнетреугольному виду)
	for (int i = 0; i < n; i++) {
		// Поиск ведущего элемента в столбце (для устойчивости)
		int max_row = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > abs(A[max_row][i])) {
				max_row = k;
			}
		}

		// Обмен строк, если нужно
		if (max_row != i) {
			swap(A[i], A[max_row]);
			swap(b[i], b[max_row]);
		}

		// Проверка на нулевой диагональный элемент
		//if (abs(A[i][i]) < 1e-10) {
			//throw runtime_error("Matrix is singular or ill-conditioned");
		//}

		
		double div = A[i][i];
		for (int j = i; j < n; j++) {
			A[i][j] /= div;
		}
		b[i] /= div;

		
		for (int k = i + 1; k < n; k++) {
			double factor = A[k][i];
			b[k] -= factor * b[i];
			for (int j = i; j < n; j++) {
				A[k][j] -= factor * A[i][j];
			}
		}
	}

	
	for (int i = n - 1; i >= 0; i--) {
		for (int k = i - 1; k >= 0; k--) {
			double factor = A[k][i];
			b[k] -= factor * b[i];
			A[k][i] = 0; 
		}
	}
	int k = 0;
	for (size_t x = 0; x < ROWS; x++) {
		for (size_t y = 0; y < COLS; y++)
		{
			P_out[x][y].x = b[k];
			P_out[x][y].y = b[k + 1];
			k += 2;
		}
	}
}


void inverse() {
	vector<vector<double>> Glob(2 * NUM, vector<double>(2 * ROWS * COLS));
	vector<vector<double>> transposedGlob(2 * ROWS * COLS, vector<double>(2 * NUM));

	for (size_t i = 0; i < NUM; i++)
	{
		int k = 0;
		for (size_t y = 0; y < ROWS; y++)
		{
			for (size_t x = 0; x < COLS; x++)
			{
				cell cell_center;
				cell_center.x = start.x + double((x + 0.5));
				cell_center.y = start.y - double((y + 0.5));
				
				//P ->-V
				double Local[2][2];
				//r = abs(vec(receiver, cell_center))
				double r = sqrt(pow((cell_center.x - receiver_positions[i].x), 2) + pow((cell_center.y - receiver_positions[i].y), 2));

				Local[0][0] = (((3 * pow((receiver_positions[i].x - cell_center.x), 2)) / (r * r)) - 1);
				Local[0][1] = ((3 * (receiver_positions[i].x - cell_center.x) * (receiver_positions[i].y - cell_center.y)) / (r * r));
				Local[1][0] = ((3 * (receiver_positions[i].x - cell_center.x) * (receiver_positions[i].y - cell_center.y)) / (r * r));
				Local[1][1] = (((3 * pow((receiver_positions[i].y - cell_center.y), 2)) / (r * r)) - 1);

				Glob[i * 2][k] = 1 * coeff / (r * r * r) * Local[0][0];
				Glob[i * 2 + 1][k] = 1 * coeff / (r * r * r) * Local[1][0];
				Glob[i * 2][k + 1] = 1 * coeff / (r * r * r) * Local[0][1];
				Glob[i * 2 + 1][k + 1] = 1 * coeff / (r * r * r) * Local[1][1];
				k += 2;
			}
		}
	}

	//транспонирование матрицы (Lt)
	for (size_t i = 0; i < Glob.size(); i++) {
		for (size_t j = 0; j < Glob[0].size(); j++) {
			transposedGlob[j][i] = Glob[i][j];
		}
	}

	//произведение матриц L*Lt
	int maxDim = 2 * COLS * ROWS;
	int minDim = 2 * NUM;

	vector<vector<double>> A(maxDim, vector<double>(maxDim));
	for (int i = 0; i < A.size(); i++)
	{
		fill(A[i].begin(), A[i].end(), 0);
	}

	for (size_t i = 0; i < maxDim; i++) {
		for (size_t j = 0; j < maxDim; j++) {
			for (size_t k = 0; k < minDim; k++)
			{
				A[i][j] += transposedGlob[i][k] * Glob[k][j];
			}
		}
	}

	for (size_t i = 0; i < maxDim; i++) {
		for (size_t j = 0; j < 1; j++) {
			for (size_t k = 0; k < minDim; k++)
			{
				if (k % 2 == 0)
					b[i] += transposedGlob[i][k] * clS[k / 2].x;
				else
					b[i] += transposedGlob[i][k] * clS[k / 2].y;
			}
		}
	}

	cout;
	SLAE(A, b);
}

double Functional() {
	
	double Sum = 0;
	for (size_t i = 0; i < CalculatedS.size(); i++)
	{
		Sum += pow((initialS[i].x - CalculatedS[i].x), 2);
		Sum += pow((initialS[i].y - CalculatedS[i].y), 2);
	}
	return Sum;
}


void Gen() { //value - значение намагниченности в неоднородности
	
	receiver_positions.resize(NUM);
	for (int i = 0; i < NUM; i++)
	{
		receiver_positions[i].x = i;
		receiver_positions[i].y = 0;

	}
	start = { 0, 0 };
	ROWS = 6;
	COLS = 9;
	NUM = COLS;
	double I = 1;
	P.resize(ROWS);
	P_out.resize(ROWS);
	for (size_t x = 0; x < ROWS; x++)
	{
		P[x].resize(COLS);
		P_out[x].resize(COLS);
		for (size_t y = 0; y < COLS; y++)
		{
			if (x >= 1 && x <= 2 && y >= 1 && y <= 2) {
				P[x][y].x = 1;
				P[x][y].y = 1;
			}
			else {
				P[x][y].x = 1e-15;
				P[x][y].y = 1e-15;
			}
		}
	}
	
	coeff = I / (4 * M_PI);
	initialS.resize(NUM); //сигналы в начальный момент времени
	clS.resize(NUM);
	CalculatedS.resize(NUM); //после вычисления p
	b.resize(2 * ROWS * COLS); //
}



int main()
{
	Gen();

	direct(false, P);
	inverse();
	direct(true, P_out);
	for (size_t i = 0; i < P.size(); i++)
	{
		for (size_t j = 0; j < P[i].size(); j++)
		{
			cout << i << " " << j << "     " << P[i][j].x << " " << P[i][j].y << endl;
		}

	}
	cout << endl;
	for (size_t i = 0; i < P_out.size(); i++)
	{
		for (size_t j = 0; j < P_out[i].size(); j++)
		{
			cout << i << " " << j << "     " << P_out[i][j].x << " " << P_out[i][j].y << endl;
		}
	}
	cout << Functional();
	writing();
}