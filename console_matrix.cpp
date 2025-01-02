#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <windows.h>
#include <cmath>
#include "best_file.cpp"

using namespace std;

void write_matrix(vector<vector<long double>>& matrix, int& n, int& m) {
	for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
	}
}

int main() {
    int testing_func;
    cout << "Function(cin number):"
    << endl << "1) Matrix trace"
    << endl << "2) Element of matrix"
    << endl << "3) Matrix addition"
    << endl << "4) Matrix multiplication by scalar"
    << endl << "5) Matrix multiplication"
    << endl << "6) Matrix determinant" << endl;
    cin >> testing_func;

    cout << "Cin matrix:" << endl;
    Matrix matrix1 = Matrix(0, 0, {}, {}, {}, "");

    if (testing_func == 1) {
        cout << "Result:" << endl;
        long double answer = matrix1.get_trace();
        cout << answer << endl;
    }else if (testing_func == 2) {
        int i, j;
        cout << "Cin i, j" << endl;
        cin >> i >> j;

        long double answer = matrix1.get_element(i, j);
        cout << "Result:" << endl;
        cout << answer;
    }else if (testing_func == 3) {
        cout << "Cin second Matrix:" << endl;
        Matrix matrix2 = Matrix(0, 0, {}, {}, {}, "");
        int n = matrix1.n;

        vector<vector<long double>> answer = sum_of_matrix(matrix1, matrix2).to_vector();
        cout << "Result:" << endl;
        write_matrix(answer, n, n);
    }else if (testing_func == 4) {
        long double rd_scalar;
        cout << "Cin scalar: ";
        cin >> rd_scalar;
        cout << endl;

        int n = matrix1.n;
        multiply_scalar(rd_scalar, matrix1);
        vector<vector<long double>> answer = matrix1.to_vector();
        cout << "Result:" << endl;
        write_matrix(answer, n, n);
    }else if (testing_func == 5) {
        cout << "Cin second Matrix:" << endl;
        Matrix matrix2 = Matrix(0, 0, {}, {}, {}, "");
        int n = matrix1.n;

        vector<vector<long double>> answer = multiply_matrix(matrix1, matrix2).to_vector();
        cout << "Result:" << endl;
        write_matrix(answer, n, n);
    }else {
        pair<long double, string> answer = get_determinant(matrix1);

        cout << "Result:" << endl;
        cout << answer.first << " " << answer.second;
    }
}
