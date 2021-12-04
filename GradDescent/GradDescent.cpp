#include <iostream>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

class Algorithm {
public:
    // ����������� ��� ���������� �������
    int M = 10;
    vector<vector<double>> A;
    vector<double> x, b, r;
    vector<double> F, F_prime;
    int cnt = 0; // ��� ����������

    Algorithm(vector<vector<double>> A, vector<double> x, vector<double> b, int M) {
        this->A = A;
        this->x = x;
        this->b = b;
        this->M = M;
    }

    // ��������� ����� ����� ����� ����������� �� ������������, ����� ������� ������� ������ ������� - ��� �����
    double norm() {
        F_prime = F_derivative_function();
        F = F_prime;
        return sqrt(vector_multiplication(F_prime, F_prime));
    }

    void update() {
        ++cnt;
        // ���������� ������ ������� ��������� �������
        // ����� ���� �� ����������� ����� ������, �� � �������� ������ ������, ������ ������ ������ + ����
        this->x = get_next_x();
    }

private:
    // ������ ������ � (��� � r) 1 � mm ������� ���������� ��� ��� ������, � �� ������ �������� 
    vector <double> multiply_matricies(vector<vector<double>> A, vector<double> x) {
        vector<double> result(M);
        for (int i = 0; i < M; ++i) {
            double sum_of_row = 0;
            for (int j = 0; j < M; ++j) {
                sum_of_row += A[i][j] * x[j];
            }
            result[i] = sum_of_row;
        }
        return result;
    }

    double vector_multiplication(vector<double> x, vector<double> y) {
        double res = 0;
        for (int i = 0, size = x.size(); i < size; ++i)
            res += x[i] * y[i];
        return res;
    }

    vector<double> F_derivative_function() {
        vector<double> res(M);
        vector<double> A_by_X = multiply_matricies(A, x);
        for (int i = 0; i < M; ++i)
            res[i] = A_by_X[i] - b[i];
        return res;
    }

    vector<double> inverse(vector<double> vec) {
        for (int i = 0; i < M; ++i) vec[i] *= -1;
        return vec;
    }

    vector<double> getR() {
        return inverse(F_derivative_function());
    }

    vector<double> get_next_x() {
        // �������� alpha = <Rk, Rk> / <Rk, ARk>
        this->r = getR();
        vector<double> Rk = this->r;
        double alpha = (vector_multiplication(Rk, Rk)) / (vector_multiplication(Rk, multiply_matricies(A, Rk)));
        vector<double> F_prime = F_derivative_function();
        vector<double> new_x(M);

        for (int i = 0; i < M; ++i) {
            // ������� ��� ���������� X[k+1]
            //F_prime[i] = x[i] - F_prime[i] * alpha;
            new_x[i] = x[i] - F_prime[i] * alpha;
        }
        return new_x;
    }

    // F(x) = 1/2 * <Ax, x> - <b, x>
    double F_function_calculate() {
        // ������� �������� � �� �. ��������� � ������ 1 � M - ������� ������
        vector<double> _Axb_ = multiply_matricies(A, x);
        double first_term = vector_multiplication(_Axb_, x) * 0.5;
        double second_term = vector_multiplication(b, x);
        return first_term - second_term;
    }
    // ��� ���������� �� ������������ (?)

};

int main()
{
    system("chcp 1251");
    int M = 0;
    vector<vector<double>> vec(M);
    vector<double> row(M);
    vector<double> b(M);
    double eps = 0.00001;

    cout << "������� ����������� ������� � � �";
    cin >> M;
    
    cout << "��������������� ������� �������� ������� (1 ������, ����� 2 ������ � ��)";
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            cin >> row[j];
        }
        vec[i] = row;
    }

    cout << "������� �������� b[i]";
    for (int i = 0; i < M; ++i)
        cin >> b[i];

    cout << "������� �������� eps: ";
    cin >> eps;
    
    // ������������, ��� ���������� ����� � � = ���-�� ����� / �������� � �������
    vector<double> x(M);
    cout << "������� ��������������� �������� ������� �[0]";
    for (int i = 0; i < M; ++i) cin >> x[i];

    // �������� ��������
    Algorithm algorithm(vec, x, b, M);
    for (int k = 0; k < 1000 // ����� ����������� ������ ��� ������������. � ���������� k < 1000 ������ ���� ������ � ������������, ���� x >= eps
        ; ++k) {

        if (algorithm.norm() < eps) break;
        algorithm.update();
    }
    cout << "����������� �������:\n";
    for (auto sol : algorithm.x)
        cout << sol << "\n";
    
}
