#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<cctype>
#include <tuple> 
#include <iomanip> 
#include <cmath>
#include <algorithm>
using namespace std;

void printMatrix(const vector<vector<float>>& A, int N) {
    cout << fixed << setprecision(6);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "----------------------" << endl;
}

void printMatrixLU(const vector<vector<float>>& matrix, int N, const string& name) {
    cout << "\nMacierz " << name << ":" << endl;
    cout << fixed << setprecision(6);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(12) << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printVector(const vector<float>& vec, const string& name) {
    cout << "\nWektor " << name << ":" << endl;
    cout << fixed << setprecision(6);
    for (size_t i = 0; i < vec.size(); i++) {
        cout << setw(12) << vec[i] << endl;
    }
    cout << endl;
}

tuple<vector<vector<float>>, vector<float>, int> read_file(string file) {
    ifstream myfile(file);
    if (!myfile.is_open()) {
        cerr << "Nie mozna otworzyc pliku!" << endl;
        return make_tuple(vector<vector<float>>(), vector<float>(), 0);
    }

    vector<vector<float>> A;
    vector<float> b;
    int N = 0;
    string line;
    string tmp;
    float value;

    cout << "Plik otwarty, wczytywanie danych..." << endl;

    while (getline(myfile, line)) {
        stringstream ss(line);
        ss >> tmp;

        if (tmp == "N") {
            ss >> tmp;  // Skip '='
            ss >> N;    // Wczytaj wartoœæ N
            cout << "Wczytano N = " << N << endl;
            break;
        }
    }

    bool reading_b = false;
    while (getline(myfile, line)) {
        stringstream ss(line);
        ss >> tmp;

        if (tmp == "b:") {
            reading_b = true;
            continue;
        }

        if (reading_b) {
            stringstream line_stream(line);
            while (line_stream >> value) {
                b.push_back(value);
            }

            cout << "Wczytano b: ";
            for (size_t i = 0; i < b.size(); i++) {
                cout << b[i] << " ";
            }
            cout << endl;
            break;
        }
    }

    bool found_A = false;
    while (getline(myfile, line)) {
        stringstream ss(line);
        ss >> tmp;

        if (tmp == "A:") {
            found_A = true;
            continue;
        }

        if (found_A) {
            vector<float> row;
            stringstream line_stream(line);
            while (line_stream >> value) {
                row.push_back(value);
            }
            if (!row.empty()) {
                A.push_back(row);
            }
        }
    }

    cout << "N = " << N << endl;
    cout << "Wczytano macierz A o rozmiarze: " << A.size() << "x" << A[0].size() << endl;
    cout << "Wczytano wektor b o rozmiarze: " << b.size() << endl;

    return make_tuple(A, b, N);
}

tuple<vector<vector<float>>, vector<vector<float>>, vector<int>> luDecomposition(vector<vector<float>> A, int N) {
    vector<vector<float>> L(N, vector<float>(N, 0.0f));
    vector<vector<float>> U = A; 
    vector<int> P(N);

    for (int i = 0; i < N; i++) {
        P[i] = i;
    }

    for (int i = 0; i < N; i++) {
        L[i][i] = 1.0f;
    }

    cout << "\n=== ROZKLAD LU - KOLEJNE ITERACJE ===" << endl;

    for (int k = 0; k < N - 1; k++) {
        cout << "\n--- Iteracja " << (k + 1) << " ---" << endl;

        int max_row = k;
        float max_val = abs(U[k][k]);

        for (int i = k + 1; i < N; i++) {
            if (abs(U[i][k]) > max_val) {
                max_val = abs(U[i][k]);
                max_row = i;
            }
        }

        if (max_row != k) {
            cout << "Zamiana wierszy " << k << " i " << max_row << endl;
            swap(U[k], U[max_row]);
            swap(P[k], P[max_row]);

            for (int j = 0; j < k; j++) {
                swap(L[k][j], L[max_row][j]);
            }
        }

        if (abs(U[k][k]) < 1e-10f) {
            cout << "Uwaga: Element glowny bardzo maly: " << U[k][k] << endl;
        }

        for (int i = k + 1; i < N; i++) {
            float multiplier = U[i][k] / U[k][k];
            L[i][k] = multiplier;

            for (int j = k; j < N; j++) {
                U[i][j] = U[i][j] - multiplier * U[k][j];
            }
        }

        cout << "Po iteracji " << (k + 1) << ":" << endl;
        printMatrixLU(L, N, "L");
        printMatrixLU(U, N, "U");
    }

    return make_tuple(L, U, P);
}

// Rozwiazywanie Lz = Pb metoda forward substitution
vector<float> forwardSubstitution(const vector<vector<float>>& L, const vector<float>& Pb, int N) {
    vector<float> z(N, 0.0f);

    for (int i = 0; i < N; i++) {
        z[i] = Pb[i];
        for (int j = 0; j < i; j++) {
            z[i] -= L[i][j] * z[j];
        }
        z[i] /= L[i][i];
    }

    return z;
}

// Rozwiazywanie Ux = z metoda backward substitution
vector<float> backwardSubstitution(const vector<vector<float>>& U, const vector<float>& z, int N) {
    vector<float> x(N, 0.0f);

    for (int i = N - 1; i >= 0; i--) {
        x[i] = z[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

// Sprawdzenie poprawnosci rozwiazania
void checkSolution(const vector<vector<float>>& A, const vector<float>& x, const vector<float>& b, int N) {
    cout << "\n=== SPRAWDZENIE POPRAWNOSCI ROZWIAZANIA ===" << endl;
    vector<float> Ax(N, 0.0f);

    // Oblicz A*x
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    // Porownaj z b
    cout << "Porownanie A*x z b:" << endl;
    cout << fixed << setprecision(6);
    cout << setw(5) << "i" << setw(15) << "A*x[i]" << setw(15) << "b[i]" << setw(15) << "roznica" << endl;
    cout << string(50, '-') << endl;

    float max_error = 0.0f;
    for (int i = 0; i < N; i++) {
        float error = abs(Ax[i] - b[i]);
        max_error = max(max_error, error);
        cout << setw(5) << i << setw(15) << Ax[i] << setw(15) << b[i] << setw(15) << error << endl;
    }
    cout << "\nMaksymalny blad: " << max_error << endl;
}

int main() {
    vector<string> file_names = {"LU_gr1INO 1.txt", "gauss_elimination_gr1IO_A.txt", "gauss_elimination_gr1IO_B.txt", "gauss_elimination_gr1IO_C.txt"};

    for (const string file_name : file_names) {
        // Wczytaj dane
        vector<vector<float>> A;
        vector<float> b;
        int N;

        tie(A, b, N) = read_file(file_name);

        if (N == 0) {
            cerr << "Blad wczytywania danych!" << endl;
            return -1;
        }

        cout << "\n=== DANE WEJŒCIOWE ===" << endl;
        cout << "Rozmiar ukladu: " << N << "x" << N << endl;

        // Wyswietl macierz A i wektor b
        cout << "\nMacierz wspolczynnikow A:" << endl;
        printMatrixLU(A, N, "A");

        cout << "Wektor wyrazow wolnych b:" << endl;
        for (int i = 0; i < N; i++) {
            cout << setw(8) << b[i] << " ";
        }
        cout << endl << endl;

        // Wykonaj rozklad LU
        vector<vector<float>> L, U;
        vector<int> P;

        tie(L, U, P) = luDecomposition(A, N);

        cout << "\n=== WYNIKI ROZKLADU LU ===" << endl;
        printMatrixLU(L, N, "L (finalna)");
        printMatrixLU(U, N, "U (finalna)");

        // Zastosuj permutacje do wektora b
        vector<float> Pb(N);
        for (int i = 0; i < N; i++) {
            Pb[i] = b[P[i]];
        }

        cout << "Wektor permutacji P: ";
        for (int i = 0; i < N; i++) {
            cout << P[i] << " ";
        }
        cout << endl;

        printVector(Pb, "Pb (po permutacji)");

        // Rozwiaz Lz = Pb
        vector<float> z = forwardSubstitution(L, Pb, N);
        printVector(z, "z (rozwiazanie Lz = Pb)");

        // Rozwiaz Ux = z
        vector<float> x = backwardSubstitution(U, z, N);
        printVector(x, "x (rozwiazanie finalne)");

        // Sprawdz poprawnosc
        checkSolution(A, x, b, N);
    }
    

    return 0;
}