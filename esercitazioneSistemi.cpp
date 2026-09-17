#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace std;

static const int MAXN = 120; // sufficiente per n=34 (e margine)

// Utility stampa 
static void header(const char* s) {
    cout << "\n==================== " << s << " ====================\n";
}

// Norme 
static float normaInfMatrice(const float A[MAXN][MAXN], int n) {
    float mx = 0.0f;
    for (int i = 0; i < n; ++i) {
        float s = 0.0f;
        for (int j = 0; j < n; ++j) s += (float)fabs(A[i][j]);
        if (s > mx) mx = s;
    }
    return mx;
}

static float normaInfVettore(const float v[MAXN], int n) {
    float mx = 0.0f;
    for (int i = 0; i < n; ++i) {
        float a = (float)fabs(v[i]);
        if (a > mx) mx = a;
    }
    return mx;
}

static float errInf(const float x[MAXN], const float y[MAXN], int n) {
    float mx = 0.0f;
    for (int i = 0; i < n; ++i) {
        float d = (float)fabs(x[i] - y[i]);
        if (d > mx) mx = d;
    }
    return mx;
}

// Copie 
static void copiaMatrice(const float A[MAXN][MAXN], float B[MAXN][MAXN], int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            B[i][j] = A[i][j];
}

static void copiaVettore(const float a[MAXN], float b[MAXN], int n) {
    for (int i = 0; i < n; ++i) b[i] = a[i];
}

// b = A * xbar con xbar = (1,...,1)^T 
static void calcolaTermineNoto(const float A[MAXN][MAXN], float b[MAXN], int n) {
    for (int i = 0; i < n; ++i) {
        float s = 0.0f;
        for (int j = 0; j < n; ++j) s += A[i][j]; // somma riga (xbar=1)
        b[i] = s;
    }
}

// Perturbazione: delta b = ||b||inf * (-0.01, 0.01, -0.01, ...) 
static void perturbaB(const float b[MAXN], float btilde[MAXN], int n) {
    float nb = normaInfVettore(b, n);
    for (int i = 0; i < n; ++i) {
        float segno = (i % 2 == 0) ? -0.01f : 0.01f;
        btilde[i] = b[i] + nb * segno;
    }
}

// Residuo r = b - A x (norma infinito) 
static float residuoInf(const float A[MAXN][MAXN], const float x[MAXN], const float b[MAXN], int n) {
    float mx = 0.0f;
    for (int i = 0; i < n; ++i) {
        float Ax = 0.0f;
        for (int j = 0; j < n; ++j) Ax += A[i][j] * x[j];
        float ri = b[i] - Ax;
        float ari = (float)fabs(ri);
        if (ari > mx) mx = ari;
    }
    return mx;
}

// Gauss con pivoting parziale 
static void swapRows(float A[MAXN][MAXN], float b[MAXN], int r1, int r2, int n) {
    if (r1 == r2) return;
    for (int j = 0; j < n; ++j) {
        float tmp = A[r1][j];
        A[r1][j] = A[r2][j];
        A[r2][j] = tmp;
    }
    float tb = b[r1];
    b[r1] = b[r2];
    b[r2] = tb;
}

static bool gaussPivot(float A[MAXN][MAXN], float b[MAXN], float x[MAXN], int n) {
    // Eliminazione
    for (int k = 0; k < n - 1; ++k) {
        int piv = k;
        float mx = (float)fabs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            float v = (float)fabs(A[i][k]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (mx == 0.0f) return false;
        swapRows(A, b, k, piv, n);

        for (int i = k + 1; i < n; ++i) {
            float m = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }
    if ((float)fabs(A[n-1][n-1]) == 0.0f) return false;

    // Sostituzione all’indietro
    for (int i = n - 1; i >= 0; --i) {
        float s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return true;
}

// Costruzione matrici 
static void costruisciA1(float A[MAXN][MAXN]) {
    for (int i=0;i<MAXN;++i) for (int j=0;j<MAXN;++j) A[i][j]=0.0f;
    A[0][0]= 3; A[0][1]= 1; A[0][2]=-1; A[0][3]= 0;
    A[1][0]= 0; A[1][1]= 7; A[1][2]=-3; A[1][3]= 0;
    A[2][0]= 0; A[2][1]=-3; A[2][2]= 9; A[2][3]=-2;
    A[3][0]= 0; A[3][1]= 0; A[3][2]= 4; A[3][3]=-10;
}

static void costruisciA2(float A[MAXN][MAXN]) {
    for (int i=0;i<MAXN;++i) for (int j=0;j<MAXN;++j) A[i][j]=0.0f;
    A[0][0]= 2; A[0][1]= 4; A[0][2]=-2; A[0][3]= 0;
    A[1][0]= 1; A[1][1]= 3; A[1][2]= 0; A[1][3]= 1;
    A[2][0]= 3; A[2][1]=-1; A[2][2]= 1; A[2][3]= 2;
    A[3][0]= 0; A[3][1]=-1; A[3][2]= 2; A[3][3]= 1;
}

// Pascal n x n: recurrence (0-based)
static void costruisciPascal(float A[MAXN][MAXN], int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = 0.0f;

    for (int i = 0; i < n; ++i) A[i][0] = 1.0f;
    for (int j = 0; j < n; ++j) A[0][j] = 1.0f;

    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            A[i][j] = A[i-1][j] + A[i][j-1];
}

// Tridiagonale: 2 diag, -1 sotto/sopra
static void costruisciTridiagonale(float A[MAXN][MAXN], int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) A[i][j] = 2.0f;
            else if (abs(i - j) == 1) A[i][j] = -1.0f;
            else A[i][j] = 0.0f;
        }
    }
}

// Test su una matrice 
static void runMatrix(const char* nome, const float Aorig[MAXN][MAXN], int n) {
    float x_atteso[MAXN];
    for (int i = 0; i < n; ++i) x_atteso[i] = 1.0f;

    float normaA = normaInfMatrice(Aorig, n);

    // b = A*xbar
    float b0[MAXN];
    calcolaTermineNoto(Aorig, b0, n);
    float normaB = normaInfVettore(b0, n);

    // Risoluzione Ax=b0
    float A1[MAXN][MAXN], b1[MAXN], x[MAXN];
    copiaMatrice(Aorig, A1, n);
    copiaVettore(b0, b1, n);

    bool ok = gaussPivot(A1, b1, x, n);

    cout << "\n--- " << nome << " (n=" << n << ") ---\n";
    cout << "||A||_inf = " << (double)normaA << "\n";
    cout << "||b||_inf = " << (double)normaB << "\n";

    if (!ok) {
        cout << "Gauss (pivoting) fallito: pivot nullo / matrice singolare.\n";
        return;
    }

    float errSol = errInf(x, x_atteso, n);
    float res = residuoInf(Aorig, x, b0, n);

    cout << "Errore ||x - x_atteso||_inf = " << (double)errSol << "\n";
    cout << "Residuo ||b - A x||_inf      = " << (double)res << "\n";

    // Perturbazione b~
    float btilde[MAXN];
    perturbaB(b0, btilde, n);

    // delta b e sua norma
    float db[MAXN];
    for (int i = 0; i < n; ++i) db[i] = btilde[i] - b0[i];
    float normDb = normaInfVettore(db, n);

    // Risoluzione A x~ = b~
    float A2[MAXN][MAXN], b2[MAXN], xtilde[MAXN];
    copiaMatrice(Aorig, A2, n);
    copiaVettore(btilde, b2, n);

    ok = gaussPivot(A2, b2, xtilde, n);
    if (!ok) {
        cout << "Gauss (pivoting) su (b+delta b) fallito.\n";
        return;
    }

    float errPert = errInf(x, xtilde, n);
    float res2 = residuoInf(Aorig, xtilde, btilde, n);

    cout << "||delta b||_inf               = " << (double)normDb << "\n";
    cout << "Errore ||x - x~||_inf          = " << (double)errPert << "\n";
    cout << "Residuo ||b~ - A x~||_inf       = " << (double)res2 << "\n";

    // Fattori di amplificazione
    float alpha = 0.0f;
    if (normDb != 0.0f) alpha = errPert / normDb;

    float normX = normaInfVettore(x, n);
    float beta = 0.0f;
    if (normDb != 0.0f && normaB != 0.0f && normX != 0.0f) {
        beta = (errPert / normX) / (normDb / normaB);
    }

    cout << "Amplificazione assoluta  alpha = ||x-x~|| / ||db|| = " << (double)alpha << "\n";
    cout << "Amplificazione relativa   beta = (||x-x~||/||x||) / (||db||/||b||) = " << (double)beta << "\n";
}

int main() {
    // Matricola 5669424 => d1=2, d0=4 => nT = 10*(d1+1)+d0 = 34
    const int d0 = 4;
    const int d1 = 2;
    const int nT = 10 * (d1 + 1) + d0; // 34
    const int nP = 10;                // Pascal: n=10 

    cout << scientific << setprecision(16);

    header("SISTEMI LINEARI");

    float A1[MAXN][MAXN], A2[MAXN][MAXN], P[MAXN][MAXN], T[MAXN][MAXN];
    costruisciA1(A1);
    costruisciA2(A2);
    costruisciPascal(P, nP);
    costruisciTridiagonale(T, nT);

    header("RISULTATI");
    runMatrix("Matrice A1 (4x4)", A1, 4);
    runMatrix("Matrice A2 (4x4)", A2, 4);
    runMatrix("Matrice Pascal",   P,  nP);
    runMatrix("Matrice Tridiagonale", T, nT);

    cout << "\nFine.\n";
    return 0;
}

