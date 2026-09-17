#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace std;

//  Funzioni di supporto
static void header(const char* s) {
    cout << "\n--- " << s << " ---\n";
}

static void stampa_err(double approx, double esatto) {
    double ea = fabs(approx - esatto);
    cout << "approx=" << approx << "  exact=" << esatto
         << "  abs_err=" << ea;
    if (esatto != 0.0) cout << "  rel_err=" << ea / fabs(esatto);
    cout << "\n";
}

// ESERCIZIO 1 
void esercizio1() {
    int d0 = 4; // ultima cifra matricola
    int d1 = 2; // penultima cifra matricola

    // b molto grande per provocare perdita di significato
    double b = (d1 + 1) * pow(10.0, 20.0);
    double c = -b;

    header("ESERCIZIO 1");
    cout << scientific << setprecision(16);

    for (int i = 0; i <= 6; ++i) {
        double a = (d0 + 1) * pow(10.0, (double)i);

        double res1 = (a + b) + c;
        double res2 = a + (b + c); // qui b+c fa 0

        double diff = fabs(res1 - res2);

        cout << "i=" << i << "  a=" << a << "\n";
        cout << " (a+b)+c = " << res1 << "\n";
        cout << " a+(b+c) = " << res2 << "\n";
        cout << " diff    = " << diff;

        // errore relativo rispetto al valore corretto (a)
        if (a != 0.0) cout << "  (rel wrt a) = " << diff / fabs(a);
        cout << "\n\n";
    }
}

//  ESERCIZIO 2 
// Taylor di exp(x) di grado N
double taylor_exp(double x, int N) {
    double somma = 1.0, termine = 1.0;
    for (int n = 1; n <= N; ++n) {
        termine *= x / (double)n;
        somma += termine;
    }
    return somma;
}

void esercizio2() {
    header("ESERCIZIO 2");
    cout << scientific << setprecision(16);

    double x_values[4] = {0.5, 30.0, -0.5, -30.0};
    int N_values[5] = {3, 10, 50, 100, 150};

    // ---- Algoritmo 1 ----
    cout << "ALGORITMO 1 (Taylor diretto)\n";
    for (int ix = 0; ix < 4; ++ix) {
        double x = x_values[ix];
        cout << "\nx = " << x << "\n";
        for (int in = 0; in < 5; ++in) {
            int N = N_values[in];
            double approx = taylor_exp(x, N);
            double esatto = exp(x);

            cout << " N=" << setw(3) << N << "  ";
            stampa_err(approx, esatto);
        }
    }

    // ---- Algoritmo 2 ----
    // exp(-x) = 1/exp(x), approssimando exp(x) con Taylor su +x
    cout << "\nALGORITMO 2 (exp(-x) approx 1 / fN(x))\n";
    double x_pos[2] = {0.5, 30.0};
    for (int ix = 0; ix < 2; ++ix) {
        double x = x_pos[ix];
        cout << "\nx = -" << x << "\n";
        for (int in = 0; in < 5; ++in) {
            int N = N_values[in];
            double denom = taylor_exp(x, N);     // approx exp(x)
            double approx;

            // evito divisioni per numeri troppo piccoli
            if (fabs(denom) < 1e-300) approx = numeric_limits<double>::infinity();
            else approx = 1.0 / denom;

            double esatto = exp(-x);

            cout << " N=" << setw(3) << N << "  ";
            stampa_err(approx, esatto);
        }
    }
}

// ESERCIZIO 3 
void esercizio3() {
    header("ESERCIZIO 3");
    cout << scientific << setprecision(16);

    // float
    float eps_f = 1.0f;
    int d_f = 0;
    while (1.0f + eps_f / 2.0f > 1.0f) {
        eps_f /= 2.0f;
        d_f++;
    }

    // double
    double eps_d = 1.0;
    int d_d = 0;
    while (1.0 + eps_d / 2.0 > 1.0) {
        eps_d /= 2.0;
        d_d++;
    }

    cout << "float : eps=" << eps_f << "  d=" << d_f
         << "  std_epsilon=" << numeric_limits<float>::epsilon() << "\n";
    cout << "double: eps=" << eps_d << "  d=" << d_d
         << "  std_epsilon=" << numeric_limits<double>::epsilon() << "\n";
}

int main() {
    int scelta = 0;
    do {
        cout << "\n--- MENU ---\n";
        cout << "1. Esercizio 1 (aritmetica di macchina)\n";
        cout << "2. Esercizio 2 (Taylor exp e errori)\n";
        cout << "3. Esercizio 3 (epsilon di macchina)\n";
        cout << "0. Esci\n";
        cout << "Scegli un esercizio: ";
        cin >> scelta;

        switch (scelta) {
            case 1: esercizio1(); break;
            case 2: esercizio2(); break;
            case 3: esercizio3(); break;
            case 0: cout << "Uscita...\n"; break;
            default: cout << "Scelta non valida!\n";
        }
    } while (scelta != 0);

    return 0;
}

