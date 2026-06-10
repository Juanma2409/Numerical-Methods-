#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

// ============================================================
// Funcion del problema
// f(x) = x*sin(x^2)
// ============================================================
long double f(long double x) {
    return x * sin(x * x);
}

// ============================================================
// Derivadas exactas
// f'(x)  = sin(x^2) + 2x^2 cos(x^2)
// f''(x) = 6x cos(x^2) - 4x^3 sin(x^2)
// ============================================================
long double df_exacta(long double x) {
    return sin(x * x) + 2.0L * x * x * cos(x * x);
}

long double d2f_exacta(long double x) {
    return 6.0L * x * cos(x * x) - 4.0L * x * x * x * sin(x * x);
}

// ============================================================
// Diferencias finitas
// ============================================================
long double diferencia_adelante(long double (*func)(long double), long double x, long double h) {
    return (func(x + h) - func(x)) / h;
}

long double diferencia_atras(long double (*func)(long double), long double x, long double h) {
    return (func(x) - func(x - h)) / h;
}

long double diferencia_centrada(long double (*func)(long double), long double x, long double h) {
    return (func(x + h) - func(x - h)) / (2.0L * h);
}

long double segunda_derivada(long double (*func)(long double), long double x, long double h) {
    return (func(x + h) - 2.0L * func(x) + func(x - h)) / (h * h);
}

// ============================================================
// Utilidades de salida
// ============================================================
void linea(char c = '=', int n = 60) {
    cout << string(n, c) << '\n';
}

void mostrar_valor(const string& etiqueta, long double valor) {
    cout << left << setw(28) << etiqueta << ": "
         << right << setw(20) << valor << '\n';
}

int main() {
    cout << fixed << setprecision(15);

    // Datos del problema
    long double x0 = 0.0L;
    long double h  = 0.01L;

    // Valores exactos
    long double exacta_primera = df_exacta(x0);
    long double exacta_segunda = d2f_exacta(x0);

    // Valores de la función evaluados
    long double fxmh = f(x0 - h);
    long double fx0   = f(x0);
    long double fxph  = f(x0 + h);

    // Aproximaciones
    long double D_adelante = diferencia_adelante(f, x0, h);
    long double D_atras    = diferencia_atras(f, x0, h);
    long double D_centrada = diferencia_centrada(f, x0, h);
    long double D2         = segunda_derivada(f, x0, h);

    // Errores absolutos
    long double err_adelante = fabsl(exacta_primera - D_adelante);
    long double err_atras    = fabsl(exacta_primera - D_atras);
    long double err_centrada = fabsl(exacta_primera - D_centrada);
    long double err_D2       = fabsl(exacta_segunda - D2);

    // Salida
    linea('=');
    cout << "      DIFERENCIAS FINITAS - METODOS NUMERICOS\n";
    linea('=');

    cout << "\nFuncion analizada:\n";
    cout << "f(x) = x*sin(x^2)\n";

    cout << "\nDatos del problema:\n";
    mostrar_valor("Punto de evaluacion x0", x0);
    mostrar_valor("Paso h", h);

    linea('-');
    cout << "VALORES DE LA FUNCION\n";
    linea('-');
    mostrar_valor("f(x0 - h)", fxmh);
    mostrar_valor("f(x0)", fx0);
    mostrar_valor("f(x0 + h)", fxph);

    linea('-');
    cout << "DERIVADAS EXACTAS\n";
    linea('-');
    mostrar_valor("f'(x0)", exacta_primera);
    mostrar_valor("f''(x0)", exacta_segunda);

    linea('=');
    cout << "PRIMERA DERIVADA\n";
    linea('=');
    mostrar_valor("Diferencia hacia adelante", D_adelante);
    mostrar_valor("Error absoluto", err_adelante);
    cout << '\n';
    mostrar_valor("Diferencia hacia atras", D_atras);
    mostrar_valor("Error absoluto", err_atras);
    cout << '\n';
    mostrar_valor("Diferencia centrada", D_centrada);
    mostrar_valor("Error absoluto", err_centrada);

    linea('=');
    cout << "SEGUNDA DERIVADA\n";
    linea('=');
    mostrar_valor("Aproximacion", D2);
    mostrar_valor("Error absoluto", err_D2);

    linea('=');
    cout << "FIN DEL PROGRAMA\n";
    linea('=');

    return 0;
}