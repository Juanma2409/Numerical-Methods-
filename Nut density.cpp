#include <iostream>
#include <cmath>
using namespace std;

int main() {
    // Datos experimentales
    double D = 17.2;   // mm (distancia entre caras)
    double d = 8.5;    // mm (diámetro interno)
    double h = 6.4;    // mm (altura)
    double m = 12.3;   // g (masa)

    // Incertidumbres instrumentales
    double deltaD = 0.05;  // mm
    double deltad = 0.05;  // mm
    double deltah = 0.05;  // mm
    double deltam = 0.05;  // g

    // Constantes
    double pi = M_PI;
    double sqrt3 = sqrt(3);

    // Volumen (mm³)
    double V_hex = (3 * sqrt3 / 8) * pow(D, 2) * h;
    double V_cil = pi * pow(d, 2) / 4 * h;
    double V = V_hex - V_cil;

    // Convertir a cm³
    double V_cm3 = V / 1000.0;

    // Derivadas parciales
    double dV_dD = (3 * sqrt3 / 4) * D * h;
    double dV_dd = -(pi / 2) * d * h;
    double dV_dh = (3 * sqrt3 / 8) * pow(D, 2) - (pi / 4) * pow(d, 2);

    // Error del volumen
    double deltaV = sqrt(
        pow(dV_dD * deltaD, 2) +
        pow(dV_dd * deltad, 2) +
        pow(dV_dh * deltah, 2)
    );

    // Convertir error a cm³
    double deltaV_cm3 = deltaV / 1000.0;

    // Densidad
    double densidad = m / V_cm3;

    // Derivadas para densidad
    double dDen_dm = 1 / V_cm3;
    double dDen_dV = -m / pow(V_cm3, 2);

    // Error de densidad
    double deltaDen = sqrt(
        pow(dDen_dm * deltam, 2) +
        pow(dDen_dV * deltaV_cm3, 2)
    );

    // Resultados
    cout << "RESULTADOS EXPERIMENTALES\n";
    cout << "=========================================\n";
    cout << "Volumen = (" << V_cm3 << " ± " << deltaV_cm3 << ") cm^3\n";
    cout << "Densidad = (" << densidad << " ± " << deltaDen << ") g/cm^3\n";
    cout << "=========================================\n";
    cout << "Error relativo del volumen = " << (deltaV_cm3 / V_cm3) * 100 << " %\n";
    cout << "Error relativo de la densidad = " << (deltaDen / densidad) * 100 << " %\n";
    cout << "=========================================\n";

    return 0;
}