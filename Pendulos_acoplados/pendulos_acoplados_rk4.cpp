/*
 * ============================================================================
 * SIMULADOR CIENTÍFICO: DOS PÉNDULOS SIMPLES ACOPLADOS POR RESORTE
 * Método: Runge-Kutta de 4º orden (RK4)
 * Curso: Métodos Numéricos
 * ============================================================================
 *
 * DESCRIPCIÓN FÍSICA:
 * Dos péndulos de longitud b y masa m, suspendidos de puntos separados
 * una distancia b (en reposo). Conectados por un resorte de constante k.
 * El sistema tiene 2 grados de libertad: θ₁(t) y θ₂(t).
 *
 * ENERGÍAS (del problema):
 *   T = (1/2) m b² (θ̇₁² + θ̇₂²)
 *   U = m g b (1 - cos θ₁) + m g b (1 - cos θ₂) 
 *       + (1/2) k b² (sin θ₁ - sin θ₂)²
 *
 * ECUACIONES DE MOVIMIENTO (derivadas del Lagrangiano L = T - U):
 *   θ̈₁ = - (g/b) sin θ₁ - (k/m) (sin θ₁ - sin θ₂) cos θ₁
 *   θ̈₂ = - (g/b) sin θ₂ + (k/m) (sin θ₁ - sin θ₂) cos θ₂
 *
 * Estas ecuaciones son NO LINEALES. RK4 las resuelve numéricamente
 * sin necesidad de aproximaciones de pequeño ángulo.
 *
 * ============================================================================
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include <ios>  // para numeric_limits

#ifdef _WIN32
    #include <windows.h>   // Para configurar UTF-8 en consola de Windows
#endif

using namespace std;

// Estructura para parámetros físicos
struct Parametros {
    double m;   // masa (kg)
    double b;   // longitud del péndulo (m)
    double k;   // constante del resorte (N/m)
    double g;   // aceleración gravedad (m/s²)
};

// ====================== FUNCIÓN DE DERIVADAS (f(t, y)) ======================
/*
 * y = [θ₁, ω₁, θ₂, ω₂]
 * dydt = [ω₁, θ̈₁, ω₂, θ̈₂]
 */
vector<double> derivadas(const vector<double>& y, const Parametros& p) {
    vector<double> dydt(4, 0.0);
    
    double theta1 = y[0];
    double omega1 = y[1];
    double theta2 = y[2];
    double omega2 = y[3];
    
    double sin1 = sin(theta1);
    double cos1 = cos(theta1);
    double sin2 = sin(theta2);
    double cos2 = cos(theta2);
    double diff_sin = sin1 - sin2;
    
    // Ecuaciones de movimiento derivadas del Lagrangiano
    double accel1 = -(p.g / p.b) * sin1 - (p.k / p.m) * diff_sin * cos1;
    double accel2 = -(p.g / p.b) * sin2 + (p.k / p.m) * diff_sin * cos2;
    
    dydt[0] = omega1;   // dθ₁/dt = ω₁
    dydt[1] = accel1;   // dω₁/dt = θ̈₁
    dydt[2] = omega2;   // dθ₂/dt = ω₂
    dydt[3] = accel2;   // dω₂/dt = θ̈₂
    
    return dydt;
}

// ====================== PASO DE RUNGE-KUTTA 4 ======================
/*
 * Implementación clásica de RK4 (método explícito de 4 etapas).
 * Es de orden 4: error local O(h^5), global O(h^4).
 * Muy estable y preciso para problemas de mecánica.
 */
vector<double> rk4_step(const vector<double>& y, double h, const Parametros& p) {
    auto k1 = derivadas(y, p);
    
    vector<double> y2(4);
    for(int i = 0; i < 4; i++) y2[i] = y[i] + (h/2.0) * k1[i];
    auto k2 = derivadas(y2, p);
    
    vector<double> y3(4);
    for(int i = 0; i < 4; i++) y3[i] = y[i] + (h/2.0) * k2[i];
    auto k3 = derivadas(y3, p);
    
    vector<double> y4(4);
    for(int i = 0; i < 4; i++) y4[i] = y[i] + h * k3[i];
    auto k4 = derivadas(y4, p);
    
    vector<double> y_new(4);
    for(int i = 0; i < 4; i++) {
        y_new[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    return y_new;
}

// ====================== CÁLCULO DE ENERGÍA TOTAL ======================
double calcular_energia(const vector<double>& y, const Parametros& p) {
    double T = 0.5 * p.m * p.b * p.b * (y[1]*y[1] + y[3]*y[3]);           // Cinética
    double U_grav = p.m * p.g * p.b * ( (1 - cos(y[0])) + (1 - cos(y[2])) );
    double U_resorte = 0.5 * p.k * p.b * p.b * pow(sin(y[0]) - sin(y[2]), 2);
    return T + U_grav + U_resorte;
}

// ====================== FRECUENCIAS ANALÍTICAS (PEQUEÑO ÁNGULO) ======================
/*
 * Para comparación: en aproximación de pequeño ángulo el sistema es lineal.
 * Frecuencias de los modos normales:
 *   ω_sym  = sqrt(g/b)
 *   ω_asym = sqrt( g/b + 2k/m )
 */
void imprimir_frecuencias_analiticas(const Parametros& p) {
    double omega_sym = sqrt(p.g / p.b);
    double omega_asym = sqrt(p.g / p.b + 2.0 * p.k / p.m);
    
    cout << "\n=== COMPARACIÓN CON SOLUCIÓN ANALÍTICA (pequeño ángulo) ===\n";
    cout << "Frecuencia modo simétrico (ω_sym):   " << omega_sym << " rad/s\n";
    cout << "Frecuencia modo antisimétrico (ω_asym): " << omega_asym << " rad/s\n";
    cout << "Período modo simétrico:  " << 2*M_PI/omega_sym << " s\n";
    cout << "Período modo antisimétrico: " << 2*M_PI/omega_asym << " s\n";
}

// ====================== PROGRAMA PRINCIPAL ======================
int main() {
    // ====================== CONFIGURACIÓN UTF-8 EN WINDOWS ======================
    #ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
        SetConsoleCP(CP_UTF8);
    #endif

    cout << fixed << setprecision(8);
    cout << "\n";
    cout << "╔════════════════════════════════════════════════════════════════════════════╗\n";
    cout << "║   SIMULADOR DE PÉNDULOS ACOPLADOS - RUNGE-KUTTA 4º ORDEN (RK4)            ║\n";
    cout << "║   Proyecto Final - Métodos Numéricos                                       ║\n";
    cout << "╚════════════════════════════════════════════════════════════════════════════╝\n\n";
    
    Parametros p;
    
    // ====================== ENTRADA DE PARÁMETROS ======================
    cout << ">>> INGRESO DE PARÁMETROS FÍSICOS\n";
    cout << "Masa de cada péndulo m (kg)           : "; cin >> p.m;
    cout << "Longitud de cada péndulo (m)    : "; cin >> p.b;
    cout << "Constante del resorte k (N/m)     : "; cin >> p.k;
    cout << "Gravedad g (m/s², usualmente 9.81): "; cin >> p.g;
    
    // ====================== CONDICIONES INICIALES ======================
    cout << "\n>>> CONDICIONES INICIALES\n";
    cout << "Seleccione el tipo de simulación:\n";
    cout << "  1. Modo simétrico     (θ1 = θ0,  θ2 = θ0)   - Ambos péndulos en fase\n";
    cout << "  2. Modo antisimétrico (θ1 = θ0,  θ2 = -θ0)  - En oposición de fase\n";
    cout << "  3. Personalizado      (ingresar θ1 y θ2 por separado)\n";
    cout << "Opción: ";
    
    int modo;
    cin >> modo;
    
    double theta1_0 = 0.0, theta2_0 = 0.0;
    double theta0;
    
    if (modo == 1) {
        cout << "Ángulo inicial común θ0 (rad, ej. 0.3): "; cin >> theta0;
        theta1_0 = theta0;
        theta2_0 = theta0;
        cout << "→ Modo SIMÉTRICO seleccionado.\n";
    } 
    else if (modo == 2) {
        cout << "Ángulo inicial |θ0| (rad, ej. 0.3): "; cin >> theta0;
        theta1_0 = theta0;
        theta2_0 = -theta0;
        cout << "→ Modo ANTISIMÉTRICO seleccionado.\n";
    } 
    else {
        cout << "Ángulo inicial θ1 (rad): "; cin >> theta1_0;
        cout << "Ángulo inicial θ2 (rad): "; cin >> theta2_0;
        cout << "→ Condiciones PERSONALIZADAS.\n";
    }
    
    // ====================== PARÁMETROS DE INTEGRACIÓN ======================
    cout << "\n>>> PARÁMETROS DE INTEGRACIÓN NUMÉRICA\n";
    double t_max, h;
    cout << "Tiempo total de simulación t_max (s, ej. 30): "; cin >> t_max;
    cout << "Paso de tiempo h (s, recomendado 0.001 - 0.01): "; cin >> h;
    
    // ====================== SIMULACIÓN ======================
    vector<double> y = {theta1_0, 0.0, theta2_0, 0.0};  // Estado inicial: [θ1, ω1, θ2, ω2]
    
    ofstream datos("datos_simulacion.dat");
    datos << "# t(s)   theta1(rad)   omega1(rad/s)   theta2(rad)   omega2(rad/s)   E_total(J)\n";
    
    double E_inicial = calcular_energia(y, p);
    double E_max = E_inicial;
    double E_min = E_inicial;
    
    cout << "\n>>> INICIANDO SIMULACIÓN RK4...\n";
    cout << "Energía inicial del sistema: " << E_inicial << " J\n\n";
    
    int paso = 0;
    int pasos_totales = static_cast<int>(t_max / h);
    
    for (double t = 0.0; t < t_max + 1e-9; t += h) {
        // Guardar datos cada paso (o cada N pasos si se quiere reducir tamaño)
        double E = calcular_energia(y, p);
        if (E > E_max) E_max = E;
        if (E < E_min) E_min = E;
        
        datos << t << "  " 
              << y[0] << "  " << y[1] << "  " 
              << y[2] << "  " << y[3] << "  " 
              << E << "\n";
        
        // Mostrar progreso cada ~5% o cada segundo
        if (paso % max(1, pasos_totales / 20) == 0) {
            cout << "  t = " << setw(8) << t << " s   |   θ1 = " << setw(10) << y[0] 
                 << " rad   |   θ2 = " << setw(10) << y[2] << " rad   |   E = " << E << "\n";
        }
        
        y = rk4_step(y, h, p);
        paso++;
    }
    
    datos.close();
    
    double error_energia = fabs(E_max - E_min) / E_inicial * 100.0;
    
    cout << "\n╔════════════════════════════════════════════════════════════════════════════╗\n";
    cout << "║                           RESULTADOS DE LA SIMULACIÓN                       ║\n";
    cout << "╚════════════════════════════════════════════════════════════════════════════╝\n";
    cout << "Paso de tiempo utilizado (h): " << h << " s\n";
    cout << "Número de pasos RK4: " << paso << "\n";
    cout << "Energía inicial:     " << E_inicial << " J\n";
    cout << "Energía máxima:      " << E_max << " J\n";
    cout << "Energía mínima:      " << E_min << " J\n";
    cout << "Deriva energética:   " << error_energia << " %  (debe ser muy pequeña)\n";
    
    imprimir_frecuencias_analiticas(p);
    
    cout << "\nArchivos generados:\n";
    cout << "  • datos_simulacion.dat   → Datos completos para graficar\n";
    cout << "\nPara generar las gráficas ejecuta el script de Python:\n";
    cout << "  python3 graficar_pendulos.py\n\n";
    
    cout << "¡Simulación completada exitosamente!\n";
    
    // ====================== PAUSA PARA WINDOWS (evita que se cierre la consola) ======================
    cout << "\n\n════════════════════════════════════════════════════════════════════════════\n";
    cout << "Presiona ENTER para cerrar esta ventana...\n";
    cin.ignore(numeric_limits<streamsize>::max(), '\n');  // limpia buffer
    cin.get();                                            // espera ENTER
    
    return 0;
}
