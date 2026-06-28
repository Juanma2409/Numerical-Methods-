# SIMULADOR DE PÉNDULOS ACOPLADOS POR RESORTE
## Proyecto Final - Métodos Numéricos | Runge-Kutta 4º Orden

### Descripción del Problema
Dos péndulos simples de longitud **b** y masa **m**, conectados por un resorte de constante **k**.  
El sistema es **no lineal** y se resuelve numéricamente con **RK4**.

### Estructura del Proyecto

| Archivo                        | Descripción                                      |
|--------------------------------|--------------------------------------------------|
| `pendulos_acoplados_rk4.cpp`   | Programa principal en C++ (RK4 + entrada de usuario) |
| `graficar_pendulos.py`         | Script de Python para generar gráficas científicas |
| `datos_simulacion.dat`         | Datos de salida generados por la simulación (después de ejecutar) |
| `README_SIMULADOR.md`          | Este archivo de instrucciones                    |

### Cómo Compilar y Ejecutar

**1. Compilar el programa C++**
```bash
g++ -o pendulos_rk4 pendulos_acoplados_rk4.cpp -lm
```

**2. Ejecutar la simulación**
```bash
./pendulos_rk4
```

El programa te pedirá:
- Parámetros físicos (m, b, k, g)
- Tipo de simulación (Modo Simétrico, Antisimétrico o Personalizado)
- Tiempo total y paso de integración `h`

**3. Generar las gráficas**
```bash
python3 graficar_pendulos.py
```

Se generarán automáticamente 3 imágenes PNG de alta calidad:
- `grafica_1_angulos_energia_fase.png`
- `grafica_2_analisis_detallado.png`
- `grafica_3_retratos_fase_combinados.png`

### Salida del Programa C++

El programa imprime en pantalla:
- Progreso de la simulación
- Energía inicial, máxima, mínima y deriva porcentual (¡importante para validar precisión!)
- Frecuencias analíticas de los modos normales (comparación con aproximación lineal)

### Archivos de Datos Generados

- `datos_simulacion.dat`: contiene 6 columnas:
  1. `t` (tiempo)
  2. `theta1`
  3. `omega1`
  4. `theta2`
  5. `omega2`
  6. `E_total` (energía total del sistema)

### Aspectos Didácticos Incluidos

1. **Derivación completa** de las ecuaciones de movimiento a partir del Lagrangiano (comentada en el código).
2. **Implementación explícita de RK4** paso a paso (sin librerías externas).
3. **Validación numérica** mediante conservación de energía.
4. **Comparación** con solución analítica de pequeño ángulo (frecuencias de modos normales).
5. **Múltiples modos** de excitación inicial (simétrico / antisimétrico).
6. **Gráficas científicas** listas para el informe (ángulos, fase, energía, deriva).
7. **Animación** de los resultados.

