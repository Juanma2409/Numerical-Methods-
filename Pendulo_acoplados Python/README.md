# Péndulos Acoplados — Análisis de Caos en Tiempo Real

Simulación interactiva de dos péndulos acoplados por un resorte, con visualización en tiempo real del retrato de fase y la evolución energética del sistema. Implementado en Python con Pygame y un integrador RK4 escrito desde cero.

---

## El problema físico

Dos péndulos planos de masa `m` y longitud `L` cuelgan de pivotes fijos separados una distancia de 1.5 m. Entre las masas hay un resorte horizontal de constante `k` que transmite fuerza cuando los ángulos difieren. Las ecuaciones de movimiento son:

```
θ̈₁ = -(g/L)·sin(θ₁) - (k/m)·(sin θ₁ − sin θ₂)·cos θ₁
θ̈₂ = -(g/L)·sin(θ₂) + (k/m)·(sin θ₁ − sin θ₂)·cos θ₂
```

El sistema tiene 4 grados de libertad: `[θ₁, ω₁, θ₂, ω₂]`. Para amplitudes pequeñas el comportamiento es periódico y predecible. Para amplitudes grandes el acoplamiento no lineal puede llevar al **caos determinista**: dos condiciones iniciales arbitrariamente cercanas divergen exponencialmente con el tiempo.

---

## Integrador numérico — RK4

El integrador es una implementación explícita del método de Runge-Kutta de 4° orden, sin usar ninguna librería de ODEs:

```python
def rk4(y, k, h):
    k1 = derivs(y, k)
    k2 = derivs(y + h/2 * k1, k)
    k3 = derivs(y + h/2 * k2, k)
    k4 = derivs(y + h * k3, k)
    return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
```

El paso de integración es `h = 0.004 s` y se ejecutan 10 pasos por frame, lo que da una velocidad de simulación de `0.04 s / frame` a 60 FPS.

---

## Paneles de visualización

La ventana está dividida en 3 paneles:

```
┌─────────────────────┬──────────────────────┐
│   ANIMACIÓN 2D      │  RETRATOS DE FASE    │
│   péndulos+resorte  │  θ₁/ω₁  y  θ₂/ω₂   │
├─────────────────────┴──────────────────────┤
│         ENERGÍA  Ep(t)  vs  E_total(t)     │
└────────────────────────────────────────────┘
```

**Panel superior izquierdo — Animación 2D**

Muestra los dos péndulos (azul y rojo) con sus estelas de trayectoria y el resorte que los conecta. En la esquina inferior se muestran el tiempo simulado, la constante del resorte actual, el error de energía relativo `ΔE/E₀` (indicador de precisión numérica) y el estado de pausa.

**Panel superior derecho — Retratos de fase**

Muestra `θ vs ω` para cada péndulo. El ángulo se normaliza siempre al rango `[-π, π]`, de modo que el retrato permanece confinado incluso cuando el péndulo da vueltas completas. El punto blanco indica la posición actual del sistema en el espacio de fases.

- Sistema **regular**: el retrato forma curvas cerradas y ordenadas.
- Sistema **caótico**: el retrato llena el espacio de forma errática sin patrón aparente.

**Panel inferior — Energía Ep(t) vs E_total(t)**

Muestra en el mismo gráfico:
- **Ep(t)** en púrpura: energía potencial (gravitacional + elástica), que oscila continuamente.
- **E_total(t)** en cian: energía mecánica total, que debería mantenerse constante si el integrador no acumula error. Una línea prácticamente plana confirma que RK4 conserva bien la energía; si se inclina notablemente, indica acumulación de error numérico.

Los valores actuales de `Ep`, `Ek` y `E_total` se muestran en tiempo real en la esquina.

---

## Requisitos

```
python >= 3.9
pygame >= 2.0
numpy  >= 1.21
```

Instalar dependencias:

```bash
pip install pygame numpy
```

---

## Cómo ejecutar

```bash
python caos_pendulos_pygame1.py
```

La ventana se abre a 1920×1080. Se puede redimensionar libremente con el ratón o usar pantalla completa con `F`.

---

## Controles

| Tecla | Acción |
|---|---|
| `ESPACIO` | Pausar / reanudar la simulación |
| `R` | Reiniciar con las mismas condiciones iniciales |
| `1` | Condición inicial: modo simétrico — regular (`θ₁=θ₂=0.30 rad`) |
| `2` | Condición inicial: modo antisimétrico — regular (`θ₁=0.50, θ₂=−0.50 rad`) |
| `3` | Condición inicial: caótico leve (`θ₁=0.90, θ₂=−0.40 rad`) |
| `4` | Condición inicial: caótico fuerte (`θ₁=1.40, θ₂=−1.00 rad`) |
| `↑` | Aumentar constante del resorte `k` en +0.2 N/m |
| `↓` | Disminuir constante del resorte `k` en −0.2 N/m (mínimo 0) |
| `F` | Alternar pantalla completa |
| `S` | Guardar captura de pantalla PNG |
| `Q` / `ESC` | Salir |

---

## Parámetros internos

Los siguientes valores se pueden modificar directamente al inicio del archivo para explorar comportamientos distintos:

| Variable | Valor por defecto | Descripción |
|---|---|---|
| `M` | `1.0 kg` | Masa de cada péndulo |
| `B` | `1.0 m` | Longitud de cada péndulo |
| `K` | `1.5 N/m` | Constante del resorte inicial |
| `G` | `9.81 m/s²` | Gravedad |
| `H_STEP` | `0.004 s` | Paso de integración RK4 |
| `STEPS_FRAME` | `10` | Pasos RK4 por frame (controla velocidad de simulación) |
| `TRAIL_LEN` | `220` | Longitud de la estela de las masas |
| `PHASE_LEN` | `6000` | Puntos almacenados en el retrato de fase |

---

## Qué observar

**Transición al caos con `k`:** empezar con preset `1` (regular) e ir subiendo `k` con `↑`. Al aumentar el acoplamiento, el sistema puede volverse caótico aunque las condiciones iniciales sean pequeñas.

**Sensibilidad a condiciones iniciales:** cargar el preset `4` y comparar visualmente cómo el retrato de fase nunca repite el mismo trazo, a diferencia del preset `1` donde las curvas se superponen perfectamente ciclo a ciclo.

**Error numérico:** observar el panel de energía en sesiones largas. Un buen integradormantendrá la línea `E_total` prácticamente horizontal. Si se baja `H_STEP` a `0.001` el error disminuye; si se sube a `0.02` se vuelve visible la deriva.

---

## Estructura del código

```
caos_pendulos_pygame1.py
│
├── derivs()              # ecuaciones de movimiento (lado derecho del sistema)
├── rk4()                 # integrador RK4 de 4 etapas (implementación propia)
├── wrap_angle()          # normaliza θ a [-π, π] para el retrato de fase
├── energia_potencial()   # Ep gravitacional + elástica
├── energia_cinetica()    # Ek rotacional
├── energia()             # E_total = Ep + Ek
│
├── class Sim             # estado completo de la simulación y buffers de historial
│   ├── __init__()
│   ├── reset()
│   └── step()            # avanza un paso RK4 y actualiza todos los buffers
│
├── draw_panel()          # dibuja fondo y título de un panel
├── draw_spring()         # dibuja el resorte como sinusoide a lo largo del segmento
├── draw_plot()           # gráfica de línea genérica con rango configurable
│
└── main()                # bucle principal: eventos → física → render
```
