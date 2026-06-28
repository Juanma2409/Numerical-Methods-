#!/usr/bin/env python3
"""
============================================================================
SCRIPT DE VISUALIZACIÓN CIENTÍFICA
Péndulos Acoplados por Resorte - Simulación RK4
============================================================================
Genera gráficas profesionales:
  - Ángulos vs tiempo
  - Espacios de fase 
  - Conservación de energía
  - Comparación de modos 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

# Configuración de estilo profesional
rcParams['figure.figsize'] = (14, 10)
rcParams['font.size'] = 11
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 14
rcParams['legend.fontsize'] = 10
rcParams['lines.linewidth'] = 1.8
rcParams['grid.alpha'] = 0.3

def cargar_datos(archivo='datos_simulacion.dat'):
    """Carga los datos generados por el programa C++"""
    try:
        data = np.loadtxt(archivo, comments='#')
        t = data[:, 0]
        theta1 = data[:, 1]
        omega1 = data[:, 2]
        theta2 = data[:, 3]
        omega2 = data[:, 4]
        E = data[:, 5]
        return t, theta1, omega1, theta2, omega2, E
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo '{archivo}'")
        print("Primero ejecuta el programa C++ para generar los datos.")
        exit(1)

def crear_graficas(t, theta1, omega1, theta2, omega2, E):
    """Genera todas las gráficas científicas"""
    
    # ==================== FIGURA 1: ÁNGULOS Y ENERGÍA ====================
    fig1 = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(3, 2, figure=fig1, hspace=0.35, wspace=0.3)
    
    # Subplot 1: Ángulos vs Tiempo
    ax1 = fig1.add_subplot(gs[0, :])
    ax1.plot(t, theta1, label=r'$\theta_1(t)$', color='#1f77b4')
    ax1.plot(t, theta2, label=r'$\theta_2(t)$', color='#d62728')
    ax1.set_xlabel('Tiempo $t$ (s)')
    ax1.set_ylabel('Ángulo (rad)')
    ax1.set_title('Evolución Temporal de los Ángulos', fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True)
    ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Subplot 2: Energía Total
    ax2 = fig1.add_subplot(gs[1, :])
    ax2.plot(t, E, color='#2ca02c', label='Energía Total $E(t)$')
    E_mean = np.mean(E)
    ax2.axhline(E_mean, color='red', linestyle='--', alpha=0.7, 
                label=f'Energía media = {E_mean:.6f} J')
    ax2.set_xlabel('Tiempo $t$ (s)')
    ax2.set_ylabel('Energía (J)')
    ax2.set_title('Conservación de la Energía (Validación del Método RK4)', fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True)
    
    # Subplot 3 y 4: Espacios de fase
    ax3 = fig1.add_subplot(gs[2, 0])
    ax3.plot(theta1, omega1, color='#1f77b4', alpha=0.85)
    ax3.set_xlabel(r'$\theta_1$ (rad)')
    ax3.set_ylabel(r'$\omega_1$ (rad/s)')
    ax3.set_title('Retrato de Fase - Péndulo 1', fontweight='bold')
    ax3.grid(True)
    ax3.axhline(0, color='gray', linestyle='-', alpha=0.4)
    ax3.axvline(0, color='gray', linestyle='-', alpha=0.4)
    
    ax4 = fig1.add_subplot(gs[2, 1])
    ax4.plot(theta2, omega2, color='#d62728', alpha=0.85)
    ax4.set_xlabel(r'$\theta_2$ (rad)')
    ax4.set_ylabel(r'$\omega_2$ (rad/s)')
    ax4.set_title('Retrato de Fase - Péndulo 2', fontweight='bold')
    ax4.grid(True)
    ax4.axhline(0, color='gray', linestyle='-', alpha=0.4)
    ax4.axvline(0, color='gray', linestyle='-', alpha=0.4)
    
    plt.suptitle('Simulación de Péndulos Acoplados por Resorte - RK4', 
                 fontsize=16, fontweight='bold', y=0.98)
    plt.savefig('grafica_1_angulos_energia_fase.png', dpi=300, bbox_inches='tight')
    print("✓ Guardada: grafica_1_angulos_energia_fase.png")
    
    # ==================== FIGURA 2: COMPARACIÓN Y DETALLES ====================
    fig2, axes = plt.subplots(2, 2, figsize=(14, 9))
    
    # Diferencia de ángulos (oscilación del resorte)
    axes[0,0].plot(t, theta1 - theta2, color='#9467bd')
    axes[0,0].set_xlabel('Tiempo $t$ (s)')
    axes[0,0].set_ylabel(r'$\theta_1 - \theta_2$ (rad)')
    axes[0,0].set_title('Diferencia Angular (Estiramiento del Resorte)', fontweight='bold')
    axes[0,0].grid(True)
    
    # Velocidades angulares
    axes[0,1].plot(t, omega1, label=r'$\omega_1(t)$', color='#1f77b4')
    axes[0,1].plot(t, omega2, label=r'$\omega_2(t)$', color='#d62728')
    axes[0,1].set_xlabel('Tiempo $t$ (s)')
    axes[0,1].set_ylabel('Velocidad angular (rad/s)')
    axes[0,1].set_title('Velocidades Angulares', fontweight='bold')
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    # Zoom en energía 

    
    axes[1,0].plot(t, E - E_mean, color='#2ca02c')
    axes[1,0].set_xlabel('Tiempo $t$ (s)')
    axes[1,0].set_ylabel(r'$E(t) - \langle E \rangle$ (J)')
    axes[1,0].set_title('Deriva de Energía (debe ser ~0)', fontweight='bold')
    axes[1,0].grid(True)
    
    # Histograma de energía (distribución)
    axes[1,1].hist(E, bins=50, color='#ff7f0e', alpha=0.7, edgecolor='black')
    axes[1,1].axvline(E_mean, color='red', linestyle='--', linewidth=2, label='Media')
    axes[1,1].set_xlabel('Energía (J)')
    axes[1,1].set_ylabel('Frecuencia')
    axes[1,1].set_title('Distribución de la Energía Total', fontweight='bold')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.suptitle('Análisis Detallado del Sistema Acoplado', 
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig('grafica_2_analisis_detallado.png', dpi=300, bbox_inches='tight')
    print("✓ Guardada: grafica_2_analisis_detallado.png")
    
    # ==================== FIGURA 3: RETRATOS DE FASE COMBINADOS ====================
    fig3, ax = plt.subplots(figsize=(10, 8))
    ax.plot(theta1, omega1, label=r'Péndulo 1 ($\theta_1, \omega_1$)', color='#1f77b4', alpha=0.8)
    ax.plot(theta2, omega2, label=r'Péndulo 2 ($\theta_2, \omega_2$)', color='#d62728', alpha=0.8)
    ax.set_xlabel('Ángulo (rad)')
    ax.set_ylabel('Velocidad angular (rad/s)')
    ax.set_title('Retratos de Fase Superpuestos\n(Sistema Acoplado)', fontweight='bold', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True)
    ax.axhline(0, color='gray', linestyle='-', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='-', alpha=0.5)
    plt.savefig('grafica_3_retratos_fase_combinados.png', dpi=300, bbox_inches='tight')
    print("✓ Guardada: grafica_3_retratos_fase_combinados.png")
    
    print("\n¡Todas las gráficas han sido generadas exitosamente!")

if __name__ == "__main__":
    print("Cargando datos de simulación...")
    t, theta1, omega1, theta2, omega2, E = cargar_datos()
    print(f"Datos cargados: {len(t)} puntos de simulación")
    print(f"Rango de tiempo: {t[0]:.2f} - {t[-1]:.2f} s")
    print(f"Energía media: {np.mean(E):.8f} J | Desviación estándar: {np.std(E):.2e} J\n")
    
    crear_graficas(t, theta1, omega1, theta2, omega2, E)
