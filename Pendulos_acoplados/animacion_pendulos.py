#!/usr/bin/env python3
"""
============================================================================
ANIMACIÓN 2D - PÉNDULOS ACOPLADOS POR RESORTE
============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, FancyBboxPatch
import matplotlib.patches as mpatches

# ====================== CONFIGURACIÓN ======================
L = 1.0                    # Longitud del péndulo (¡debe coincidir con tu simulación C++!)
d = 1.0                    # Distancia horizontal entre puntos de suspensión
TRAIL_LENGTH = 40          # Longitud de la estela (más = más largo el rastro)
DOWNSAMPLE = 30             # Usar 1 de cada DOWNSAMPLE frames (1 = todos, 2 = más rápido)
FPS = 30                   # Frames por segundo de la animación

def cargar_datos(archivo='datos_simulacion.dat'):
    try:
        data = np.loadtxt(archivo, comments='#')
        t = data[:, 0]
        theta1 = data[:, 1]
        theta2 = data[:, 3]
        return t, theta1, theta2
    except FileNotFoundError:
        print(f"Error: No se encontró '{archivo}'")
        exit(1)

def draw_spring(ax, x1, y1, x2, y2, num_coils=12, amplitude=0.04, color='#2ca02c', linewidth=2.0):
    """Dibuja un resorte realista con forma de onda sinusoidal"""
    dx = x2 - x1
    dy = y2 - y1
    length = np.sqrt(dx**2 + dy**2)
    if length < 0.01:
        return ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth)[0]
    
    # Vector unitario
    ux, uy = dx / length, dy / length
    # Vector perpendicular
    px, py = -uy, ux
    
    # Puntos a lo largo del resorte
    s = np.linspace(0, length, 200)
    wave = amplitude * np.sin(2 * np.pi * num_coils * s / length)
    
    x = x1 + ux * s + px * wave
    y = y1 + uy * s + py * wave
    
    return ax.plot(x, y, color=color, linewidth=linewidth)[0]

def crear_animacion_avanzada(t, theta1, theta2, guardar=False, 
                             formato='gif', nombre_archivo='animacion_pendulos_calidad.gif'):
    """
    Crea una animación con estelas y resorte realista.
    """
    # Downsample para mejor rendimiento
    idx = np.arange(0, len(t), DOWNSAMPLE)
    t = t[idx]
    theta1 = theta1[idx]
    theta2 = theta2[idx]
    
    # Calcular posiciones
    x1 = L * np.sin(theta1)
    y1 = -L * np.cos(theta1)
    x2 = d + L * np.sin(theta2)
    y2 = -L * np.cos(theta2)
    
    # Figura y ejes
    fig, ax = plt.subplots(figsize=(12, 9), dpi=120)
    ax.set_xlim(-1.8, d + 1.8)
    ax.set_ylim(-1.6, 0.6)
    ax.set_aspect('equal')
    ax.set_facecolor('#0d1b2a')           # Fondo oscuro elegante
    fig.patch.set_facecolor('#0d1b2a')
    
    # Techo
    ax.plot([-0.8, d + 0.8], [0, 0], color='#415a77', linewidth=12, solid_capstyle='round', zorder=1)
    
    # Puntos de suspensión
    ax.add_patch(Circle((0, 0), 0.04, color='#778da9', zorder=5))
    ax.add_patch(Circle((d, 0), 0.04, color='#778da9', zorder=5))
    
    # Estelas (trails)
    trail1, = ax.plot([], [], color='#1f77b4', alpha=0.25, linewidth=1.5)
    trail2, = ax.plot([], [], color='#d62728', alpha=0.25, linewidth=1.5)
    
    # Péndulos (líneas)
    rod1, = ax.plot([], [], color='#1f77b4', linewidth=3.5, solid_capstyle='round', zorder=6)
    rod2, = ax.plot([], [], color='#d62728', linewidth=3.5, solid_capstyle='round', zorder=6)
    
    # Masas
    mass1 = Circle((0, 0), 0.09, color='#1f77b4', ec='white', linewidth=2, zorder=10)
    mass2 = Circle((0, 0), 0.09, color='#d62728', ec='white', linewidth=2, zorder=10)
    ax.add_patch(mass1)
    ax.add_patch(mass2)
    
    # Etiquetas en las masas
    label1 = ax.text(0, 0, 'm₁', color='white', ha='center', va='center', fontsize=9, fontweight='bold', zorder=11)
    label2 = ax.text(0, 0, 'm₂', color='white', ha='center', va='center', fontsize=9, fontweight='bold', zorder=11)
    
    # Resorte (se actualiza dinámicamente)
    spring_line, = ax.plot([], [], color='#2ca02c', linewidth=2.2, zorder=7)
    
    # Título y texto
    ax.set_title('Animación 2D \nPéndulos Acoplados por Resorte — RK4', 
                 color='white', fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel('Posición X (m)', color='white', fontsize=11)
    ax.set_ylabel('Posición Y (m)', color='white', fontsize=11)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('white')
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=13,
                        color='white', bbox=dict(boxstyle='round,pad=0.4', 
                        facecolor='#1b263b', edgecolor='#415a77', alpha=0.9))
    
    # Variables para las estelas
    trail_x1, trail_y1 = [], []
    trail_x2, trail_y2 = [], []
    
    def init():
        rod1.set_data([], [])
        rod2.set_data([], [])
        mass1.center = (0, 0)
        mass2.center = (0, 0)
        spring_line.set_data([], [])
        trail1.set_data([], [])
        trail2.set_data([], [])
        label1.set_position((0, 0))
        label2.set_position((0, 0))
        time_text.set_text('')
        return rod1, rod2, mass1, mass2, spring_line, trail1, trail2, label1, label2, time_text
    
    def animate(i):
        # Actualizar péndulo 1
        rod1.set_data([0, x1[i]], [0, y1[i]])
        mass1.center = (x1[i], y1[i])
        label1.set_position((x1[i], y1[i]))
        
        # Actualizar péndulo 2
        rod2.set_data([d, x2[i]], [0, y2[i]])
        mass2.center = (x2[i], y2[i])
        label2.set_position((x2[i], y2[i]))
        
        # Actualizar resorte (versión mejorada con onda)
        dx = x2[i] - x1[i]
        dy = y2[i] - y1[i]
        length = np.sqrt(dx**2 + dy**2)
        if length > 0.05:
            num_points = 180
            s = np.linspace(0, 1, num_points)
            wave = 0.035 * np.sin(2 * np.pi * 14 * s)   # 14 vueltas
            # Rotar la onda según la dirección del resorte
            angle = np.arctan2(dy, dx)
            x_wave = x1[i] + s * dx - wave * np.sin(angle)
            y_wave = y1[i] + s * dy + wave * np.cos(angle)
            spring_line.set_data(x_wave, y_wave)
        else:
            spring_line.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        
        # Actualizar estelas
        trail_x1.append(x1[i])
        trail_y1.append(y1[i])
        trail_x2.append(x2[i])
        trail_y2.append(y2[i])
        
        if len(trail_x1) > TRAIL_LENGTH:
            trail_x1.pop(0)
            trail_y1.pop(0)
            trail_x2.pop(0)
            trail_y2.pop(0)
        
        trail1.set_data(trail_x1, trail_y1)
        trail2.set_data(trail_x2, trail_y2)
        
        # Tiempo
        time_text.set_text(f'Tiempo = {t[i]:.2f} s')
        
        return rod1, rod2, mass1, mass2, spring_line, trail1, trail2, label1, label2, time_text
    
    # Crear animación
    ani = FuncAnimation(fig, animate, frames=len(t), init_func=init,
                        interval=1000 / FPS, blit=False, repeat=True,
                        cache_frame_data=False)
    
    plt.tight_layout()
    
    if guardar:
        print(f"Guardando animación como {nombre_archivo}...")
        if formato == 'gif':
            ani.save(nombre_archivo, writer='pillow', fps=FPS, dpi=120)
        elif formato == 'mp4':
            ani.save(nombre_archivo, writer='ffmpeg', fps=FPS, dpi=120, 
                     extra_args=['-vcodec', 'libx264'])
        print(f"¡Animación guardada exitosamente: {nombre_archivo}")
    else:
        print("Mostrando animación ...")
        plt.show()
    
    return ani

if __name__ == "__main__":
    print("="*60)
    print("ANIMACIÓN 2D - PÉNDULOS ACOPLADOS")
    print("="*60)
    
    t, theta1, theta2 = cargar_datos()
    print(f"Datos cargados: {len(t)} frames (después de downsample: {len(t)//DOWNSAMPLE})")
    
    # === CONFIGURACIÓN ===
    # Cambia estos valores según lo que quieras:
    GUARDAR = False                    # True = guardar archivo, False = mostrar en pantalla
    FORMATO = 'gif'                    # 'gif' o 'mp4' (mp4 requiere ffmpeg instalado)
    NOMBRE = 'animacion_pendulos.gif'
    
    crear_animacion_avanzada(t, theta1, theta2, 
                             guardar=GUARDAR, 
                             formato=FORMATO, 
                             nombre_archivo=NOMBRE)
