#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pendulos Acoplados - Analisis de Caos en Tiempo Real
Pygame + RK4 - Layout 3 paneles:
  [ANIMACION 2D]      [RETRATOS DE FASE]
  [   ENERGIA Ep(t) vs E_total(t)       ]
"""

import sys
import math
import numpy as np
from collections import deque
import pygame

# ===================================================
#  CONFIGURACION
# ===================================================

W, H  = 1920, 1080
FPS   = 60

M  = 1.0
B  = 1.0
K  = 1.5
G  = 9.81

H_STEP      = 0.004
STEPS_FRAME = 10

TRAIL_LEN  = 220
PHASE_LEN  = 6000
XY_LEN     = 5000

PRESETS = {
    '1': ('Modo simetrico (regular)',      np.array([ 0.30, 0.0,  0.30, 0.0])),
    '2': ('Modo antisimetrico (regular)',  np.array([ 0.50, 0.0, -0.50, 0.0])),
    '3': ('Caotico leve',                  np.array([ 0.90, 0.0, -0.40, 0.2])),
    '4': ('Caotico fuerte',                np.array([ 1.40, 0.0, -1.00, 0.5])),
}
PRESET_DEFAULT = '3'

# ===================================================
#  PALETA DE COLORES
# ===================================================

BG        = (13,  27,  42)
PANEL_BG  = (27,  38,  59)
GRID_C    = (65,  90, 119)
WHITE     = (224, 230, 240)
MUTED     = (119, 141, 169)
BLUE      = ( 31, 119, 180)
RED       = (214,  39,  40)
GREEN     = ( 44, 160,  44)
ORANGE    = (255, 127,  14)
PURPLE    = (148, 103, 189)
PINK      = (227,  87, 194)
YELLOW    = (255, 215,   0)
CEILING   = ( 65,  90, 119)
PIVOT_C   = (119, 141, 169)

# ===================================================
#  FISICA - RK4
# ===================================================

def derivs(y, k):
    th1, w1, th2, w2 = y
    s1 = math.sin(th1); c1 = math.cos(th1)
    s2 = math.sin(th2); c2 = math.cos(th2)
    ds = s1 - s2
    a1 = -(G / B) * s1 - (k / M) * ds * c1
    a2 = -(G / B) * s2 + (k / M) * ds * c2
    return np.array([w1, a1, w2, a2])

def rk4(y, k, h):
    k1 = derivs(y, k)
    k2 = derivs(y + h/2 * k1, k)
    k3 = derivs(y + h/2 * k2, k)
    k4 = derivs(y + h * k3,   k)
    return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

def wrap_angle(th):
    """Normaliza theta al rango [-pi, pi] para el retrato de fase."""
    return ((th + math.pi) % (2 * math.pi)) - math.pi

def energia_potencial(y, k):
    th1, _, th2, _ = y
    Ug = M * G * B * ((1 - math.cos(th1)) + (1 - math.cos(th2)))
    Us = 0.5 * k * B**2 * (math.sin(th1) - math.sin(th2))**2
    return Ug + Us

def energia_cinetica(y):
    _, w1, _, w2 = y
    return 0.5 * M * B**2 * (w1**2 + w2**2)

def energia(y, k):
    return energia_cinetica(y) + energia_potencial(y, k)

# ===================================================
#  ESTADO DE SIMULACION
# ===================================================

class Sim:
    def __init__(self, y0, k):
        self.k   = k
        self.y0  = y0.copy()
        self.Y   = y0.copy()
        self.t   = 0.0
        self.E0  = energia(y0, k)
        self.paused = False

        self.trail1  = deque(maxlen=TRAIL_LEN)
        self.trail2  = deque(maxlen=TRAIL_LEN)
        self.ph1x    = deque(maxlen=PHASE_LEN)
        self.ph1y    = deque(maxlen=PHASE_LEN)
        self.ph2x    = deque(maxlen=PHASE_LEN)
        self.ph2y    = deque(maxlen=PHASE_LEN)

        # Energia vs tiempo
        self.en_t   = deque(maxlen=XY_LEN)
        self.en_pot = deque(maxlen=XY_LEN)
        self.en_tot = deque(maxlen=XY_LEN)

    def reset(self, y0=None, k=None):
        if k  is not None: self.k  = k
        if y0 is not None: self.y0 = y0.copy()
        self.__init__(self.y0, self.k)

    def step(self):
        self.Y = rk4(self.Y, self.k, H_STEP)
        self.t += H_STEP

        th1, w1, th2, w2 = self.Y

        # Cartesianas: pivote en (+-0.75, 0), masa cuelga hacia abajo (y negativo = abajo)
        x1 = -0.75 + B * math.sin(th1)
        y1 = -B * math.cos(th1)
        x2 =  0.75 + B * math.sin(th2)
        y2 = -B * math.cos(th2)

        self.trail1.append((x1, y1))
        self.trail2.append((x2, y2))
        # Retrato de fase con angulo normalizado a [-pi, pi]
        self.ph1x.append(wrap_angle(th1)); self.ph1y.append(w1)
        self.ph2x.append(wrap_angle(th2)); self.ph2y.append(w2)

        # Energia
        self.en_t.append(self.t)
        self.en_pot.append(energia_potencial(self.Y, self.k))
        self.en_tot.append(energia(self.Y, self.k))



# ===================================================
#  DIBUJO - HELPERS
# ===================================================

def draw_panel(surf, rect, title, font_sm):
    x, y, pw, ph = rect
    pygame.draw.rect(surf, PANEL_BG, rect)
    pygame.draw.rect(surf, GRID_C, rect, 1)
    t = font_sm.render(title, True, WHITE)
    surf.blit(t, (x + 8, y + 6))
    return pygame.Rect(x+1, y+28, pw-2, ph-29)


def world_to_screen(wx, wy, cx, cy, scale):
    # wx, wy en metros. cy es la linea del pivote (techo).
    # wy positivo = hacia abajo en pantalla
    return int(cx + wx * scale), int(cy + wy * scale)


def draw_spring(surf, p1, p2, color, n_coils=13, amp=5):
    x1, y1 = p1; x2, y2 = p2
    dx, dy = x2-x1, y2-y1
    ln = math.hypot(dx, dy)
    if ln < 2: return
    ux, uy = dx/ln, dy/ln
    px, py = -uy, ux
    pts = []
    for i in range(181):
        s = i / 180
        wave = amp * math.sin(2 * math.pi * n_coils * s)
        pts.append((x1 + ux*s*ln + px*wave,
                    y1 + uy*s*ln + py*wave))
    pygame.draw.lines(surf, color, False, [(int(x), int(y)) for x, y in pts], 2)


def draw_plot(surf, area, xs, ys, color, y_zero_line=False,
              x_range=None, y_range=None, lw=1):
    if len(xs) < 2: return
    ax, ay, aw, ah = area.x, area.y, area.width, area.height
    xmin = x_range[0] if x_range else min(xs)
    xmax = x_range[1] if x_range else max(xs)
    if xmax == xmin: xmax = xmin + 1
    ys_arr = list(ys)
    ymin = y_range[0] if y_range else min(ys_arr)
    ymax = y_range[1] if y_range else max(ys_arr)
    if ymax == ymin: ymin -= 0.5; ymax += 0.5

    def to_px(x, y):
        px = ax + int((x - xmin) / (xmax - xmin) * aw)
        py = ay + ah - int((y - ymin) / (ymax - ymin) * ah)
        return px, py

    if y_zero_line:
        y0s = ay + ah - int((0 - ymin) / (ymax - ymin) * ah)
        if ay <= y0s <= ay + ah:
            pygame.draw.line(surf, GRID_C, (ax, y0s), (ax+aw, y0s), 1)

    xs_list = list(xs)
    pts = [to_px(xs_list[i], ys_arr[i]) for i in range(len(xs_list))]
    pts = [(max(ax, min(ax+aw, p[0])), max(ay, min(ay+ah, p[1]))) for p in pts]
    if len(pts) >= 2:
        pygame.draw.lines(surf, color, False, pts, lw)


def draw_axes(surf, area, color=GRID_C):
    pygame.draw.rect(surf, (10, 20, 35), area)
    pygame.draw.rect(surf, color, area, 1)


def text_at(surf, font, txt, pos, color=WHITE):
    surf.blit(font.render(txt, True, color), pos)




# ===================================================
#  MAIN
# ===================================================

def main():
    pygame.init()
    screen = pygame.display.set_mode((W, H), pygame.RESIZABLE)
    pygame.display.set_caption("Pendulos Acoplados - Analisis de Caos")
    clock = pygame.time.Clock()
    fullscreen = False

    font    = pygame.font.SysFont('monospace', 18, bold=False)
    font_sm = pygame.font.SysFont('monospace', 16, bold=True)
    font_lg = pygame.font.SysFont('monospace', 20, bold=True)

    y0 = PRESETS[PRESET_DEFAULT][1]
    sim = Sim(y0, K)
    preset_name = PRESETS[PRESET_DEFAULT][0]

    screenshot_n = [0]
    running = True

    while running:
        sw, sh = screen.get_size()
        HW, HH = sw // 2, sh // 2
        panels = {
            'anim':  (0,   0,   HW, HH),
            'phase': (HW,  0,   HW, HH),
            'xy':    (0,   HH,  sw, HH),   # energia ocupa toda la fila inferior
        }

        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                if ev.key in (pygame.K_q, pygame.K_ESCAPE):
                    running = False
                elif ev.key == pygame.K_SPACE:
                    sim.paused = not sim.paused
                elif ev.key == pygame.K_r:
                    sim.reset()
                elif ev.key == pygame.K_f:
                    fullscreen = not fullscreen
                    if fullscreen:
                        screen = pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
                    else:
                        screen = pygame.display.set_mode((W, H), pygame.RESIZABLE)
                elif ev.key == pygame.K_s:
                    screenshot_n[0] += 1
                    fname = f'caos_pendulos_{screenshot_n[0]:02d}.png'
                    pygame.image.save(screen, fname)
                    print(f'Screenshot: {fname}')
                elif ev.key == pygame.K_UP:
                    sim.reset(k=round(sim.k + 0.2, 1))
                elif ev.key == pygame.K_DOWN:
                    sim.reset(k=round(max(0.0, sim.k - 0.2), 1))
                elif ev.unicode in PRESETS:
                    key = ev.unicode
                    preset_name = PRESETS[key][0]
                    sim.reset(y0=PRESETS[key][1])

        if not sim.paused:
            for _ in range(STEPS_FRAME):
                sim.step()

        screen.fill(BG)

        # ===================================================
        #  PANEL A - Animacion 2D
        # ===================================================
        inner_a = draw_panel(screen, panels['anim'], '  ANIMACION 2D', font_sm)
        ax, ay, aw, ah = inner_a

        # Pivote en la parte SUPERIOR, pendulos CUELGAN hacia abajo
        cx = ax + aw // 2
        cy = ay + int(ah * 0.18)    # linea del techo (pivot)
        scale = min(aw, ah) * 0.40

        # Posiciones pantalla de pivotes
        p1x, p1y = world_to_screen(-0.75, 0, cx, cy, scale)
        p2x, p2y = world_to_screen( 0.75, 0, cx, cy, scale)

        th1, w1, th2, w2 = sim.Y
        # masa cuelga: y_mundo = B*cos(th) positivo -> hacia abajo en pantalla
        b1x, b1y = world_to_screen(-0.75 + B*math.sin(th1),  B*math.cos(th1), cx, cy, scale)
        b2x, b2y = world_to_screen( 0.75 + B*math.sin(th2),  B*math.cos(th2), cx, cy, scale)

        # Soporte (techo)
        pygame.draw.line(screen, CEILING, (p1x-50, p1y), (p2x+50, p1y), 10)
        for px_s in range(p1x-50, p2x+55, 16):
            pygame.draw.line(screen, GRID_C, (px_s, p1y), (px_s-8, p1y-9), 2)

        # Estelas masa 1
        if len(sim.trail1) > 1:
            for i in range(1, len(sim.trail1)):
                alpha_i = int(200 * i / len(sim.trail1))
                wx1, wy1 = sim.trail1[i-1]
                wx2, wy2 = sim.trail1[i]
                # trail almacena y negativo; para pantalla necesitamos -wy -> cuelga
                sx1, sy1 = world_to_screen(wx1, -wy1, cx, cy, scale)
                sx2, sy2 = world_to_screen(wx2, -wy2, cx, cy, scale)
                c = (max(0, BLUE[0]-80+alpha_i//3),
                     max(0, BLUE[1]-80+alpha_i//3),
                     min(255, BLUE[2]+alpha_i//3))
                pygame.draw.line(screen, c, (sx1,sy1), (sx2,sy2), 1)

        # Estelas masa 2
        if len(sim.trail2) > 1:
            for i in range(1, len(sim.trail2)):
                wx1, wy1 = sim.trail2[i-1]
                wx2, wy2 = sim.trail2[i]
                sx1, sy1 = world_to_screen(wx1, -wy1, cx, cy, scale)
                sx2, sy2 = world_to_screen(wx2, -wy2, cx, cy, scale)
                alpha_i = int(200 * i / len(sim.trail2))
                c = (min(255, RED[0]+alpha_i//4), max(0, RED[1]), max(0, RED[2]))
                pygame.draw.line(screen, c, (sx1,sy1), (sx2,sy2), 1)

        # Resorte entre masas
        draw_spring(screen, (b1x, b1y), (b2x, b2y), GREEN)

        # Varillas (del pivote a la masa)
        pygame.draw.line(screen, BLUE, (p1x, p1y), (b1x, b1y), 3)
        pygame.draw.line(screen, RED,  (p2x, p2y), (b2x, b2y), 3)

        # Pivotes
        pygame.draw.circle(screen, PIVOT_C, (p1x, p1y), 7)
        pygame.draw.circle(screen, PIVOT_C, (p2x, p2y), 7)

        # Masas
        pygame.draw.circle(screen, BLUE,  (b1x, b1y), 15)
        pygame.draw.circle(screen, WHITE, (b1x, b1y), 15, 2)
        pygame.draw.circle(screen, RED,   (b2x, b2y), 15)
        pygame.draw.circle(screen, WHITE, (b2x, b2y), 15, 2)

        lbl = font_sm.render('m1', True, WHITE)
        screen.blit(lbl, (b1x - lbl.get_width()//2, b1y - lbl.get_height()//2))
        lbl = font_sm.render('m2', True, WHITE)
        screen.blit(lbl, (b2x - lbl.get_width()//2, b2y - lbl.get_height()//2))

        # Info panel A
        E  = energia(sim.Y, sim.k)
        dE = abs(E - sim.E0) / sim.E0 * 100 if sim.E0 else 0
        y_info = ay + ah - 96
        text_at(screen, font, f't = {sim.t:8.2f} s',          (ax+6, y_info),    WHITE)
        text_at(screen, font, f'k = {sim.k:.1f} N/m  [UP/DN]',(ax+6, y_info+24), MUTED)
        text_at(screen, font, f'dE/E0 = {dE:.5f}%',           (ax+6, y_info+48), GREEN)
        text_at(screen, font_lg, 'PAUSADO' if sim.paused else 'corriendo',
                (ax+6, y_info+72), YELLOW if sim.paused else GREEN)

        preset_lbl = font_sm.render(f'CI: {preset_name}', True, MUTED)
        screen.blit(preset_lbl, (ax+6, ay+2))

        # ===================================================
        #  PANEL B - Retratos de fase
        # ===================================================
        inner_b = draw_panel(screen, panels['phase'], '  RETRATOS DE FASE', font_sm)
        bx, by, bw, bh = inner_b

        half = bw // 2
        area_ph1 = pygame.Rect(bx,        by, half-2, bh)
        area_ph2 = pygame.Rect(bx+half+2, by, half-2, bh)

        draw_axes(screen, area_ph1)
        draw_axes(screen, area_ph2)

        # Rango fijo en X para que el retrato no se expanda con vueltas completas
        PHASE_XR = (-math.pi - 0.05, math.pi + 0.05)

        if len(sim.ph1x) > 2:
            yr = (min(sim.ph1y)-0.1, max(sim.ph1y)+0.1)
            draw_plot(screen, area_ph1, sim.ph1x, sim.ph1y, BLUE,
                      y_zero_line=True, x_range=PHASE_XR, y_range=yr)
            if yr[1] > yr[0]:
                dpx = area_ph1.x + int((sim.ph1x[-1]-PHASE_XR[0])/(PHASE_XR[1]-PHASE_XR[0])*area_ph1.width)
                dpy = area_ph1.y+area_ph1.height - int((sim.ph1y[-1]-yr[0])/(yr[1]-yr[0])*area_ph1.height)
                pygame.draw.circle(screen, WHITE, (dpx, dpy), 4)

        if len(sim.ph2x) > 2:
            yr = (min(sim.ph2y)-0.1, max(sim.ph2y)+0.1)
            draw_plot(screen, area_ph2, sim.ph2x, sim.ph2y, RED,
                      y_zero_line=True, x_range=PHASE_XR, y_range=yr)
            if yr[1] > yr[0]:
                dpx = area_ph2.x + int((sim.ph2x[-1]-PHASE_XR[0])/(PHASE_XR[1]-PHASE_XR[0])*area_ph2.width)
                dpy = area_ph2.y+area_ph2.height - int((sim.ph2y[-1]-yr[0])/(yr[1]-yr[0])*area_ph2.height)
                pygame.draw.circle(screen, WHITE, (dpx, dpy), 4)

        text_at(screen, font_sm, 'th1 in [-pi,pi] vs w1', (bx+4,       by+2), BLUE)
        text_at(screen, font_sm, 'th2 in [-pi,pi] vs w2', (bx+half+6,  by+2), RED)
        # etiquetas eje X
        for area_ph in [area_ph1, area_ph2]:
            for lbl_txt, frac in [('-pi', 0.02), ('0', 0.49), ('+pi', 0.93)]:
                lx = area_ph.x + int(frac * area_ph.width)
                text_at(screen, font, lbl_txt, (lx, area_ph.y + area_ph.height - 20), MUTED)

        # ===================================================
        #  PANEL C - Energia potencial vs total (fila inferior completa)
        # ===================================================
        inner_d = draw_panel(screen, panels['xy'],
                             '  ENERGIA: Ep(t)  vs  E_total(t)', font_sm)
        dx2, dy2, dw, dh = inner_d

        area_en = pygame.Rect(dx2, dy2, dw, dh)
        draw_axes(screen, area_en)

        if len(sim.en_t) > 2:
            pot_list = list(sim.en_pot)
            tot_list = list(sim.en_tot)
            all_e = pot_list + tot_list
            ymin_e = min(all_e) - 0.05
            ymax_e = max(all_e) + 0.05
            xr_e   = (0, max(sim.en_t) + 0.1)

            # Ep(t) — energia potencial (oscila)
            draw_plot(screen, area_en, sim.en_t, sim.en_pot, PURPLE,
                      x_range=xr_e, y_range=(ymin_e, ymax_e), lw=2)
            # E_total(t) — deberia ser casi constante (RK4 conserva bien)
            draw_plot(screen, area_en, sim.en_t, sim.en_tot, (23, 190, 207),
                      x_range=xr_e, y_range=(ymin_e, ymax_e), lw=3)

            ep_cur = sim.en_pot[-1]
            et_cur = sim.en_tot[-1]
            ek_cur = et_cur - ep_cur
            text_at(screen, font_lg, f'Ep = {ep_cur:.4f} J', (dx2+6, dy2+dh-78), PURPLE)
            text_at(screen, font_lg, f'Ek = {ek_cur:.4f} J', (dx2+6, dy2+dh-52), ORANGE)
            text_at(screen, font_lg, f'E  = {et_cur:.4f} J (total)', (dx2+6, dy2+dh-26), (23, 190, 207))

        text_at(screen, font_sm, '-- Ep potencial', (dx2+6,      dy2+2), PURPLE)
        text_at(screen, font_sm, '-- E total',      (dx2+dw//2,  dy2+2), (23, 190, 207))
        text_at(screen, font,    't (s)',            (dx2+dw-65,  dy2+dh-22), MUTED)

        # ===================================================
        #  BARRA DE AYUDA
        # ===================================================
        help_y = sh - 22
        ayuda = ('ESPACIO pausar  |  R reiniciar  |  1-4 presets  |  UP/DN k+-0.2  |  F pantalla completa  |  S screenshot  |  Q salir')
        text_at(screen, font, ayuda, (10, help_y), GRID_C)

        # Frecuencias analiticas
        f_sym  = math.sqrt(G/B) / (2*math.pi)
        f_asym = math.sqrt(G/B + 2*sim.k/M) / (2*math.pi)
        freq_str = (f'f_sim={f_sym:.3f}Hz  f_asim={f_asym:.3f}Hz  '
                    f'm={M}kg  b={B}m  k={sim.k}N/m  g={G}')
        lbl = font.render(freq_str, True, MUTED)
        screen.blit(lbl, (sw - lbl.get_width() - 8, 4))

        pygame.display.flip()
        clock.tick(FPS)

    pygame.quit()
    sys.exit()


if __name__ == '__main__':
    main()
