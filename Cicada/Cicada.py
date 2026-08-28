# ============================================================
# ANALISIS Y LIMPIEZA DE AUDIO DE CHICHARRAS (CIGARRAS)
# FFT / STFT + REDUCCION DE RUIDO + ESPECTROGRAMA
# ============================================================

import os
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt

from scipy import signal
from scipy.ndimage import gaussian_filter


# ============================================================
# 1. CONFIGURACION
# ============================================================

# Cambia este nombre por el nombre de tu archivo
ARCHIVO_AUDIO = "Cicada.wav"

# Archivos de salida
ARCHIVO_LIMPIO = "chicharras_limpio.wav"
ESPECTROGRAMA_ORIGINAL = "espectrograma_original.png"
ESPECTROGRAMA_LIMPIO = "espectrograma_limpio.png"
ESPECTROGRAMA_COMPARACION = "comparacion_espectrogramas.png"

# Duracion utilizada para estimar el ruido inicial
# Si al principio del audio hay bastante silencio/ruido ambiental,
# este valor funciona bien.
SEGUNDOS_RUIDO = 10

# Parametros STFT
N_FFT = 4096
HOP = 1024
VENTANA = "hann"

# Intensidad de reduccion del ruido
FACTOR_RUIDO = 1.5

# Piso de ganancia:
# evita que una frecuencia sea eliminada completamente.
PISO_GANANCIA = 0.08

# Suavizado de la mascara de reduccion de ruido
SUAVIZADO_FRECUENCIA = 1.0
SUAVIZADO_TIEMPO = 1.0


# ============================================================
# 2. CARGAR AUDIO
# ============================================================

print("Cargando audio...")

audio, fs = sf.read(ARCHIVO_AUDIO)

print(f"Frecuencia de muestreo: {fs} Hz")
print(f"Numero de canales: {audio.ndim}")
print(f"Duracion: {len(audio) / fs / 60:.2f} minutos")


# ============================================================
# 3. CONVERTIR A MONO
# ============================================================

if audio.ndim == 2:
    audio = np.mean(audio, axis=1)

print("Audio convertido a mono.")


# ============================================================
# 4. NORMALIZACION
# ============================================================

max_val = np.max(np.abs(audio))

if max_val > 0:
    audio = audio / max_val

print("Audio normalizado.")


# ============================================================
# 5. ESTIMAR EL PERFIL DEL RUIDO
# ============================================================

print("Estimando el ruido...")

muestras_ruido = int(SEGUNDOS_RUIDO * fs)

segmento_ruido = audio[:muestras_ruido]


frecuencias_ruido, tiempos_ruido, Z_ruido = signal.stft(
    segmento_ruido,
    fs=fs,
    window=VENTANA,
    nperseg=N_FFT,
    noverlap=N_FFT - HOP,
    boundary="zeros"
)

magnitud_ruido = np.abs(Z_ruido)

# Usamos la mediana temporal como estimacion del ruido
perfil_ruido = np.median(magnitud_ruido, axis=1)

print("Perfil de ruido calculado.")


# ============================================================
# 6. FUNCION DE REDUCCION DE RUIDO
# ============================================================

def limpiar_audio(audio_segmento, fs, perfil_ruido):

    f, t, Z = signal.stft(
        audio_segmento,
        fs=fs,
        window=VENTANA,
        nperseg=N_FFT,
        noverlap=N_FFT - HOP,
        boundary="zeros"
    )

    magnitud = np.abs(Z)

    # Evitar division por cero
    epsilon = 1e-10

    # Relacion aproximada señal/ruido
    SNR = magnitud / (FACTOR_RUIDO * perfil_ruido[:, np.newaxis] + epsilon)

    # Mascara tipo Wiener
    mascara = SNR**2 / (1 + SNR**2)

    # Evitar eliminar completamente las frecuencias
    mascara = np.maximum(mascara, PISO_GANANCIA)

    # Suavizar la mascara
    mascara = gaussian_filter(
        mascara,
        sigma=(SUAVIZADO_FRECUENCIA, SUAVIZADO_TIEMPO)
    )

    # Aplicar la mascara
    Z_limpio = Z * mascara

    # Reconstruccion mediante transformada inversa
    _, audio_limpio = signal.istft(
        Z_limpio,
        fs=fs,
        window=VENTANA,
        nperseg=N_FFT,
        noverlap=N_FFT - HOP,
        input_onesided=True,
        boundary=True
    )

    return audio_limpio


# ============================================================
# 7. PROCESAR EL AUDIO POR BLOQUES
# ============================================================

print("Iniciando limpieza del audio...")

# Para evitar consumir demasiada memoria trabajamos por bloques.
DURACION_BLOQUE = 60  # segundos

muestras_bloque = int(DURACION_BLOQUE * fs)

resultado = []

numero_bloques = int(np.ceil(len(audio) / muestras_bloque))

for i in range(numero_bloques):

    inicio = i * muestras_bloque
    fin = min((i + 1) * muestras_bloque, len(audio))

    bloque = audio[inicio:fin]

    print(
        f"Procesando bloque {i + 1}/{numero_bloques} "
        f"({inicio/fs/60:.1f} - {fin/fs/60:.1f} min)"
    )

    bloque_limpio = limpiar_audio(
        bloque,
        fs,
        perfil_ruido
    )

    # Ajustar longitud
    bloque_limpio = bloque_limpio[:len(bloque)]

    resultado.append(bloque_limpio)


# ============================================================
# 8. UNIR LOS BLOQUES
# ============================================================

audio_limpio = np.concatenate(resultado)

# Ajustar longitud final
audio_limpio = audio_limpio[:len(audio)]


# ============================================================
# 9. NORMALIZAR AUDIO LIMPIO
# ============================================================

max_limpio = np.max(np.abs(audio_limpio))

if max_limpio > 0:
    audio_limpio = audio_limpio / max_limpio * 0.95


# ============================================================
# 10. GUARDAR AUDIO LIMPIO
# ============================================================

sf.write(
    ARCHIVO_LIMPIO,
    audio_limpio,
    fs
)

print("\nAudio limpio guardado en:")
print(ARCHIVO_LIMPIO)


# ============================================================
# 11. FUNCION PARA CREAR ESPECTROGRAMA
# ============================================================

def crear_espectrograma(
    audio,
    fs,
    nombre_archivo,
    titulo,
    max_hz=None,
    duracion_maxima=None
):

    # Si se desea analizar solamente una parte
    if duracion_maxima is not None:

        muestras = int(duracion_maxima * fs)
        audio = audio[:muestras]

    # Espectrograma
    f, t, Sxx = signal.spectrogram(
        audio,
        fs=fs,
        window="hann",
        nperseg=2048,
        noverlap=1536,
        scaling="density",
        mode="magnitude"
    )

    # Convertir a dB
    Sxx_dB = 20 * np.log10(Sxx + 1e-12)

    # Limitar frecuencia
    if max_hz is not None:

        indices = f <= max_hz

        f = f[indices]
        Sxx_dB = Sxx_dB[indices, :]

    plt.figure(figsize=(14, 7))

    plt.pcolormesh(
        t,
        f,
        Sxx_dB,
        shading="auto"
    )

    plt.colorbar(label="Magnitud (dB)")

    plt.xlabel("Tiempo (s)")
    plt.ylabel("Frecuencia (Hz)")
    plt.title(titulo)

    plt.tight_layout()

    plt.savefig(
        nombre_archivo,
        dpi=300
    )

    plt.show()


# ============================================================
# 12. ESPECTROGRAMA DEL AUDIO ORIGINAL
# ============================================================

print("\nGenerando espectrograma del audio original...")

crear_espectrograma(
    audio,
    fs,
    ESPECTROGRAMA_ORIGINAL,
    "Espectrograma - Audio original",
    max_hz=20000,
    duracion_maxima=300
)


# ============================================================
# 13. ESPECTROGRAMA DEL AUDIO LIMPIO
# ============================================================

print("Generando espectrograma del audio limpio...")

crear_espectrograma(
    audio_limpio,
    fs,
    ESPECTROGRAMA_LIMPIO,
    "Espectrograma - Audio limpio",
    max_hz=20000,
    duracion_maxima=300
)


# ============================================================
# 14. COMPARACION ORIGINAL VS LIMPIO
# ============================================================

def comparar_espectrogramas(
    original,
    limpio,
    fs,
    duracion=300,
    max_hz=20000
):

    muestras = min(
        len(original),
        int(duracion * fs)
    )

    original = original[:muestras]
    limpio = limpio[:muestras]

    f1, t1, S1 = signal.spectrogram(
        original,
        fs=fs,
        window="hann",
        nperseg=2048,
        noverlap=1536,
        mode="magnitude"
    )

    f2, t2, S2 = signal.spectrogram(
        limpio,
        fs=fs,
        window="hann",
        nperseg=2048,
        noverlap=1536,
        mode="magnitude"
    )

    S1 = 20 * np.log10(S1 + 1e-12)
    S2 = 20 * np.log10(S2 + 1e-12)

    indice1 = f1 <= max_hz
    indice2 = f2 <= max_hz

    f1 = f1[indice1]
    f2 = f2[indice2]

    S1 = S1[indice1]
    S2 = S2[indice2]

    plt.figure(figsize=(15, 10))

    # --------------------------------------------------------
    # ORIGINAL
    # --------------------------------------------------------

    plt.subplot(2, 1, 1)

    plt.pcolormesh(
        t1,
        f1,
        S1,
        shading="auto"
    )

    plt.colorbar(label="Magnitud (dB)")

    plt.ylabel("Frecuencia (Hz)")
    plt.title("Audio original")

    # --------------------------------------------------------
    # LIMPIO
    # --------------------------------------------------------

    plt.subplot(2, 1, 2)

    plt.pcolormesh(
        t2,
        f2,
        S2,
        shading="auto"
    )

    plt.colorbar(label="Magnitud (dB)")

    plt.xlabel("Tiempo (s)")
    plt.ylabel("Frecuencia (Hz)")
    plt.title("Audio después de reducción de ruido")

    plt.tight_layout()

    plt.savefig(
        ESPECTROGRAMA_COMPARACION,
        dpi=300
    )

    plt.show()


comparar_espectrogramas(
    audio,
    audio_limpio,
    fs,
    duracion=300,
    max_hz=20000
)


print("\n==========================================")
print("PROCESAMIENTO FINALIZADO")
print("==========================================")
print(f"Audio limpio: {ARCHIVO_LIMPIO}")
print(f"Espectrograma original: {ESPECTROGRAMA_ORIGINAL}")
print(f"Espectrograma limpio: {ESPECTROGRAMA_LIMPIO}")
print(f"Comparacion: {ESPECTROGRAMA_COMPARACION}")