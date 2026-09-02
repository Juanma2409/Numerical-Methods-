
/*
CICADA ANALYZER - FINAL ENGINE
================================

One executable for:
  Phase 2: noise estimation + Wiener denoising + streaming ISTFT/OLA
  Phase 3: temporal/spectral features
  Phase 4: adaptive multi-feature event detection
  Multiresolution global spectrogram for Python
  Full cleaned spectrogram magnitude + phase (optional)
  Event table + event statistics
  Machine-readable JSON metadata

Design goals:
  - C++17
  - No <filesystem>
  - No std::clamp
  - Windows/MinGW friendly
  - No whole-recording RAM load
  - Streaming/block-based WAV processing
  - Exact input duration/sample count in clean WAV
  - Scientific data separated from visual summaries
  - Heuristic event detection explicitly documented as such

Default usage:
  Cicada.exe

Useful examples:
  Cicada.exe --help
  Cicada.exe --nfft 4096 --hop 1024
  Cicada.exe --noise 10 --noise-factor 1.5 --floor 0.08
  Cicada.exe --event-threshold 2.5 --event-min-ms 100 --event-merge-ms 250
  Cicada.exe --no-phase

Output (default results):
  audio_clean.wav
  noise_profile.f32
  spectrogram_clean_magnitude.f32
  spectrogram_clean_phase.f32     (unless --no-phase)
  spectrogram_global.f32
  features.f32
  events.csv
  event_features.csv
  metadata_final.json
  validation_summary.json

Notes:
  - dBFS is digital full-scale reference, not dB SPL.
  - "dominant frequency" is not "fundamental frequency".
  - F0 is intentionally NOT guessed automatically in this engine.
    Frequency tracking and harmonic inference are left to Python/event-level
    analysis where confidence and user inspection can be applied.
  - Event detection is a configurable signal-processing heuristic, not
    species identification.
*/

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace cicada {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double EPS = 1.0e-12;

static float f32(double x) {
    return static_cast<float>(x);
}

static double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static double robust_clip_z(double z) {
    if (!std::isfinite(z)) return 0.0;
    if (z < 0.0) return 0.0;
    if (z > 10.0) return 10.0;
    return z;
}

struct WavInfo {
    uint16_t audio_format = 0;   // 1 PCM, 3 IEEE float
    uint16_t channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;
    uint16_t block_align = 0;
    uint64_t data_offset = 0;
    uint64_t data_bytes = 0;
};

class WavReader {
public:
    explicit WavReader(const std::string& path)
        : file_(path.c_str(), std::ios::binary) {

        if (!file_)
            throw std::runtime_error(
                "No se pudo abrir el WAV: " + path
            );

        parse_header();

        file_.clear();
        file_.seekg(
            static_cast<std::streamoff>(info_.data_offset),
            std::ios::beg
        );
    }

    const WavInfo& info() const {
        return info_;
    }

    uint64_t total_frames() const {
        if (info_.block_align == 0) return 0;
        return info_.data_bytes / info_.block_align;
    }

    size_t read_mono(
        std::vector<float>& mono,
        size_t max_frames
    ) {
        if (max_frames == 0) return 0;

        const size_t bytes =
            max_frames * static_cast<size_t>(info_.block_align);

        raw_.resize(bytes);

        file_.read(
            reinterpret_cast<char*>(raw_.data()),
            static_cast<std::streamsize>(bytes)
        );

        const std::streamsize got = file_.gcount();

        if (got <= 0) return 0;

        const size_t frames =
            static_cast<size_t>(got) /
            static_cast<size_t>(info_.block_align);

        mono.resize(frames);

        const size_t bps =
            static_cast<size_t>(info_.bits_per_sample / 8);

        for (size_t i = 0; i < frames; ++i) {
            double sum = 0.0;

            for (uint16_t ch = 0; ch < info_.channels; ++ch) {
                const uint8_t* p =
                    raw_.data() +
                    i * static_cast<size_t>(info_.block_align) +
                    ch * bps;

                sum += decode_sample(p);
            }

            mono[i] = static_cast<float>(
                sum / static_cast<double>(info_.channels)
            );
        }

        return frames;
    }

private:
    std::ifstream file_;
    WavInfo info_;
    std::vector<uint8_t> raw_;

    static uint16_t u16(const uint8_t* p) {
        return static_cast<uint16_t>(p[0]) |
               (static_cast<uint16_t>(p[1]) << 8);
    }

    static uint32_t u32(const uint8_t* p) {
        return static_cast<uint32_t>(p[0]) |
               (static_cast<uint32_t>(p[1]) << 8) |
               (static_cast<uint32_t>(p[2]) << 16) |
               (static_cast<uint32_t>(p[3]) << 24);
    }

    static bool id_is(const uint8_t* p, const char* id) {
        return std::memcmp(p, id, 4) == 0;
    }

    void parse_header() {
        uint8_t riff[12];

        file_.read(
            reinterpret_cast<char*>(riff),
            sizeof(riff)
        );

        if (file_.gcount() !=
                static_cast<std::streamsize>(sizeof(riff)) ||
            !id_is(riff, "RIFF") ||
            !id_is(riff + 8, "WAVE")) {

            throw std::runtime_error(
                "El archivo no parece ser un WAV RIFF valido."
            );
        }

        bool got_fmt = false;
        bool got_data = false;

        while (file_) {
            uint8_t hdr[8];

            file_.read(
                reinterpret_cast<char*>(hdr),
                sizeof(hdr)
            );

            if (file_.gcount() !=
                    static_cast<std::streamsize>(sizeof(hdr)))
                break;

            const uint32_t size = u32(hdr + 4);

            const std::streamoff chunk_start =
                file_.tellg();

            if (id_is(hdr, "fmt ")) {

                if (size < 16)
                    throw std::runtime_error(
                        "Chunk fmt invalido."
                    );

                std::vector<uint8_t> fmt(size);

                file_.read(
                    reinterpret_cast<char*>(fmt.data()),
                    static_cast<std::streamsize>(size)
                );

                if (file_.gcount() !=
                        static_cast<std::streamsize>(size)) {

                    throw std::runtime_error(
                        "Chunk fmt truncado."
                    );
                }

                info_.audio_format =
                    u16(fmt.data());

                info_.channels =
                    u16(fmt.data() + 2);

                info_.sample_rate =
                    u32(fmt.data() + 4);

                info_.block_align =
                    u16(fmt.data() + 12);

                info_.bits_per_sample =
                    u16(fmt.data() + 14);

                got_fmt = true;
            }
            else if (id_is(hdr, "data")) {

                info_.data_offset =
                    static_cast<uint64_t>(chunk_start);

                info_.data_bytes =
                    static_cast<uint64_t>(size);

                file_.seekg(
                    static_cast<std::streamoff>(size),
                    std::ios::cur
                );

                got_data = true;
            }
            else {
                file_.seekg(
                    static_cast<std::streamoff>(size),
                    std::ios::cur
                );
            }

            if (size & 1U)
                file_.seekg(1, std::ios::cur);

            if (got_fmt && got_data)
                break;
        }

        if (!got_fmt || !got_data)
            throw std::runtime_error(
                "WAV sin chunks fmt/data validos."
            );

        if (info_.channels == 0 ||
            info_.sample_rate == 0 ||
            info_.block_align == 0) {

            throw std::runtime_error(
                "Parametros WAV invalidos."
            );
        }

        if (info_.audio_format != 1 &&
            info_.audio_format != 3) {

            throw std::runtime_error(
                "Solo se admiten WAV PCM (1) o IEEE float (3)."
            );
        }

        if (info_.audio_format == 1 &&
            !(info_.bits_per_sample == 8 ||
              info_.bits_per_sample == 16 ||
              info_.bits_per_sample == 24 ||
              info_.bits_per_sample == 32)) {

            throw std::runtime_error(
                "PCM: profundidad de bits no soportada."
            );
        }

        if (info_.audio_format == 3 &&
            info_.bits_per_sample != 32) {

            throw std::runtime_error(
                "Float WAV: solo float32."
            );
        }
    }

    double decode_sample(const uint8_t* p) const {

        if (info_.audio_format == 3) {
            float v;
            std::memcpy(&v, p, sizeof(float));

            double x = static_cast<double>(v);

            if (x < -1.0) x = -1.0;
            if (x > 1.0) x = 1.0;

            return x;
        }

        switch (info_.bits_per_sample) {

            case 8:
                return (
                    static_cast<double>(p[0]) -
                    128.0
                ) / 128.0;

            case 16: {
                const int16_t v =
                    static_cast<int16_t>(u16(p));

                return static_cast<double>(v) /
                       32768.0;
            }

            case 24: {
                int32_t v =
                    static_cast<int32_t>(p[0]) |
                    (static_cast<int32_t>(p[1]) << 8) |
                    (static_cast<int32_t>(p[2]) << 16);

                if (v & 0x00800000)
                    v |= static_cast<int32_t>(0xFF000000);

                return static_cast<double>(v) /
                       8388608.0;
            }

            case 32: {
                const int32_t v =
                    static_cast<int32_t>(u32(p));

                return static_cast<double>(v) /
                       2147483648.0;
            }

            default:
                throw std::runtime_error(
                    "Bits por muestra no soportados."
                );
        }
    }
};

class WavMono16Writer {
public:
    WavMono16Writer(
        const std::string& path,
        uint32_t sample_rate,
        uint64_t expected_frames
    )
        : file_(path.c_str(), std::ios::binary),
          sample_rate_(sample_rate),
          expected_frames_(expected_frames) {

        if (!file_)
            throw std::runtime_error(
                "No se pudo crear WAV de salida: " + path
            );

        write_header_placeholder();
    }

    void write_samples(
        const std::vector<float>& x,
        size_t count
    ) {
        if (count > x.size())
            throw std::runtime_error(
                "Cantidad de muestras invalida."
            );

        int16_t buf[4096];

        size_t pos = 0;

        while (pos < count) {
            const size_t n =
                std::min<size_t>(
                    4096,
                    count - pos
                );

            for (size_t i = 0; i < n; ++i) {
                double v =
                    static_cast<double>(x[pos + i]);

                if (v < -1.0) v = -1.0;
                if (v > 1.0) v = 1.0;

                if (v >= 1.0) {
                    buf[i] = 32767;
                }
                else if (v <= -1.0) {
                    buf[i] = -32768;
                }
                else {
                    buf[i] = static_cast<int16_t>(
                        std::lrint(v * 32767.0)
                    );
                }
            }

            file_.write(
                reinterpret_cast<const char*>(buf),
                static_cast<std::streamsize>(
                    n * sizeof(int16_t)
                )
            );

            if (!file_)
                throw std::runtime_error(
                    "Error escribiendo WAV de salida."
                );

            samples_written_ += n;
            pos += n;
        }
    }

    uint64_t samples_written() const {
        return samples_written_;
    }

    void close() {
        if (closed_) return;

        if (samples_written_ != expected_frames_)
            throw std::runtime_error(
                "El WAV de salida no tiene exactamente "
                "la cantidad esperada de muestras."
            );

        const uint32_t data_bytes =
            static_cast<uint32_t>(
                samples_written_ * sizeof(int16_t)
            );

        file_.seekp(40, std::ios::beg);
        write_u32(data_bytes);

        const uint32_t riff_size =
            36u + data_bytes;

        file_.seekp(4, std::ios::beg);
        write_u32(riff_size);

        file_.flush();
        file_.close();

        closed_ = true;
    }

    ~WavMono16Writer() {
        try {
            if (!closed_) close();
        }
        catch (...) {
        }
    }

private:
    std::ofstream file_;
    uint32_t sample_rate_;
    uint64_t expected_frames_;
    uint64_t samples_written_ = 0;
    bool closed_ = false;

    void write_u16(uint16_t v) {
        file_.put(static_cast<char>(v & 0xFF));
        file_.put(static_cast<char>((v >> 8) & 0xFF));
    }

    void write_u32(uint32_t v) {
        file_.write(
            reinterpret_cast<const char*>(&v),
            sizeof(uint32_t)
        );

        // RIFF is little-endian. Correct byte writing is explicit
        // in the placeholder header below; this function is only used
        // after seek and on Windows/MinGW little-endian targets.
    }

    void write_header_placeholder() {
        // PCM mono 16-bit RIFF/WAVE, 44-byte header.
        file_.write("RIFF", 4);
        write_le32(0);

        file_.write("WAVE", 4);

        file_.write("fmt ", 4);
        write_le32(16);
        write_le16(1);       // PCM
        write_le16(1);       // mono
        write_le32(sample_rate_);
        write_le32(sample_rate_ * 2); // byte rate
        write_le16(2);       // block align
        write_le16(16);      // bits

        file_.write("data", 4);
        write_le32(0);
    }

    void write_le16(uint16_t v) {
        file_.put(static_cast<char>(v & 0xFF));
        file_.put(static_cast<char>((v >> 8) & 0xFF));
    }

    void write_le32(uint32_t v) {
        file_.put(static_cast<char>(v & 0xFF));
        file_.put(static_cast<char>((v >> 8) & 0xFF));
        file_.put(static_cast<char>((v >> 16) & 0xFF));
        file_.put(static_cast<char>((v >> 24) & 0xFF));
    }
};

struct Complex {
    double re = 0.0;
    double im = 0.0;
};

static Complex c_add(Complex a, Complex b) {
    return {a.re + b.re, a.im + b.im};
}

static Complex c_sub(Complex a, Complex b) {
    return {a.re - b.re, a.im - b.im};
}

static Complex c_mul(Complex a, Complex b) {
    return {
        a.re * b.re - a.im * b.im,
        a.re * b.im + a.im * b.re
    };
}

static Complex c_conj(Complex a) {
    return {a.re, -a.im};
}

void fft(
    std::vector<Complex>& a,
    bool inverse = false
) {
    const size_t n = a.size();

    if (n < 2 ||
        (n & (n - 1)) != 0) {

        throw std::runtime_error(
            "NFFT debe ser potencia de 2."
        );
    }

    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;

        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }

        j ^= bit;

        if (i < j)
            std::swap(a[i], a[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {

        const double sign =
            inverse ? 1.0 : -1.0;

        const double angle =
            sign * 2.0 * PI /
            static_cast<double>(len);

        const Complex wlen{
            std::cos(angle),
            std::sin(angle)
        };

        for (size_t i = 0; i < n; i += len) {

            Complex w{1.0, 0.0};

            const size_t half = len >> 1;

            for (size_t j = 0; j < half; ++j) {

                const Complex u =
                    a[i + j];

                const Complex v =
                    c_mul(
                        a[i + j + half],
                        w
                    );

                a[i + j] =
                    c_add(u, v);

                a[i + j + half] =
                    c_sub(u, v);

                w = c_mul(w, wlen);
            }
        }
    }

    if (inverse) {
        const double inv_n =
            1.0 / static_cast<double>(n);

        for (size_t i = 0; i < n; ++i) {
            a[i].re *= inv_n;
            a[i].im *= inv_n;
        }
    }
}

static std::vector<float> hann_window(size_t n) {
    std::vector<float> w(n);

    if (n == 1) {
        w[0] = 1.0f;
        return w;
    }

    for (size_t i = 0; i < n; ++i) {
        w[i] = static_cast<float>(
            0.5 *
            (
                1.0 -
                std::cos(
                    2.0 * PI *
                    static_cast<double>(i) /
                    static_cast<double>(n - 1)
                )
            )
        );
    }

    return w;
}

static bool directory_exists(const std::string& path) {
#ifdef _WIN32
    const DWORD attr =
        GetFileAttributesA(path.c_str());

    return attr != INVALID_FILE_ATTRIBUTES &&
           (attr & FILE_ATTRIBUTE_DIRECTORY) != 0;
#else
    std::ifstream f(path.c_str());
    (void)f;
    return false;
#endif
}

static void ensure_directory(const std::string& path) {
    if (directory_exists(path))
        return;

#ifdef _WIN32
    if (CreateDirectoryA(path.c_str(), NULL))
        return;

    const DWORD e = GetLastError();

    if (e == ERROR_ALREADY_EXISTS &&
        directory_exists(path))
        return;

    throw std::runtime_error(
        "No se pudo crear la carpeta: " + path
    );
#else
    throw std::runtime_error(
        "Esta compilacion esta orientada a Windows."
    );
#endif
}

static std::string executable_directory() {
#ifdef _WIN32
    char buffer[MAX_PATH];

    const DWORD len =
        GetModuleFileNameA(
            NULL,
            buffer,
            MAX_PATH
        );

    if (len == 0)
        throw std::runtime_error(
            "No se pudo determinar la ubicacion del ejecutable."
        );

    std::string p(buffer, len);

    const size_t pos =
        p.find_last_of("\\/");

    if (pos == std::string::npos)
        return ".";

    return p.substr(0, pos);
#else
    return ".";
#endif
}

struct Config {
    std::string input;
    std::string output;

    size_t nfft = 2048;
    size_t hop = 1024;
    double block_seconds = 20.0;

    double noise_seconds = 10.0;
    double noise_factor = 1.5;
    double gain_floor = 0.08;

    double fmin_hz = 100.0;
    double fmax_hz = 12000.0;
    double rolloff_fraction = 0.85;

    size_t global_time_bins = 4096;
    size_t global_freq_bins = 512;

    double event_threshold = 2.5;
    double event_min_ms = 100.0;
    double event_merge_ms = 250.0;

    bool save_phase = true;
};

static bool parse_size(
    const char* s,
    size_t& value
) {
    try {
        const unsigned long long x =
            std::stoull(s);

        if (x == 0)
            return false;

        value = static_cast<size_t>(x);
        return true;
    }
    catch (...) {
        return false;
    }
}

static bool parse_double(
    const char* s,
    double& value
) {
    try {
        value = std::stod(s);
        return std::isfinite(value);
    }
    catch (...) {
        return false;
    }
}

static void usage() {
    std::cout
        << "CICADA ANALYZER - FINAL ENGINE\n\n"
        << "Uso:\n"
        << "  Cicada.exe\n"
        << "  Cicada.exe --input otro.wav\n\n"
        << "Entrada/salida:\n"
        << "  --input <archivo.wav>\n"
        << "  --output <carpeta>\n\n"
        << "STFT:\n"
        << "  --nfft <entero>             default 2048\n"
        << "  --hop <entero>              default 1024\n"
        << "  --block <segundos>          default 20\n\n"
        << "Reduccion Wiener:\n"
        << "  --noise <segundos>          default 10\n"
        << "  --noise-factor <x>          default 1.5\n"
        << "  --floor <x>                 default 0.08\n\n"
        << "Analisis espectral:\n"
        << "  --fmin <Hz>                 default 100\n"
        << "  --fmax <Hz>                 default 12000\n"
        << "  --rolloff <0-1>             default 0.85\n\n"
        << "Espectrograma global:\n"
        << "  --global-time-bins <N>      default 4096\n"
        << "  --global-freq-bins <N>      default 512\n\n"
        << "Deteccion de eventos:\n"
        << "  --event-threshold <x>       default 2.5\n"
        << "  --event-min-ms <ms>         default 100\n"
        << "  --event-merge-ms <ms>       default 250\n\n"
        << "Salida:\n"
        << "  --no-phase                  no guardar fase STFT completa\n"
        << "  --help\n\n"
        << "El motor no asigna automaticamente frecuencia fundamental (F0):\n"
        << "frecuencia dominante != frecuencia fundamental.\n";
}

struct FeatureFrame {
    float time_seconds = 0.0f;

    float rms = 0.0f;
    float rms_dbfs = -240.0f;
    float peak = 0.0f;
    float crest_factor = 0.0f;

    float dominant_frequency_hz = 0.0f;
    float spectral_centroid_hz = 0.0f;
    float spectral_bandwidth_hz = 0.0f;
    float spectral_rolloff_hz = 0.0f;
    float spectral_flatness = 0.0f;
    float spectral_flux = 0.0f;

    float band_0_500_ratio = 0.0f;
    float band_500_1000_ratio = 0.0f;
    float band_1000_2000_ratio = 0.0f;
    float band_2000_4000_ratio = 0.0f;
    float band_4000_8000_ratio = 0.0f;
    float band_8000_12000_ratio = 0.0f;

    float spectral_power = 0.0f;
};

static const int FEATURE_COUNT = 18;

static void write_feature(
    std::ofstream& out,
    const FeatureFrame& f
) {
    const float values[FEATURE_COUNT] = {
        f.time_seconds,
        f.rms,
        f.rms_dbfs,
        f.peak,
        f.crest_factor,
        f.dominant_frequency_hz,
        f.spectral_centroid_hz,
        f.spectral_bandwidth_hz,
        f.spectral_rolloff_hz,
        f.spectral_flatness,
        f.spectral_flux,
        f.band_0_500_ratio,
        f.band_500_1000_ratio,
        f.band_1000_2000_ratio,
        f.band_2000_4000_ratio,
        f.band_4000_8000_ratio,
        f.band_8000_12000_ratio,
        f.spectral_power
    };

    out.write(
        reinterpret_cast<const char*>(values),
        sizeof(values)
    );

    if (!out)
        throw std::runtime_error(
            "Error escribiendo features.f32."
        );
}

static double spectral_band_power(
    const std::vector<double>& mag,
    uint32_t fs,
    size_t nfft,
    double lo,
    double hi,
    double total_power
) {
    if (total_power <= EPS)
        return 0.0;

    const double df =
        static_cast<double>(fs) /
        static_cast<double>(nfft);

    size_t k0 =
        static_cast<size_t>(
            std::ceil(lo / df)
        );

    size_t k1 =
        static_cast<size_t>(
            std::floor(hi / df)
        );

    const size_t max_bin =
        nfft / 2;

    if (k0 > max_bin)
        return 0.0;

    if (k1 > max_bin)
        k1 = max_bin;

    if (k1 < k0)
        return 0.0;

    double p = 0.0;

    for (size_t k = k0; k <= k1; ++k) {

        double q =
            mag[k] * mag[k];

        if (k != 0 && k != max_bin)
            q *= 2.0;

        p += q;
    }

    return clamp01(p / total_power);
}

static FeatureFrame calculate_features(
    const std::vector<float>& frame,
    const std::vector<float>& window,
    std::vector<Complex>& spectrum,
    std::vector<double>& magnitude,
    const std::vector<double>& previous_norm,
    bool have_previous,
    std::vector<double>& current_norm,
    uint32_t fs,
    size_t nfft,
    uint64_t frame_start,
    double fmin,
    double fmax,
    double rolloff
) {
    double sum_sq = 0.0;
    double peak_abs = 0.0;

    for (size_t i = 0; i < nfft; ++i) {

        const double x =
            static_cast<double>(frame[i]);

        sum_sq += x * x;

        peak_abs =
            std::max(
                peak_abs,
                std::abs(x)
            );

        spectrum[i] = {
            x * static_cast<double>(window[i]),
            0.0
        };
    }

    fft(spectrum, false);

    const size_t bins =
        nfft / 2 + 1;

    for (size_t k = 0; k < bins; ++k) {

        magnitude[k] =
            std::hypot(
                spectrum[k].re,
                spectrum[k].im
            );
    }

    const double rms =
        std::sqrt(
            sum_sq /
            static_cast<double>(nfft)
        );

    const double df =
        static_cast<double>(fs) /
        static_cast<double>(nfft);

    const double nyquist =
        static_cast<double>(fs) / 2.0;

    const double lo =
        std::max(
            0.0,
            std::min(
                fmin,
                nyquist
            )
        );

    const double hi =
        std::max(
            lo,
            std::min(
                fmax,
                nyquist
            )
        );

    size_t k0 =
        static_cast<size_t>(
            std::ceil(lo / df)
        );

    size_t k1 =
        static_cast<size_t>(
            std::floor(hi / df)
        );

    if (k0 >= bins) k0 = bins - 1;
    if (k1 >= bins) k1 = bins - 1;
    if (k1 < k0) k1 = k0;

    double sum_mag = 0.0;
    double sum_f_mag = 0.0;
    double sum_f2_mag = 0.0;
    double total_power = 0.0;

    double max_mag = -1.0;
    size_t dominant_bin = k0;

    for (size_t k = k0; k <= k1; ++k) {

        const double m = magnitude[k];
        const double f =
            static_cast<double>(k) * df;

        sum_mag += m;
        sum_f_mag += f * m;
        sum_f2_mag += f * f * m;

        double p = m * m;

        if (k != 0 && k != nfft / 2)
            p *= 2.0;

        total_power += p;

        if (m > max_mag) {
            max_mag = m;
            dominant_bin = k;
        }
    }

    const double centroid =
        sum_mag > EPS
        ? sum_f_mag / sum_mag
        : 0.0;

    double variance = 0.0;

    if (sum_mag > EPS) {

        variance =
            sum_f2_mag / sum_mag -
            centroid * centroid;

        if (variance < 0.0)
            variance = 0.0;
    }

    const double bandwidth =
        std::sqrt(variance);

    const double target =
        total_power *
        clamp01(rolloff);

    double cumulative = 0.0;
    double rolloff_hz = lo;

    for (size_t k = k0; k <= k1; ++k) {

        double p =
            magnitude[k] *
            magnitude[k];

        if (k != 0 && k != nfft / 2)
            p *= 2.0;

        cumulative += p;

        if (cumulative >= target) {

            rolloff_hz =
                static_cast<double>(k) * df;

            break;
        }
    }

    const double n_analysis =
        static_cast<double>(k1 - k0 + 1);

    double log_sum = 0.0;
    double arithmetic_sum = 0.0;

    for (size_t k = k0; k <= k1; ++k) {

        const double m =
            std::max(
                magnitude[k],
                1.0e-20
            );

        log_sum += std::log(m);
        arithmetic_sum += m;
    }

    const double geom =
        std::exp(
            log_sum /
            std::max(
                n_analysis,
                1.0
            )
        );

    const double arith =
        arithmetic_sum /
        std::max(
            n_analysis,
            1.0
        );

    const double flatness =
        arith > EPS
        ? clamp01(geom / arith)
        : 0.0;

    current_norm.assign(bins, 0.0);

    double norm_sq = 0.0;

    for (size_t k = k0; k <= k1; ++k)
        norm_sq +=
            magnitude[k] *
            magnitude[k];

    const double spec_norm =
        std::sqrt(norm_sq);

    if (spec_norm > EPS) {

        for (size_t k = k0; k <= k1; ++k)
            current_norm[k] =
                magnitude[k] /
                spec_norm;
    }

    double flux = 0.0;

    if (have_previous &&
        previous_norm.size() == current_norm.size()) {

        double q = 0.0;

        for (size_t k = k0; k <= k1; ++k) {

            const double d =
                current_norm[k] -
                previous_norm[k];

            q += d * d;
        }

        flux = std::sqrt(q);
    }

    const double spectral_power =
        total_power /
        static_cast<double>(nfft * nfft);

    const double rms_dbfs =
        20.0 *
        std::log10(
            std::max(
                rms,
                1.0e-12
            )
        );

    const double crest =
        rms > EPS
        ? peak_abs / rms
        : 0.0;

    const double b0 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            0.0,
            500.0,
            total_power
        );

    const double b1 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            500.0,
            1000.0,
            total_power
        );

    const double b2 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            1000.0,
            2000.0,
            total_power
        );

    const double b3 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            2000.0,
            4000.0,
            total_power
        );

    const double b4 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            4000.0,
            8000.0,
            total_power
        );

    const double b5 =
        spectral_band_power(
            magnitude,
            fs,
            nfft,
            8000.0,
            12000.0,
            total_power
        );

    return {
        static_cast<float>(
            (
                static_cast<double>(frame_start) +
                0.5 * static_cast<double>(nfft)
            ) /
            static_cast<double>(fs)
        ),

        f32(rms),
        f32(rms_dbfs),
        f32(peak_abs),
        f32(crest),

        f32(
            static_cast<double>(dominant_bin) *
            df
        ),

        f32(centroid),
        f32(bandwidth),
        f32(rolloff_hz),
        f32(flatness),
        f32(flux),

        f32(b0),
        f32(b1),
        f32(b2),
        f32(b3),
        f32(b4),
        f32(b5),

        f32(spectral_power)
    };
}

struct RobustStats {
    double median = 0.0;
    double mad = 0.0;
    double scale = 1.0;
};

static double median_copy(
    const std::vector<double>& values
) {
    if (values.empty())
        return 0.0;

    std::vector<double> x(values);

    const size_t n = x.size();
    const size_t mid = n / 2;

    std::nth_element(
        x.begin(),
        x.begin() + static_cast<std::ptrdiff_t>(mid),
        x.end()
    );

    if ((n & 1U) != 0)
        return x[mid];

    const double hi = x[mid];

    std::nth_element(
        x.begin(),
        x.begin() +
            static_cast<std::ptrdiff_t>(mid - 1),
        x.end()
    );

    const double lo = x[mid - 1];

    return 0.5 * (lo + hi);
}

static RobustStats robust_statistics(
    const std::vector<double>& x
) {
    if (x.empty())
        return {};

    const double med =
        median_copy(x);

    std::vector<double> deviations(x.size());

    for (size_t i = 0; i < x.size(); ++i)
        deviations[i] =
            std::abs(x[i] - med);

    const double mad =
        median_copy(deviations);

    // 1.4826 * MAD is a robust estimate of sigma for Gaussian noise.
    const double scale =
        std::max(
            1.4826 * mad,
            1.0e-9
        );

    return {
        med,
        mad,
        scale
    };
}

struct Event {
    uint64_t first_frame = 0;
    uint64_t last_frame = 0;

    double start_seconds = 0.0;
    double end_seconds = 0.0;
    double duration_seconds = 0.0;

    double peak = 0.0;
    double rms = 0.0;
    double rms_dbfs = -240.0;

    double dominant_frequency_hz = 0.0;
    double dominant_frequency_median_hz = 0.0;
    double dominant_frequency_slope_hz_per_s = 0.0;

    double centroid_hz = 0.0;
    double bandwidth_hz = 0.0;
    double rolloff_hz = 0.0;
    double flatness = 0.0;
    double flux = 0.0;

    double band_0_500 = 0.0;
    double band_500_1000 = 0.0;
    double band_1000_2000 = 0.0;
    double band_2000_4000 = 0.0;
    double band_4000_8000 = 0.0;
    double band_8000_12000 = 0.0;

    double spectral_power = 0.0;

    double attack_seconds = 0.0;
    double decay_seconds = 0.0;

    double detection_score = 0.0;
};

static double median_feature(
    const std::vector<FeatureFrame>& f,
    size_t first,
    size_t last,
    int which
) {
    if (first > last ||
        last >= f.size())
        return 0.0;

    std::vector<double> x;
    x.reserve(last - first + 1);

    for (size_t i = first; i <= last; ++i) {

        double v = 0.0;

        switch (which) {
            case 0: v = f[i].rms_dbfs; break;
            case 1: v = f[i].dominant_frequency_hz; break;
            case 2: v = f[i].spectral_centroid_hz; break;
            case 3: v = f[i].spectral_bandwidth_hz; break;
            case 4: v = f[i].spectral_rolloff_hz; break;
            case 5: v = f[i].spectral_flatness; break;
            case 6: v = f[i].spectral_flux; break;
            case 7: v = f[i].band_0_500_ratio; break;
            case 8: v = f[i].band_500_1000_ratio; break;
            case 9: v = f[i].band_1000_2000_ratio; break;
            case 10: v = f[i].band_2000_4000_ratio; break;
            case 11: v = f[i].band_4000_8000_ratio; break;
            case 12: v = f[i].band_8000_12000_ratio; break;
            default: return 0.0;
        }

        if (std::isfinite(v))
            x.push_back(v);
    }

    return median_copy(x);
}

static Event summarize_event(
    const std::vector<FeatureFrame>& features,
    size_t first,
    size_t last,
    double score
) {
    Event e;

    e.first_frame = first;
    e.last_frame = last;

    e.start_seconds =
        static_cast<double>(
            features[first].time_seconds
        );

    e.end_seconds =
        static_cast<double>(
            features[last].time_seconds
        );

    // Extend each frame center to a time interval approximation.
    e.duration_seconds =
        std::max(
            0.0,
            e.end_seconds -
            e.start_seconds
        );

    e.detection_score = score;

    e.rms_dbfs =
        median_feature(
            features,
            first,
            last,
            0
        );

    e.dominant_frequency_median_hz =
        median_feature(
            features,
            first,
            last,
            1
        );

    e.dominant_frequency_hz =
        e.dominant_frequency_median_hz;

    e.centroid_hz =
        median_feature(
            features,
            first,
            last,
            2
        );

    e.bandwidth_hz =
        median_feature(
            features,
            first,
            last,
            3
        );

    e.rolloff_hz =
        median_feature(
            features,
            first,
            last,
            4
        );

    e.flatness =
        median_feature(
            features,
            first,
            last,
            5
        );

    e.flux =
        median_feature(
            features,
            first,
            last,
            6
        );

    e.band_0_500 =
        median_feature(
            features,
            first,
            last,
            7
        );

    e.band_500_1000 =
        median_feature(
            features,
            first,
            last,
            8
        );

    e.band_1000_2000 =
        median_feature(
            features,
            first,
            last,
            9
        );

    e.band_2000_4000 =
        median_feature(
            features,
            first,
            last,
            10
        );

    e.band_4000_8000 =
        median_feature(
            features,
            first,
            last,
            11
        );

    e.band_8000_12000 =
        median_feature(
            features,
            first,
            last,
            12
        );

    double sum_rms_sq = 0.0;
    double peak = 0.0;
    double sum_power = 0.0;

    for (size_t i = first; i <= last; ++i) {

        const double r =
            features[i].rms;

        sum_rms_sq += r * r;

        peak =
            std::max(
                peak,
                static_cast<double>(
                    features[i].peak
                )
            );

        sum_power +=
            features[i].spectral_power;
    }

    const size_t count =
        last - first + 1;

    e.rms =
        std::sqrt(
            sum_rms_sq /
            std::max<size_t>(1, count)
        );

    e.peak = peak;

    e.spectral_power =
        sum_power /
        static_cast<double>(
            std::max<size_t>(1, count)
        );

    // Frequency slope by ordinary least squares.
    if (count >= 2) {

        double sum_t = 0.0;
        double sum_f = 0.0;
        double sum_tt = 0.0;
        double sum_tf = 0.0;

        for (size_t i = first; i <= last; ++i) {

            const double t =
                features[i].time_seconds;

            const double f =
                features[i].dominant_frequency_hz;

            sum_t += t;
            sum_f += f;
            sum_tt += t * t;
            sum_tf += t * f;
        }

        const double n =
            static_cast<double>(count);

        const double denom =
            n * sum_tt -
            sum_t * sum_t;

        if (std::abs(denom) > EPS)
            e.dominant_frequency_slope_hz_per_s =
                (
                    n * sum_tf -
                    sum_t * sum_f
                ) / denom;
    }

    // Attack / decay from RMS envelope inside the event.
    double max_rms = 0.0;
    size_t max_index = first;

    for (size_t i = first; i <= last; ++i) {

        if (features[i].rms > max_rms) {
            max_rms =
                features[i].rms;

            max_index = i;
        }
    }

    if (max_rms > EPS) {

        const double threshold10 =
            0.10 * max_rms;

        const double threshold90 =
            0.90 * max_rms;

        size_t attack_start = first;

        for (size_t i = first;
             i <= max_index;
             ++i) {

            if (features[i].rms >= threshold10) {
                attack_start = i;
                break;
            }
        }

        size_t attack90 = attack_start;

        for (size_t i = attack_start;
             i <= max_index;
             ++i) {

            if (features[i].rms >= threshold90) {
                attack90 = i;
                break;
            }
        }

        e.attack_seconds =
            std::max(
                0.0,
                static_cast<double>(
                    features[attack90].time_seconds -
                    features[attack_start].time_seconds
                )
            );

        size_t decay90 = max_index;

        for (size_t i = max_index;
             i <= last;
             ++i) {

            if (features[i].rms <= threshold90) {
                decay90 = i;
                break;
            }
        }

        size_t decay10 = decay90;

        for (size_t i = decay90;
             i <= last;
             ++i) {

            if (features[i].rms <= threshold10) {
                decay10 = i;
                break;
            }
        }

        e.decay_seconds =
            std::max(
                0.0,
                static_cast<double>(
                    features[decay10].time_seconds -
                    features[decay90].time_seconds
                )
            );
    }

    return e;
}

static std::vector<Event> detect_events(
    const std::vector<FeatureFrame>& features,
    uint32_t fs,
    size_t hop,
    double threshold,
    double min_ms,
    double merge_ms
) {
    std::vector<Event> events;

    if (features.empty())
        return events;

    std::vector<double> rms_db;
    std::vector<double> band_4_8;
    std::vector<double> flux;
    std::vector<double> inv_flatness;

    rms_db.reserve(features.size());
    band_4_8.reserve(features.size());
    flux.reserve(features.size());
    inv_flatness.reserve(features.size());

    for (const FeatureFrame& f : features) {

        rms_db.push_back(
            static_cast<double>(f.rms_dbfs)
        );

        band_4_8.push_back(
            std::log10(
                std::max(
                    1.0e-12,
                    static_cast<double>(
                        f.band_4000_8000_ratio
                    )
                )
            )
        );

        flux.push_back(
            static_cast<double>(f.spectral_flux)
        );

        inv_flatness.push_back(
            1.0 -
            clamp01(
                static_cast<double>(
                    f.spectral_flatness
                )
            )
        );
    }

    const RobustStats s_rms =
        robust_statistics(rms_db);

    const RobustStats s_band =
        robust_statistics(band_4_8);

    const RobustStats s_flux =
        robust_statistics(flux);

    const RobustStats s_flat =
        robust_statistics(inv_flatness);

    std::vector<double> score(features.size(), 0.0);
    std::vector<uint8_t> active(features.size(), 0);

    for (size_t i = 0; i < features.size(); ++i) {

        const double z_rms =
            robust_clip_z(
                (
                    rms_db[i] -
                    s_rms.median
                ) /
                s_rms.scale
            );

        const double z_band =
            robust_clip_z(
                (
                    band_4_8[i] -
                    s_band.median
                ) /
                s_band.scale
            );

        const double z_flux =
            robust_clip_z(
                (
                    flux[i] -
                    s_flux.median
                ) /
                s_flux.scale
            );

        const double z_flat =
            robust_clip_z(
                (
                    inv_flatness[i] -
                    s_flat.median
                ) /
                s_flat.scale
            );

        // Multi-feature heuristic:
        // amplitude is primary, spectral concentration supports it,
        // spectral flux helps find transitions.
        const double combined =
            0.50 * z_rms +
            0.25 * z_band +
            0.15 * z_flux +
            0.10 * z_flat;

        score[i] = combined;

        if (combined >= threshold)
            active[i] = 1;
    }

    // Temporal merging: active frames are converted into contiguous episodes.
    size_t i = 0;

    while (i < active.size()) {

        while (i < active.size() &&
               active[i] == 0)
            ++i;

        if (i >= active.size())
            break;

        const size_t start = i;

        while (i < active.size() &&
               active[i] != 0)
            ++i;

        const size_t end = i - 1;

        const double duration =
            static_cast<double>(
                features[end].time_seconds -
                features[start].time_seconds
            ) +
            static_cast<double>(hop) /
            static_cast<double>(fs);

        if (duration >= min_ms / 1000.0) {

            const double gap_limit =
                merge_ms / 1000.0;

            size_t merge_end = end;

            // Merge following events if the silent gap is short.
            size_t j = end + 1;

            while (j < active.size()) {

                while (j < active.size() &&
                       active[j] == 0)
                    ++j;

                if (j >= active.size())
                    break;

                const size_t next_start = j;

                while (j < active.size() &&
                       active[j] != 0)
                    ++j;

                const size_t next_end = j - 1;

                const double gap =
                    features[next_start].time_seconds -
                    features[merge_end].time_seconds;

                if (gap <= gap_limit) {

                    merge_end = next_end;
                    continue;
                }

                break;
            }

            double max_score = 0.0;

            for (size_t k = start;
                 k <= merge_end;
                 ++k) {

                max_score =
                    std::max(
                        max_score,
                        score[k]
                    );
            }

            events.push_back(
                summarize_event(
                    features,
                    start,
                    merge_end,
                    max_score
                )
            );

            i = merge_end + 1;
        }
    }

    // Inter-event interval.
    return events;
}

static void write_metadata(
    const std::string& path,
    const WavInfo& wi,
    uint64_t total_frames,
    size_t feature_frames,
    size_t nfft,
    size_t hop,
    const Config& cfg,
    bool phase_saved,
    double dc_offset
) {
    std::ofstream out(path.c_str());

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir metadata_final.json."
        );

    const size_t bins =
        nfft / 2 + 1;

    const size_t global_count =
        cfg.global_time_bins *
        cfg.global_freq_bins;

    out
        << "{\n"
        << "  \"format\": \"cicada-final-engine-v1\",\n"
        << "  \"input_sample_rate\": "
        << wi.sample_rate << ",\n"
        << "  \"input_channels\": "
        << wi.channels << ",\n"
        << "  \"input_bits_per_sample\": "
        << wi.bits_per_sample << ",\n"
        << "  \"input_duration_seconds\": "
        << std::setprecision(15)
        << (
            static_cast<double>(
                total_frames
            ) /
            static_cast<double>(
                wi.sample_rate
            )
        )
        << ",\n"
        << "  \"mono_processing\": true,\n"
        << "  \"dc_offset_removed\": "
        << dc_offset << ",\n"
        << "  \"nfft\": "
        << nfft << ",\n"
        << "  \"hop\": "
        << hop << ",\n"
        << "  \"window\": \"hann_symmetric\",\n"
        << "  \"frequency_resolution_hz\": "
        << (
            static_cast<double>(
                wi.sample_rate
            ) /
            static_cast<double>(
                nfft
            )
        )
        << ",\n"
        << "  \"time_step_seconds\": "
        << (
            static_cast<double>(
                hop
            ) /
            static_cast<double>(
                wi.sample_rate
            )
        )
        << ",\n"
        << "  \"analysis_fmin_hz\": "
        << cfg.fmin_hz << ",\n"
        << "  \"analysis_fmax_hz\": "
        << cfg.fmax_hz << ",\n"
        << "  \"wiener_noise_seconds\": "
        << cfg.noise_seconds << ",\n"
        << "  \"wiener_noise_factor\": "
        << cfg.noise_factor << ",\n"
        << "  \"wiener_gain_floor\": "
        << cfg.gain_floor << ",\n"
        << "  \"rolloff_fraction\": "
        << cfg.rolloff_fraction << ",\n"
        << "  \"feature_frames\": "
        << feature_frames << ",\n"
        << "  \"feature_count\": "
        << FEATURE_COUNT << ",\n"
        << "  \"feature_dtype\": \"float32\",\n"
        << "  \"feature_layout\": \"row-major [frame, feature]\",\n"
        << "  \"feature_names\": [\n"
        << "    \"time_seconds\",\n"
        << "    \"rms\",\n"
        << "    \"rms_dbfs\",\n"
        << "    \"peak\",\n"
        << "    \"crest_factor\",\n"
        << "    \"dominant_frequency_hz\",\n"
        << "    \"spectral_centroid_hz\",\n"
        << "    \"spectral_bandwidth_hz\",\n"
        << "    \"spectral_rolloff_hz\",\n"
        << "    \"spectral_flatness\",\n"
        << "    \"spectral_flux\",\n"
        << "    \"band_0_500_ratio\",\n"
        << "    \"band_500_1000_ratio\",\n"
        << "    \"band_1000_2000_ratio\",\n"
        << "    \"band_2000_4000_ratio\",\n"
        << "    \"band_4000_8000_ratio\",\n"
        << "    \"band_8000_12000_ratio\",\n"
        << "    \"spectral_power\"\n"
        << "  ],\n"
        << "  \"feature_notes\": {\n"
        << "    \"rms_dbfs\": \"20*log10(RMS), digital full scale; not dB SPL\",\n"
        << "    \"dominant_frequency_hz\": \"highest-magnitude FFT bin in analysis band; not necessarily F0\",\n"
        << "    \"spectral_flux\": \"L2 change of normalized magnitude spectrum\",\n"
        << "    \"band_ratio\": \"spectral power in band divided by total analysis-band power\",\n"
        << "    \"event_detection\": \"robust multi-feature heuristic; not species identification\"\n"
        << "  },\n"
        << "  \"full_spectrogram\": {\n"
        << "    \"frequency_bins\": "
        << bins << ",\n"
        << "    \"dtype\": \"float32\",\n"
        << "    \"layout\": \"row-major [frame, frequency_bin]\",\n"
        << "    \"phase_saved\": "
        << (phase_saved ? "true" : "false")
        << "\n"
        << "  },\n"
        << "  \"global_spectrogram\": {\n"
        << "    \"time_bins\": "
        << cfg.global_time_bins << ",\n"
        << "    \"frequency_bins\": "
        << cfg.global_freq_bins << ",\n"
        << "    \"values\": "
        << global_count << ",\n"
        << "    \"dtype\": \"float32\",\n"
        << "    \"aggregation\": \"mean magnitude over time/frequency cells\"\n"
        << "  },\n"
        << "  \"event_detection\": {\n"
        << "    \"threshold\": "
        << cfg.event_threshold << ",\n"
        << "    \"minimum_duration_ms\": "
        << cfg.event_min_ms << ",\n"
        << "    \"merge_gap_ms\": "
        << cfg.event_merge_ms << ",\n"
        << "    \"score_formula\": \"0.50*z_rms + 0.25*z_band_4_8k + 0.15*z_flux + 0.10*z_inverse_flatness\",\n"
        << "    \"z_scale\": \"median and 1.4826*MAD\"\n"
        << "  },\n"
        << "  \"scientific_limits\": [\n"
        << "    \"Digital amplitude is not calibrated SPL.\",\n"
        << "    \"Dominant frequency is not automatically fundamental frequency.\",\n"
        << "    \"Event labels are acoustic episodes, not species identities.\",\n"
        << "    \"The first noise_seconds are assumed to provide a noise reference for the Wiener stage.\"\n"
        << "  ]\n"
        << "}\n";
}

static void write_validation_summary(
    const std::string& path,
    uint64_t input_frames,
    uint64_t output_frames,
    uint32_t fs,
    double input_peak,
    double output_peak,
    double input_rms,
    double output_rms,
    uint64_t clipped_input,
    uint64_t clipped_output,
    size_t feature_frames,
    size_t event_count
) {
    std::ofstream out(path.c_str());

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir validation_summary.json."
        );

    out
        << "{\n"
        << "  \"input_samples\": "
        << input_frames << ",\n"
        << "  \"output_samples\": "
        << output_frames << ",\n"
        << "  \"sample_difference\": "
        << (
            input_frames >= output_frames
            ? input_frames - output_frames
            : output_frames - input_frames
        )
        << ",\n"
        << "  \"sample_rate_hz\": "
        << fs << ",\n"
        << "  \"duration_seconds\": "
        << std::setprecision(15)
        << static_cast<double>(input_frames) /
           static_cast<double>(fs)
        << ",\n"
        << "  \"input_peak\": "
        << input_peak << ",\n"
        << "  \"output_peak\": "
        << output_peak << ",\n"
        << "  \"input_rms\": "
        << input_rms << ",\n"
        << "  \"output_rms\": "
        << output_rms << ",\n"
        << "  \"input_samples_outside_unit_range\": "
        << clipped_input << ",\n"
        << "  \"output_samples_outside_unit_range\": "
        << clipped_output << ",\n"
        << "  \"feature_frames\": "
        << feature_frames << ",\n"
        << "  \"event_count\": "
        << event_count << ",\n"
        << "  \"basic_validation_pass\": "
        << (
            input_frames == output_frames &&
            clipped_output == 0 &&
            output_peak <= 1.0 + 1.0e-6
            ? "true"
            : "false"
        )
        << "\n"
        << "}\n";
}

static void write_events_csv(
    const std::string& path,
    const std::vector<Event>& events
) {
    std::ofstream out(path.c_str());

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir events.csv."
        );

    out
        << "event_id,start_time_s,end_time_s,duration_s,"
        << "peak,rms,rms_dbfs,"
        << "dominant_frequency_hz,dominant_frequency_median_hz,"
        << "dominant_frequency_slope_hz_per_s,"
        << "spectral_centroid_hz,spectral_bandwidth_hz,"
        << "spectral_rolloff_hz,spectral_flatness,spectral_flux,"
        << "band_0_500_ratio,band_500_1000_ratio,"
        << "band_1000_2000_ratio,band_2000_4000_ratio,"
        << "band_4000_8000_ratio,band_8000_12000_ratio,"
        << "spectral_power,attack_s,decay_s,detection_score,"
        << "inter_event_interval_s\n";

    for (size_t i = 0; i < events.size(); ++i) {

        const Event& e = events[i];

        double iei = std::numeric_limits<double>::quiet_NaN();

        if (i > 0)
            iei =
                e.start_seconds -
                events[i - 1].end_seconds;

        out
            << (i + 1) << ","
            << std::setprecision(12)
            << e.start_seconds << ","
            << e.end_seconds << ","
            << e.duration_seconds << ","
            << e.peak << ","
            << e.rms << ","
            << e.rms_dbfs << ","
            << e.dominant_frequency_hz << ","
            << e.dominant_frequency_median_hz << ","
            << e.dominant_frequency_slope_hz_per_s << ","
            << e.centroid_hz << ","
            << e.bandwidth_hz << ","
            << e.rolloff_hz << ","
            << e.flatness << ","
            << e.flux << ","
            << e.band_0_500 << ","
            << e.band_500_1000 << ","
            << e.band_1000_2000 << ","
            << e.band_2000_4000 << ","
            << e.band_4000_8000 << ","
            << e.band_8000_12000 << ","
            << e.spectral_power << ","
            << e.attack_seconds << ","
            << e.decay_seconds << ","
            << e.detection_score << ","
            << iei
            << "\n";
    }
}

static void write_event_features_csv(
    const std::string& path,
    const std::vector<Event>& events
) {
    std::ofstream out(path.c_str());

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir event_features.csv."
        );

    out
        << "event_id,frame_start,frame_end\n";

    for (size_t i = 0; i < events.size(); ++i) {
        out
            << (i + 1) << ","
            << events[i].first_frame << ","
            << events[i].last_frame << "\n";
    }
}

struct RunningStats {
    double sum_sq = 0.0;
    double peak = 0.0;
    uint64_t count = 0;
    uint64_t outside = 0;

    void add(double x) {
        if (!std::isfinite(x))
            return;

        sum_sq += x * x;

        peak =
            std::max(
                peak,
                std::abs(x)
            );

        if (std::abs(x) > 1.0)
            ++outside;

        ++count;
    }

    double rms() const {
        if (count == 0) return 0.0;

        return std::sqrt(
            sum_sq /
            static_cast<double>(count)
        );
    }
};

class RingOLA {
public:
    RingOLA(
        size_t nfft,
        size_t hop,
        const std::vector<float>& window,
        WavMono16Writer& writer,
        uint64_t target_samples
    )
        : nfft_(nfft),
          hop_(hop),
          window_(window),
          writer_(writer),
          target_samples_(target_samples),
          ola_(nfft, 0.0),
          norm_(nfft, 0.0),
          output_tmp_(hop, 0.0f) {}

    void add_frame(
        const std::vector<Complex>& time_domain
    ) {
        if (time_domain.size() != nfft_)
            throw std::runtime_error(
                "Tamano IFFT incorrecto."
            );

        for (size_t i = 0; i < nfft_; ++i) {

            const double y =
                time_domain[i].re *
                static_cast<double>(
                    window_[i]
                );

            const double w =
                static_cast<double>(
                    window_[i]
                );

            ola_[i] += y;
            norm_[i] += w * w;
        }

        emit_hop();
        shift_hop();
    }

    void flush() {
        while (samples_emitted_ <
               target_samples_) {

            const size_t remaining =
                static_cast<size_t>(
                    target_samples_ -
                    samples_emitted_
                );

            const size_t n =
                std::min(
                    hop_,
                    remaining
                );

            for (size_t i = 0; i < n; ++i) {

                double y = 0.0;

                if (norm_[i] > 1.0e-12)
                    y = ola_[i] /
                        norm_[i];

                output_tmp_[i] =
                    static_cast<float>(y);
            }

            writer_.write_samples(
                output_tmp_,
                n
            );

            samples_emitted_ += n;

            shift_hop();
        }
    }

private:
    size_t nfft_;
    size_t hop_;

    const std::vector<float>& window_;

    WavMono16Writer& writer_;
    uint64_t target_samples_;

    std::vector<double> ola_;
    std::vector<double> norm_;
    std::vector<float> output_tmp_;

    uint64_t samples_emitted_ = 0;

    void emit_hop() {

        if (samples_emitted_ >= target_samples_)
            return;

        const size_t remaining =
            static_cast<size_t>(
                target_samples_ -
                samples_emitted_
            );

        const size_t n =
            std::min(
                hop_,
                remaining
            );

        for (size_t i = 0; i < n; ++i) {

            double y = 0.0;

            if (norm_[i] > 1.0e-12)
                y = ola_[i] /
                    norm_[i];

            output_tmp_[i] =
                static_cast<float>(y);
        }

        writer_.write_samples(
            output_tmp_,
            n
        );

        samples_emitted_ += n;
    }

    void shift_hop() {

        if (hop_ >= nfft_) {
            std::fill(
                ola_.begin(),
                ola_.end(),
                0.0
            );

            std::fill(
                norm_.begin(),
                norm_.end(),
                0.0
            );

            return;
        }

        const size_t remain =
            nfft_ - hop_;

        for (size_t i = 0; i < remain; ++i) {
            ola_[i] =
                ola_[i + hop_];

            norm_[i] =
                norm_[i + hop_];
        }

        for (size_t i = remain; i < nfft_; ++i) {
            ola_[i] = 0.0;
            norm_[i] = 0.0;
        }
    }
};


static void process_frames_pass(
    const std::string& input_path,
    const Config& cfg,
    const WavInfo& wi,
    const std::vector<double>& noise_profile,
    double dc_offset,
    std::ofstream& clean_mag_out,
    std::ofstream& clean_phase_out,
    std::ofstream& features_out,
    WavMono16Writer& clean_writer,
    std::vector<FeatureFrame>& features,
    std::vector<float>& global_sums,
    std::vector<uint32_t>& global_counts,
    RunningStats& input_stats,
    uint64_t& frames_processed
) {
    WavReader reader(input_path);

    const uint64_t total_frames = reader.total_frames();

    const size_t block_frames =
        std::max<size_t>(
            cfg.nfft,
            static_cast<size_t>(
                std::llround(
                    cfg.block_seconds *
                    static_cast<double>(wi.sample_rate)
                )
            )
        );

    const std::vector<float> window =
        hann_window(cfg.nfft);

    const size_t bins = cfg.nfft / 2 + 1;

    std::vector<float> block;
    std::vector<float> carry;
    std::vector<float> frame(cfg.nfft, 0.0f);
    std::vector<float> clean_frame(cfg.nfft, 0.0f);

    std::vector<Complex> spectrum(cfg.nfft);
    std::vector<Complex> time_domain(cfg.nfft);

    std::vector<double> magnitude(bins, 0.0);
    std::vector<double> previous_norm(bins, 0.0);
    std::vector<double> current_norm(bins, 0.0);

    RingOLA ola(
        cfg.nfft,
        cfg.hop,
        window,
        clean_writer,
        total_frames
    );

    const double df =
        static_cast<double>(wi.sample_rate) /
        static_cast<double>(cfg.nfft);

    const size_t freq_max_bin =
        std::min(
            bins - 1,
            static_cast<size_t>(
                std::floor(cfg.fmax_hz / df)
            )
        );

    bool have_previous = false;
    uint64_t frame_start = 0;

    auto process_one_frame = [&](const std::vector<float>& input_frame,
                                 uint64_t start_index) {
        if (input_frame.size() != cfg.nfft)
            throw std::runtime_error("Tamano de frame incorrecto.");

        // ---------------------------------------------------------
        // Analysis STFT of the original frame.
        // ---------------------------------------------------------
        for (size_t i = 0; i < cfg.nfft; ++i) {
            spectrum[i] = {
                static_cast<double>(input_frame[i]) *
                    static_cast<double>(window[i]),
                0.0
            };
        }

        fft(spectrum, false);

        for (size_t k = 0; k < bins; ++k) {
            magnitude[k] = std::hypot(
                spectrum[k].re,
                spectrum[k].im
            );
        }

        // ---------------------------------------------------------
        // Wiener mask.
        // ---------------------------------------------------------
        for (size_t k = 0; k < cfg.nfft; ++k) {
            const size_t positive_k =
                (k <= cfg.nfft / 2)
                ? k
                : cfg.nfft - k;

            const double noise =
                noise_profile[
                    std::min(
                        positive_k,
                        noise_profile.size() - 1
                    )
                ];

            const double denominator =
                cfg.noise_factor *
                std::max(noise, 1.0e-12);

            const double snr =
                magnitude[positive_k] /
                denominator;

            double mask =
                (snr * snr) /
                (1.0 + snr * snr);

            if (mask < cfg.gain_floor)
                mask = cfg.gain_floor;

            spectrum[k].re *= mask;
            spectrum[k].im *= mask;
        }

        // ---------------------------------------------------------
        // Save CLEAN one-sided spectrum.
        // ---------------------------------------------------------
        for (size_t k = 0; k < bins; ++k) {
            const float clean_mag =
                static_cast<float>(
                    std::hypot(
                        spectrum[k].re,
                        spectrum[k].im
                    )
                );

            clean_mag_out.write(
                reinterpret_cast<const char*>(&clean_mag),
                sizeof(float)
            );

            if (cfg.save_phase) {
                const float phase =
                    static_cast<float>(
                        std::atan2(
                            spectrum[k].im,
                            spectrum[k].re
                        )
                    );

                clean_phase_out.write(
                    reinterpret_cast<const char*>(&phase),
                    sizeof(float)
                );
            }
        }

        // ---------------------------------------------------------
        // Hermitian reconstruction + IFFT.
        // ---------------------------------------------------------
        for (size_t k = 1; k < cfg.nfft / 2; ++k) {
            spectrum[cfg.nfft - k] =
                c_conj(spectrum[k]);
        }

        spectrum[cfg.nfft / 2].im = 0.0;

        fft(spectrum, true);

        for (size_t i = 0; i < cfg.nfft; ++i) {
            time_domain[i] = spectrum[i];
            clean_frame[i] =
                static_cast<float>(spectrum[i].re);
        }

        // ---------------------------------------------------------
        // Streaming overlap-add reconstruction.
        // ---------------------------------------------------------
        ola.add_frame(time_domain);

        // ---------------------------------------------------------
        // Features are calculated from the CLEAN frame, not the raw
        // frame. This keeps the feature stream consistent with what
        // the user will visualize/analyze later.
        // ---------------------------------------------------------
        FeatureFrame f =
            calculate_features(
                clean_frame,
                window,
                spectrum,
                magnitude,
                previous_norm,
                have_previous,
                current_norm,
                wi.sample_rate,
                cfg.nfft,
                start_index,
                cfg.fmin_hz,
                cfg.fmax_hz,
                cfg.rolloff_fraction
            );

        features.push_back(f);
        write_feature(features_out, f);

        // ---------------------------------------------------------
        // Multiresolution global spectrogram.
        // ---------------------------------------------------------
        size_t time_bin = 0;

        if (total_frames > 0) {
            const long double fraction =
                static_cast<long double>(start_index) /
                static_cast<long double>(total_frames);

            time_bin = static_cast<size_t>(
                fraction *
                static_cast<long double>(cfg.global_time_bins)
            );

            if (time_bin >= cfg.global_time_bins)
                time_bin = cfg.global_time_bins - 1;
        }

        for (size_t gk = 0;
             gk < cfg.global_freq_bins;
             ++gk) {

            const double freq_fraction =
                cfg.global_freq_bins > 1
                ? static_cast<double>(gk) /
                  static_cast<double>(cfg.global_freq_bins - 1)
                : 0.0;

            const double target_freq =
                freq_fraction * cfg.fmax_hz;

            size_t source_k =
                static_cast<size_t>(
                    std::llround(
                        target_freq / df
                    )
                );

            source_k =
                std::min(
                    source_k,
                    freq_max_bin
                );

            const size_t cell =
                time_bin * cfg.global_freq_bins +
                gk;

            global_sums[cell] +=
                static_cast<float>(magnitude[source_k]);

            global_counts[cell] += 1;
        }

        previous_norm.swap(current_norm);
        current_norm.assign(bins, 0.0);
        have_previous = true;
    };

    while (true) {
        const size_t got =
            reader.read_mono(block, block_frames);

        if (got == 0)
            break;

        for (size_t i = 0; i < got; ++i) {
            block[i] = static_cast<float>(
                static_cast<double>(block[i]) - dc_offset
            );
            input_stats.add(
                static_cast<double>(block[i])
            );
        }

        std::vector<float> work;
        work.reserve(carry.size() + block.size());
        work.insert(work.end(), carry.begin(), carry.end());
        work.insert(work.end(), block.begin(), block.end());

        size_t pos = 0;

        while (pos + cfg.nfft <= work.size()) {
            std::copy(
                work.begin() +
                    static_cast<std::ptrdiff_t>(pos),
                work.begin() +
                    static_cast<std::ptrdiff_t>(pos + cfg.nfft),
                frame.begin()
            );

            process_one_frame(frame, frame_start);

            ++frames_processed;
            frame_start += cfg.hop;
            pos += cfg.hop;
        }

        carry.assign(
            work.begin() +
                static_cast<std::ptrdiff_t>(pos),
            work.end()
        );

        block.clear();

        const double percent =
            100.0 *
            static_cast<double>(
                std::min<uint64_t>(
                    total_frames,
                    frame_start
                )
            ) /
            std::max<double>(
                1.0,
                static_cast<double>(total_frames)
            );

        std::cerr
            << "\rProcesado: "
            << std::fixed
            << std::setprecision(1)
            << percent
            << "%    frames STFT: "
            << frames_processed
            << "   ";
    }

    // -------------------------------------------------------------
    // Final zero-padded frame. It is necessary whenever the recording
    // does not end exactly on an STFT frame boundary. The synthesis
    // stage then flushes only the original number of samples.
    // -------------------------------------------------------------
    if (!carry.empty() &&
        frame_start < total_frames) {

        std::fill(
            frame.begin(),
            frame.end(),
            0.0f
        );

        const size_t n =
            std::min(
                carry.size(),
                cfg.nfft
            );

        std::copy(
            carry.begin(),
            carry.begin() +
                static_cast<std::ptrdiff_t>(n),
            frame.begin()
        );

        process_one_frame(
            frame,
            frame_start
        );

        ++frames_processed;
    }

    ola.flush();
}


static double estimate_dc_offset(
    const std::string& input_path
) {
    WavReader reader(input_path);

    const size_t block =
        std::max<size_t>(
            4096,
            static_cast<size_t>(
                reader.info().sample_rate * 10ULL
            )
        );

    long double sum = 0.0L;
    uint64_t count = 0;

    while (true) {
        std::vector<float> x;
        const size_t got = reader.read_mono(x, block);

        if (got == 0)
            break;

        for (size_t i = 0; i < got; ++i) {
            sum += static_cast<long double>(x[i]);
            ++count;
        }
    }

    if (count == 0)
        return 0.0;

    return static_cast<double>(sum / static_cast<long double>(count));
}

static std::vector<double> estimate_noise_profile(
    const std::string& input_path,
    const WavInfo& wi,
    size_t nfft,
    size_t hop,
    double noise_seconds,
    double& actual_seconds,
    size_t& noise_frames,
    double dc_offset
) {
    WavReader reader(input_path);

    const uint64_t desired =
        std::min<uint64_t>(
            reader.total_frames(),
            static_cast<uint64_t>(
                std::llround(
                    noise_seconds *
                    static_cast<double>(
                        wi.sample_rate
                    )
                )
            )
        );

    actual_seconds =
        static_cast<double>(desired) /
        static_cast<double>(
            wi.sample_rate
        );

    const size_t bins =
        nfft / 2 + 1;

    const std::vector<float> window =
        hann_window(nfft);

    std::vector<std::vector<float>> mags;

    std::vector<float> samples;

    const size_t read_block =
        std::max<size_t>(
            nfft,
            20 * static_cast<size_t>(
                wi.sample_rate
            )
        );

    while (samples.size() < desired) {

        std::vector<float> block;

        const size_t need =
            std::min<size_t>(
                read_block,
                static_cast<size_t>(
                    desired -
                    static_cast<uint64_t>(
                        samples.size()
                    )
                )
            );

        if (need == 0)
            break;

        const size_t got =
            reader.read_mono(
                block,
                need
            );

        if (got == 0)
            break;

        for (size_t i = 0; i < block.size(); ++i) {
            block[i] = static_cast<float>(
                static_cast<double>(block[i]) - dc_offset
            );
        }

        samples.insert(
            samples.end(),
            block.begin(),
            block.end()
        );
    }

    std::vector<Complex> spectrum(nfft);
    std::vector<double> magnitude(bins);

    size_t pos = 0;

    while (pos + nfft <= samples.size()) {

        for (size_t i = 0;
             i < nfft;
             ++i) {

            spectrum[i] = {
                static_cast<double>(
                    samples[pos + i]
                ) *
                static_cast<double>(
                    window[i]
                ),
                0.0
            };
        }

        fft(
            spectrum,
            false
        );

        std::vector<float> row(bins);

        for (size_t k = 0;
             k < bins;
             ++k) {

            row[k] =
                static_cast<float>(
                    std::hypot(
                        spectrum[k].re,
                        spectrum[k].im
                    )
                );
        }

        mags.push_back(
            std::move(row)
        );

        pos += hop;
    }

    if (mags.empty())
        throw std::runtime_error(
            "No se pudieron obtener frames para estimar ruido."
        );

    noise_frames =
        mags.size();

    std::vector<double> profile(bins);

    std::vector<double> column(
        mags.size()
    );

    for (size_t k = 0;
         k < bins;
         ++k) {

        for (size_t t = 0;
             t < mags.size();
             ++t) {

            column[t] =
                static_cast<double>(
                    mags[t][k]
                );
        }

        profile[k] =
            median_copy(column);
    }

    return profile;
}

static void write_float_vector(
    const std::string& path,
    const std::vector<double>& data
) {
    std::ofstream out(
        path.c_str(),
        std::ios::binary
    );

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir perfil de ruido."
        );

    for (double x : data) {

        const float v =
            static_cast<float>(x);

        out.write(
            reinterpret_cast<const char*>(&v),
            sizeof(float)
        );
    }
}

static void write_global_spectrogram(
    const std::string& path,
    const std::vector<float>& sums,
    const std::vector<uint32_t>& counts
) {
    std::ofstream out(
        path.c_str(),
        std::ios::binary
    );

    if (!out)
        throw std::runtime_error(
            "No se pudo escribir spectrogram_global.f32."
        );

    for (size_t i = 0;
         i < sums.size();
         ++i) {

        const float v =
            counts[i] > 0
            ? sums[i] /
              static_cast<float>(
                  counts[i]
              )
            : 0.0f;

        out.write(
            reinterpret_cast<const char*>(&v),
            sizeof(float)
        );
    }
}

static void read_wav_stats(
    const std::string& path,
    RunningStats& stats
) {
    WavReader reader(path);

    const size_t block =
        std::max<size_t>(
            4096,
            static_cast<size_t>(
                reader.info().sample_rate *
                10
            )
        );

    while (true) {

        std::vector<float> x;

        const size_t got =
            reader.read_mono(
                x,
                block
            );

        if (got == 0)
            break;

        for (size_t i = 0;
             i < got;
             ++i)
            stats.add(
                static_cast<double>(
                    x[i]
                )
            );
    }
}

} // namespace cicada

int main(int argc, char** argv) {
    using namespace cicada;

    try {
        Config cfg;

        const std::string exe_dir =
            executable_directory();

        const std::string raw =
            exe_dir +
            "\\Cicada.wav";

        // The FINAL engine performs denoising itself, so the default
        // input is always the original recording. To analyze another
        // WAV (including a previously cleaned one), use --input.
        cfg.input = raw;

        cfg.output =
            exe_dir +
            "\\results";

        for (int i = 1;
             i < argc;
             ++i) {

            const std::string arg(
                argv[i]
            );

            if (arg == "--help") {
                usage();
                return 0;
            }

            if (arg == "--no-phase") {
                cfg.save_phase = false;
                continue;
            }

            if (arg == "--input" &&
                i + 1 < argc) {

                cfg.input = argv[++i];
                continue;
            }

            if (arg == "--output" &&
                i + 1 < argc) {

                cfg.output = argv[++i];
                continue;
            }

            if (arg == "--nfft" &&
                i + 1 < argc) {

                if (!parse_size(
                        argv[++i],
                        cfg.nfft))
                    throw std::runtime_error(
                        "NFFT invalido."
                    );

                continue;
            }

            if (arg == "--hop" &&
                i + 1 < argc) {

                if (!parse_size(
                        argv[++i],
                        cfg.hop))
                    throw std::runtime_error(
                        "HOP invalido."
                    );

                continue;
            }

            if (arg == "--block" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.block_seconds))
                    throw std::runtime_error(
                        "BLOCK invalido."
                    );

                continue;
            }

            if (arg == "--noise" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.noise_seconds))
                    throw std::runtime_error(
                        "NOISE invalido."
                    );

                continue;
            }

            if (arg == "--noise-factor" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.noise_factor))
                    throw std::runtime_error(
                        "NOISE FACTOR invalido."
                    );

                continue;
            }

            if (arg == "--floor" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.gain_floor))
                    throw std::runtime_error(
                        "FLOOR invalido."
                    );

                continue;
            }

            if (arg == "--fmin" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.fmin_hz))
                    throw std::runtime_error(
                        "FMIN invalido."
                    );

                continue;
            }

            if (arg == "--fmax" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.fmax_hz))
                    throw std::runtime_error(
                        "FMAX invalido."
                    );

                continue;
            }

            if (arg == "--rolloff" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.rolloff_fraction))
                    throw std::runtime_error(
                        "ROLLOFF invalido."
                    );

                continue;
            }

            if (arg == "--global-time-bins" &&
                i + 1 < argc) {

                if (!parse_size(
                        argv[++i],
                        cfg.global_time_bins))
                    throw std::runtime_error(
                        "GLOBAL TIME BINS invalido."
                    );

                continue;
            }

            if (arg == "--global-freq-bins" &&
                i + 1 < argc) {

                if (!parse_size(
                        argv[++i],
                        cfg.global_freq_bins))
                    throw std::runtime_error(
                        "GLOBAL FREQ BINS invalido."
                    );

                continue;
            }

            if (arg == "--event-threshold" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.event_threshold))
                    throw std::runtime_error(
                        "EVENT THRESHOLD invalido."
                    );

                continue;
            }

            if (arg == "--event-min-ms" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.event_min_ms))
                    throw std::runtime_error(
                        "EVENT MIN invalido."
                    );

                continue;
            }

            if (arg == "--event-merge-ms" &&
                i + 1 < argc) {

                if (!parse_double(
                        argv[++i],
                        cfg.event_merge_ms))
                    throw std::runtime_error(
                        "EVENT MERGE invalido."
                    );

                continue;
            }

            throw std::runtime_error(
                "Opcion desconocida: " + arg
            );
        }

        if ((cfg.nfft & (cfg.nfft - 1)) != 0)
            throw std::runtime_error(
                "NFFT debe ser potencia de 2."
            );

        if (cfg.hop == 0 ||
            cfg.hop > cfg.nfft)
            throw std::runtime_error(
                "HOP debe estar entre 1 y NFFT."
            );

        if (cfg.block_seconds <= 0.0)
            throw std::runtime_error(
                "BLOCK debe ser positivo."
            );

        if (cfg.noise_seconds <= 0.0)
            throw std::runtime_error(
                "NOISE debe ser positivo."
            );

        if (cfg.noise_factor <= 0.0)
            throw std::runtime_error(
                "NOISE FACTOR debe ser positivo."
            );

        if (cfg.gain_floor < 0.0 ||
            cfg.gain_floor > 1.0)
            throw std::runtime_error(
                "FLOOR debe estar entre 0 y 1."
            );

        if (cfg.fmin_hz < 0.0)
            throw std::runtime_error(
                "FMIN no puede ser negativo."
            );

        if (cfg.fmax_hz <= cfg.fmin_hz)
            throw std::runtime_error(
                "FMAX debe ser mayor que FMIN."
            );

        if (cfg.rolloff_fraction <= 0.0 ||
            cfg.rolloff_fraction > 1.0)
            throw std::runtime_error(
                "ROLLOFF debe estar entre 0 y 1."
            );

        if (cfg.global_time_bins == 0 ||
            cfg.global_freq_bins == 0)
            throw std::runtime_error(
                "Bins globales invalidos."
            );

        if (cfg.event_threshold <= 0.0)
            throw std::runtime_error(
                "EVENT THRESHOLD debe ser positivo."
            );

        if (cfg.event_min_ms < 0.0 ||
            cfg.event_merge_ms < 0.0)
            throw std::runtime_error(
                "Duraciones de eventos invalidas."
            );

        ensure_directory(cfg.output);

        WavReader initial_reader(
            cfg.input
        );

        const WavInfo wi =
            initial_reader.info();

        const uint64_t total_frames =
            initial_reader.total_frames();

        const double nyquist =
            static_cast<double>(
                wi.sample_rate
            ) / 2.0;

        const double dc_offset =
            estimate_dc_offset(cfg.input);

        cfg.fmax_hz =
            std::min(
                cfg.fmax_hz,
                nyquist
            );

        const size_t bins =
            cfg.nfft / 2 + 1;

        std::cout
            << "============================================================\n"
            << "                 CICADA ANALYZER - FINAL                    \n"
            << "============================================================\n"
            << "Entrada   : " << cfg.input << "\n"
            << "Salida    : " << cfg.output << "\n"
            << "------------------------------------------------------------\n"
            << "Fs        : " << wi.sample_rate << " Hz\n"
            << "Canales   : " << wi.channels << "\n"
            << "Bits      : " << wi.bits_per_sample << "\n"
            << "DC offset : " << std::setprecision(10) << dc_offset << "\n"
            << "Duracion  : "
            << std::fixed
            << std::setprecision(3)
            << static_cast<double>(
                total_frames
            ) /
            static_cast<double>(
                wi.sample_rate
            )
            << " s\n"
            << "------------------------------------------------------------\n"
            << "STFT      : NFFT=" << cfg.nfft
            << ", HOP=" << cfg.hop
            << ", ventana=Hann simetrica\n"
            << "Bloque    : " << cfg.block_seconds << " s\n"
            << "Ruido     : " << cfg.noise_seconds << " s\n"
            << "Factor    : " << cfg.noise_factor << "\n"
            << "Piso      : " << cfg.gain_floor << "\n"
            << "Fmin/Fmax : "
            << cfg.fmin_hz << " / "
            << cfg.fmax_hz << " Hz\n"
            << "------------------------------------------------------------\n"
            << "Detector  : multicaracteristica robusto\n"
            << "Umbral    : " << cfg.event_threshold << "\n"
            << "Min. evento: " << cfg.event_min_ms << " ms\n"
            << "Merge     : " << cfg.event_merge_ms << " ms\n"
            << "Fase STFT : "
            << (cfg.save_phase ? "guardada" : "no guardada")
            << "\n"
            << "============================================================\n\n";

        // -----------------------------------------------------
        // PASS 1: noise profile
        // -----------------------------------------------------
        std::cout
            << "1/4  Estimando perfil de ruido...\n";

        double actual_noise_seconds = 0.0;
        size_t noise_frames = 0;

        const std::vector<double> noise_profile =
            estimate_noise_profile(
                cfg.input,
                wi,
                cfg.nfft,
                cfg.hop,
                cfg.noise_seconds,
                actual_noise_seconds,
                noise_frames,
                dc_offset
            );

        const std::string noise_path =
            cfg.output +
            "\\noise_profile.f32";

        write_float_vector(
            noise_path,
            noise_profile
        );

        std::cout
            << "     Frames ruido : "
            << noise_frames
            << "\n"
            << "     Duracion real: "
            << actual_noise_seconds
            << " s\n";

        // -----------------------------------------------------
        // Output files
        // -----------------------------------------------------
        const std::string clean_wav_path =
            cfg.output +
            "\\audio_clean.wav";

        const std::string mag_path =
            cfg.output +
            "\\spectrogram_clean_magnitude.f32";

        const std::string phase_path =
            cfg.output +
            "\\spectrogram_clean_phase.f32";

        const std::string features_path =
            cfg.output +
            "\\features.f32";

        const std::string global_path =
            cfg.output +
            "\\spectrogram_global.f32";

        const std::string events_path =
            cfg.output +
            "\\events.csv";

        const std::string event_features_path =
            cfg.output +
            "\\event_features.csv";

        const std::string metadata_path =
            cfg.output +
            "\\metadata_final.json";

        const std::string validation_path =
            cfg.output +
            "\\validation_summary.json";

        std::ofstream mag_out(
            mag_path.c_str(),
            std::ios::binary
        );

        if (!mag_out)
            throw std::runtime_error(
                "No se pudo crear spectrogram_clean_magnitude.f32."
            );

        std::ofstream phase_out;

        if (cfg.save_phase) {
            phase_out.open(
                phase_path.c_str(),
                std::ios::binary
            );

            if (!phase_out)
                throw std::runtime_error(
                    "No se pudo crear spectrogram_clean_phase.f32."
                );
        }

        std::ofstream features_out(
            features_path.c_str(),
            std::ios::binary
        );

        if (!features_out)
            throw std::runtime_error(
                "No se pudo crear features.f32."
            );

        WavMono16Writer clean_writer(
            clean_wav_path,
            wi.sample_rate,
            total_frames
        );

        const size_t global_cells =
            cfg.global_time_bins *
            cfg.global_freq_bins;

        std::vector<float> global_sums(
            global_cells,
            0.0f
        );

        std::vector<uint32_t> global_counts(
            global_cells,
            0
        );

        std::vector<FeatureFrame> features;
        features.reserve(
            static_cast<size_t>(
                total_frames / cfg.hop
            ) + 2
        );

        RunningStats input_stats;
        uint64_t frames_processed = 0;

        // -----------------------------------------------------
        // PASS 2: denoise + reconstruction + features
        // -----------------------------------------------------
        std::cout
            << "\n2/4  Procesando grabacion completa...\n";

        process_frames_pass(
            cfg.input,
            cfg,
            wi,
            noise_profile,
            dc_offset,
            mag_out,
            phase_out,
            features_out,
            clean_writer,
            features,
            global_sums,
            global_counts,
            input_stats,
            frames_processed
        );

        mag_out.close();

        if (cfg.save_phase)
            phase_out.close();

        features_out.close();

        clean_writer.close();

        // Recalculate exact stats from input + output files.
        RunningStats input_exact;
        RunningStats output_exact;

        read_wav_stats(
            cfg.input,
            input_exact
        );

        read_wav_stats(
            clean_wav_path,
            output_exact
        );

        // -----------------------------------------------------
        // Global spectrogram
        // -----------------------------------------------------
        write_global_spectrogram(
            global_path,
            global_sums,
            global_counts
        );

        // -----------------------------------------------------
        // Phase 4 event detection
        // -----------------------------------------------------
        std::cout
            << "\n3/4  Detectando eventos acusticos...\n";

        const std::vector<Event> events =
            detect_events(
                features,
                wi.sample_rate,
                cfg.hop,
                cfg.event_threshold,
                cfg.event_min_ms,
                cfg.event_merge_ms
            );

        write_events_csv(
            events_path,
            events
        );

        write_event_features_csv(
            event_features_path,
            events
        );

        // -----------------------------------------------------
        // Metadata + validation
        // -----------------------------------------------------
        write_metadata(
            metadata_path,
            wi,
            total_frames,
            features.size(),
            cfg.nfft,
            cfg.hop,
            cfg,
            cfg.save_phase,
            dc_offset
        );

        write_validation_summary(
            validation_path,
            total_frames,
            output_exact.count,
            wi.sample_rate,
            input_exact.peak,
            output_exact.peak,
            input_exact.rms(),
            output_exact.rms(),
            input_exact.outside,
            output_exact.outside,
            features.size(),
            events.size()
        );

        std::cout
            << "\n4/4  Finalizando resultados...\n"
            << "\n============================================================\n"
            << "                 CICADA FINAL - OK                          \n"
            << "============================================================\n"
            << "Frecuencia de muestreo : "
            << wi.sample_rate << " Hz\n"
            << "Duracion                : "
            << static_cast<double>(
                total_frames
            ) /
            static_cast<double>(
                wi.sample_rate
            )
            << " s\n"
            << "Frames STFT             : "
            << features.size()
            << "\n"
            << "Bins de frecuencia      : "
            << bins
            << "\n"
            << "Resolucion frecuencia   : "
            << static_cast<double>(
                wi.sample_rate
            ) /
            static_cast<double>(
                cfg.nfft
            )
            << " Hz\n"
            << "Paso temporal            : "
            << static_cast<double>(
                cfg.hop
            ) /
            static_cast<double>(
                wi.sample_rate
            )
            << " s\n"
            << "------------------------------------------------------------\n"
            << "RMS entrada              : "
            << input_exact.rms()
            << "\n"
            << "RMS salida               : "
            << output_exact.rms()
            << "\n"
            << "Pico entrada             : "
            << input_exact.peak
            << "\n"
            << "Pico salida              : "
            << output_exact.peak
            << "\n"
            << "Samples entrada > |1|    : "
            << input_exact.outside
            << "\n"
            << "Samples salida > |1|     : "
            << output_exact.outside
            << "\n"
            << "Muestras salida          : "
            << output_exact.count
            << "\n"
            << "Eventos detectados       : "
            << events.size()
            << "\n"
            << "------------------------------------------------------------\n"
            << "Audio limpio             : "
            << clean_wav_path
            << "\n"
            << "Perfil de ruido          : "
            << noise_path
            << "\n"
            << "STFT magnitud limpia     : "
            << mag_path
            << "\n";

        if (cfg.save_phase)
            std::cout
                << "STFT fase limpia         : "
                << phase_path
                << "\n";

        std::cout
            << "Features                 : "
            << features_path
            << "\n"
            << "Espectrograma global     : "
            << global_path
            << "\n"
            << "Eventos                  : "
            << events_path
            << "\n"
            << "Metadatos                : "
            << metadata_path
            << "\n"
            << "Validacion               : "
            << validation_path
            << "\n"
            << "============================================================\n";

        if (output_exact.count != total_frames)
            throw std::runtime_error(
                "VALIDACION FALLIDA: la salida no conserva "
                "el numero de muestras original."
            );

        if (output_exact.outside != 0 ||
            output_exact.peak > 1.0 + 1.0e-6)
            throw std::runtime_error(
                "VALIDACION FALLIDA: clipping o valores "
                "fuera del rango digital."
            );

        return 0;
    }
    catch (const std::exception& e) {

        std::cerr
            << "\nERROR: "
            << e.what()
            << "\n";

        return 1;
    }
}
