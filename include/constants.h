#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants {
    // Gravitational constant (m^3 kg^-1 s^-2)
    constexpr float G = 6.67430e-11f;
    
    // Default simulation parameters
    constexpr float DEFAULT_SOFTENING = 1e-9f;
    constexpr float DEFAULT_DT = 0.01f;
}

#endif // CONSTANTS_H
