#ifndef CORE_REAL_COLOR_H_
#define CORE_REAL_COLOR_H_

#include <ostream>

enum class RealColor {
    RC_EMPTY,
    RC_OJAMA,
    RC_RED,
    RC_BLUE,
    RC_YELLOW,
    RC_GREEN,
    RC_PURPLE,
};

const int NUM_REAL_COLORS = 7;

constexpr int ordinal(RealColor rc) { return static_cast<int>(rc); }
constexpr RealColor intToRealColor(int x) { return static_cast<RealColor>(x); }
RealColor toRealColor(char c);
const char* toString(RealColor);
std::ostream& operator<<(std::ostream&, RealColor);

inline bool isNormalColor(RealColor rc)
{
    switch (rc) {
    case RealColor::RC_RED:
    case RealColor::RC_BLUE:
    case RealColor::RC_YELLOW:
    case RealColor::RC_GREEN:
    case RealColor::RC_PURPLE:
        return true;
    default:
        return false;
    }
}

#endif
