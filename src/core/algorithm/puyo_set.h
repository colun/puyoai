#ifndef CORE_ALGORITHM_PUYO_SET_H_
#define CORE_ALGORITHM_PUYO_SET_H_

#include <stdio.h>
#include <string>

#include "core/puyo_color.h"

class ColumnPuyoList;

class PuyoSet {
public:
    PuyoSet() :
        m_red(0),
        m_blue(0),
        m_yellow(0),
        m_green(0)
    {
    }

    PuyoSet(int red, int blue, int yellow, int green)
    {
        DCHECK(0 <= red && red < 16);
        DCHECK(0 <= blue && blue < 16);
        DCHECK(0 <= yellow && yellow < 16);
        DCHECK(0 <= green && green < 16);

        m_red = red;
        m_blue = blue;
        m_yellow = yellow;
        m_green = green;
    }

    std::string toString() const
    {
        char buf[80];
        sprintf(buf, "(%d, %d, %d, %d)", red(), blue(), yellow(), green());

        return buf;
    }

    int count() const { return m_red + m_blue + m_yellow + m_green; }

    int red() const { return m_red; }
    int blue() const { return m_blue; }
    int yellow() const { return m_yellow; }
    int green() const { return m_green; }

    void add(const ColumnPuyoList&);

    void add(PuyoSet set)
    {
        m_red += set.m_red;
        m_blue += set.m_blue;
        m_yellow += set.m_yellow;
        m_green += set.m_green;
    }

    void sub(PuyoSet set)
    {
        m_red = m_red < set.m_red ? 0 : m_red - set.m_red;
        m_blue = m_blue < set.m_blue ? 0 : m_blue - set.m_blue;
        m_yellow = m_yellow < set.m_yellow ? 0 : m_yellow - set.m_yellow;
        m_green = m_green < set.m_green ? 0 : m_green - set.m_green;
    }

    void add(PuyoColor c, int n = 1)
    {
        switch (c) {
        case PuyoColor::RED:
            DCHECK(m_red + n < 16);
            m_red += n;
            break;
        case PuyoColor::BLUE:
            DCHECK(m_blue + n < 16);
            m_blue += n;
            break;
        case PuyoColor::YELLOW:
            DCHECK(m_yellow + n < 16);
            m_yellow += n;
            break;
        case PuyoColor::GREEN:
            DCHECK(m_green + n < 16);
            m_green += n;
            break;
        default:
            DCHECK(false);
        }
    }

    friend bool operator<(const PuyoSet& lhs, const PuyoSet& rhs)
    {
        return lhs.toInt() < rhs.toInt();
    }

    friend bool operator==(const PuyoSet& lhs, const PuyoSet& rhs)
    {
        return lhs.toInt() == rhs.toInt();
    }

    friend bool operator!=(const PuyoSet& lhs, const PuyoSet& rhs)
    {
        return lhs.toInt() != rhs.toInt();
    }

private:
    int toInt() const
    {
        return m_red | (m_blue << 4) | (m_yellow << 8) | (m_green << 12);
    }

    int m_red;
    int m_blue;
    int m_yellow;
    int m_green;
};

#endif
