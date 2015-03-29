#ifndef CORE_RENSA_TRACK_RESULT_H_
#define CORE_RENSA_TRACK_RESULT_H_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "core/field_constant.h"
#include "core/position.h"

class RensaNonTrackResult {
};

// RensaChainTrackResult represents in what-th chain puyo is erased.
class RensaChainTrackResult {
public:
    RensaChainTrackResult();
    explicit RensaChainTrackResult(const std::string&);

    RensaChainTrackResult& operator=(const RensaChainTrackResult& result);

    // Nth Rensa where (x, y) is erased. 0 if not erased.
    int erasedAt(int x, int y) const { return erasedAt_[x][y]; }
    void setErasedAt(int x, int y, int nthChain) { erasedAt_[x][y] = nthChain; }

    std::string toString() const;

private:
    std::uint8_t erasedAt_[FieldConstant::MAP_WIDTH][FieldConstant::MAP_HEIGHT];
};

// RensaCoefResult represents
//  - the number of erased puyo in each chain
//  - the coefficient in each chain.
class RensaCoefResult {
public:
    RensaCoefResult() : numErased_{}, coef_{} {}

    void setCoef(int nth, int numErased, int coef)
    {
        numErased_[nth] = numErased;
        coef_[nth] = coef;
    }

    int numErased(int nth) const { return numErased_[nth]; }
    int coef(int nth) const { return coef_[nth]; }

private:
    int numErased_[20]; // numErased does not contain ojama puyos.
    int coef_[20];
};

// This tracks puyo position at "n-1"-th chain, where the puyo vanishes at n-th chain.
class RensaVanishingPositionResult {
public:
    RensaVanishingPositionResult() {
        basePuyosErasedAt_.reserve(19);
        fallingPuyosErasedAt_.reserve(19);
    }

    int size() const {
        return fallingPuyosErasedAt_.size();
    }

    // Gets the reference of position set, where the puyos vanish after falling at n-th chain.
    // The return value is valid while this result instance is valid.
    const std::vector<Position>& getReferenceFallingPuyosAt(int nthChain) const {
        return fallingPuyosErasedAt_[nthChain - 1];
    }

    // Gets the reference of position set, where the puyos vanish without falling at n-th chain.
    // The return value is valid while this result instance is valid.
    const std::vector<Position>& getReferenceBasePuyosAt(int nthChain) const {
        return basePuyosErasedAt_[nthChain - 1];
    }

    std::array<float, 2> getWeightedCenterAfterFall(int nthChain) const;

    void setFallingPuyo(int x, int yBeforeFall, int yAfterFall, int nthChain);
    void setBasePuyo(int x, int y, int nthChain);

private:
    std::vector<std::vector<Position>> basePuyosErasedAt_;
    std::vector<std::vector<Position>> fallingPuyosErasedAt_;
    std::vector<std::vector<int>> yOfFalledPuyosErasedAt_;
    void maybeResize(int nthChain);
};

#endif