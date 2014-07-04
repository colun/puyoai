#include "rensa_detector.h"

#include <iostream>
#include <set>
#include <vector>

#include "core/algorithm/column_puyo_list.h"
#include "core/algorithm/plan.h"
#include "core/algorithm/rensa_info.h"
#include "core/decision.h"
#include "core/field/core_field.h"
#include "core/kumipuyo.h"

using namespace std;

namespace {

typedef std::function<void (CoreField*, const CoreField&, int, PuyoColor, int)> SimulationCallback;

};

static inline void tryDropFire(const CoreField& originalField, SimulationCallback callback)
{
    bool visited[CoreField::MAP_WIDTH][NUM_PUYO_COLORS] {};

    for (int x = 1; x <= CoreField::WIDTH; ++x) {
        for (int y = originalField.height(x); y >= 1; --y) {
            PuyoColor c = originalField.color(x, y);

            if (!isNormalColor(c))
                continue;

            // Drop puyo on
            for (int d = -1; d <= 1; ++d) {

                if (visited[x + d][c])
                    continue;
                visited[x + d][c] = true;

                if (x + d <= 0 || CoreField::WIDTH < x + d)
                    continue;
                if (d == 0) {
                    if (originalField.color(x, y + 1) != PuyoColor::EMPTY)
                        continue;
                } else {
                    if (originalField.color(x + d, y) != PuyoColor::EMPTY)
                        continue;
                }

                CoreField f(originalField);
                int necessaryPuyos = 0;
                while (necessaryPuyos <= 4 && f.countConnectedPuyos(x, y) < 4 && f.height(x + d) <= 13) {
                    f.dropPuyoOn(x + d, c);
                    ++necessaryPuyos;
                }

                if (necessaryPuyos > 4)
                    continue;

                callback(&f, originalField, x + d, c, necessaryPuyos);
            }
        }
    }
}

static inline void tryFloatFire(const CoreField& originalField, SimulationCallback callback)
{
    for (int x = 1; x <= CoreField::WIDTH; ++x) {
        for (int y = originalField.height(x); y >= 1; --y) {
            PuyoColor c = originalField.color(x, y);

            DCHECK(c != PuyoColor::EMPTY);
            if (c == PuyoColor::OJAMA)
                continue;

            int necessaryPuyos = 4 - originalField.countConnectedPuyos(x, y);
            int restPuyos = necessaryPuyos;
            CoreField f(originalField);

            int dx = x - 1;
            // float puyo col dx
            for (; dx <= x + 1 && restPuyos > 0; ++dx) {
                if (dx <= 0 || CoreField::WIDTH < dx) {
                    continue;
                }


                // Check y
                if (dx != x) {
                    if (originalField.color(dx, y) != PuyoColor::EMPTY) {
                        continue;
                    } else { // restPuyos must be more than 0
                        f.unsafeSet(dx, y, c);
                        --restPuyos;
                    }
                }

                int dy_min = y - 1;
                // Check under y
                for (; restPuyos > 0 && dy_min > 0 && originalField.color(dx ,dy_min) == PuyoColor::EMPTY;
                     --dy_min) {
                    f.unsafeSet(dx, dy_min, c);
                    --restPuyos;
                }

                // Check over y
                for (int dy = y + 1;
                     restPuyos > 0 && dy <= 12 && originalField.color(dx ,dy) == PuyoColor::EMPTY; ++dy) {
                    f.unsafeSet(dx, dy, c);
                    --restPuyos;
                }

                // Fill ojama
                for(; dy_min > 0 && originalField.color(dx, dy_min) == PuyoColor::EMPTY; --dy_min) {
                    f.unsafeSet(dx, dy_min, PuyoColor::OJAMA);
                }

                f.recalcHeightOn(dx);
            }

            if (restPuyos <= 0) {
                callback(&f, originalField, dx, c, necessaryPuyos);
            }
        }
    }
}

void findRensas(const CoreField& field, RensaDetector::Mode mode, SimulationCallback callback)
{
    switch (mode) {
    case RensaDetector::Mode::DROP:
        tryDropFire(field, callback);
        break;
    case RensaDetector::Mode::FLOAT:
        tryFloatFire(field, callback);
        break;
    default:
        CHECK(false) << "Unknown mode : " << static_cast<int>(mode);
    }
}


std::vector<FeasibleRensaInfo>
RensaDetector::findFeasibleRensas(const CoreField& field, const KumipuyoSeq& kumipuyoSeq)
{
    std::vector<FeasibleRensaInfo> result;
    Plan::iterateAvailablePlans(field, kumipuyoSeq, kumipuyoSeq.size(), [&result](const RefPlan& plan) {
        if (plan.isRensaPlan())
            result.emplace_back(plan.rensaResult(), plan.initiatingFrames());
    });

    return result;
}

static inline void simulateInternal(CoreField* f, const CoreField& original,
                                    const ColumnPuyoList& keyPuyos, const ColumnPuyoList& firePuyos,
                                    RensaDetector::PossibleRensaCallback callback)
{
    int minHeights[CoreField::MAP_WIDTH] {
        100, original.height(1) + 1, original.height(2) + 1, original.height(3) + 1,
        original.height(4) + 1, original.height(5) + 1, original.height(6) + 1, 100,
    };

    RensaResult rensaResult = f->simulateWithMinHeights(minHeights);
    callback(rensaResult, keyPuyos, firePuyos);
}

static inline void simulateInternal(CoreField* f, const CoreField& original,
                                    const ColumnPuyoList& keyPuyos, const ColumnPuyoList& firePuyos,
                                    RensaDetector::TrackedPossibleRensaCallback callback)
{
    int minHeights[CoreField::MAP_WIDTH] {
        100, original.height(1) + 1, original.height(2) + 1, original.height(3) + 1,
        original.height(4) + 1, original.height(5) + 1, original.height(6) + 1, 100,
    };

    RensaTrackResult rensaTrackResult;
    RensaResult rensaResult = f->simulateAndTrackWithMinHeights(&rensaTrackResult, minHeights);
    callback(rensaResult, keyPuyos, firePuyos, rensaTrackResult);
}

template<typename Callback>
static void findPossibleRensasInternal(const CoreField& field,
                                       const ColumnPuyoList& keyPuyos,
                                       int leftX,
                                       int restAdded,
                                       RensaDetector::Mode mode,
                                       Callback callback)
{
    auto findRensaCallback = [keyPuyos, &callback](CoreField* f, const CoreField& originalField, int x, PuyoColor c, int n) {
        ColumnPuyoList firePuyos;
        for (int i = 0; i < n; ++i)
            firePuyos.addPuyo(x, c);
        simulateInternal(f, originalField, keyPuyos, firePuyos, callback);
    };

    findRensas(field, mode, findRensaCallback);

    if (restAdded <= 0)
        return;

    CoreField f(field);
    ColumnPuyoList puyoList(keyPuyos);

    for (int x = leftX; x <= CoreField::WIDTH; ++x) {
        if (f.height(x) >= 13)
            continue;

        for (int i = 0; i < NUM_NORMAL_PUYO_COLORS; ++i) {
            PuyoColor c = normalPuyoColorOf(i);

            f.dropPuyoOn(x, c);
            puyoList.addPuyo(x, c);

            if (f.countConnectedPuyos(x, f.height(x)) < 4)
                findPossibleRensasInternal(f, puyoList, x, restAdded - 1, mode, callback);

            f.removeTopPuyoFrom(x);
            puyoList.removeLastAddedPuyo();
        }
    }
}

std::vector<PossibleRensaInfo>
RensaDetector::findPossibleRensas(const CoreField& field, int maxKeyPuyos, Mode mode)
{
    std::vector<PossibleRensaInfo> result;
    result.reserve(100000);

    auto callback = [&result](const RensaResult& rensaResult,
                              const ColumnPuyoList& keyPuyos,
                              const ColumnPuyoList& firePuyos) {
        result.emplace_back(rensaResult, keyPuyos, firePuyos);
    };

    ColumnPuyoList puyoList;
    findPossibleRensasInternal(field, puyoList, 1, maxKeyPuyos, mode, callback);
    return result;
}

std::vector<TrackedPossibleRensaInfo>
RensaDetector::findPossibleRensasWithTracking(const CoreField& field, int maxKeyPuyos, Mode mode)
{
    std::vector<TrackedPossibleRensaInfo> result;
    result.reserve(100000);

    auto callback = [&result](const RensaResult& rensaResult, const ColumnPuyoList& keyPuyos,
                              const ColumnPuyoList& firePuyos, const RensaTrackResult& rensaTrackResult) {
        result.emplace_back(rensaResult, keyPuyos, firePuyos, rensaTrackResult);
    };

    ColumnPuyoList puyoList;
    findPossibleRensasInternal(field, puyoList, 1, maxKeyPuyos, mode, callback);
    return result;
}

void RensaDetector::iteratePossibleRensas(const CoreField& field, int maxKeyPuyos,
                                          RensaDetector::PossibleRensaCallback callback, RensaDetector::Mode mode)
{
    ColumnPuyoList puyoList;
    findPossibleRensasInternal(field, puyoList, 1, maxKeyPuyos,mode, callback);
}

void RensaDetector::iteratePossibleRensasWithTracking(const CoreField& field, int maxKeyPuyos,
                                                      RensaDetector::TrackedPossibleRensaCallback callback, RensaDetector::Mode mode)
{
    ColumnPuyoList puyoList;
    findPossibleRensasInternal(field, puyoList, 1, maxKeyPuyos, mode, callback);
}
