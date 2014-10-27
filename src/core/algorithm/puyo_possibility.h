#ifndef CORE_ALGORITHM_PUYO_POSSIBILITY_H_
#define CORE_ALGORITHM_PUYO_POSSIBILITY_H_

#include <glog/logging.h>

#include "core/algorithm/puyo_set.h"

class TsumoPossibility {
public:
    static const int MAX_K = 32;
    static const int MAX_N = 16;

    static double possibility(const PuyoSet& puyoSet, unsigned int k) {
        int a = std::min(MAX_N - 1, puyoSet.red());
        int b = std::min(MAX_N - 1, puyoSet.blue());
        int c = std::min(MAX_N - 1, puyoSet.yellow());
        int d = std::min(MAX_N - 1, puyoSet.green());

        return s_possibility[a][b][c][d][k];
    }

    static int necessaryPuyos(const PuyoSet& puyoSet, double threshold) {
        DCHECK(s_initialized) << "TsumoPossibility is not initialized.";
        DCHECK(0 <= threshold && threshold <= 1.0);

        int a = std::min(MAX_N - 1, puyoSet.red());
        int b = std::min(MAX_N - 1, puyoSet.blue());
        int c = std::min(MAX_N - 1, puyoSet.yellow());
        int d = std::min(MAX_N - 1, puyoSet.green());

        double* p = s_possibility[a][b][c][d];

        for (int k = 0; k < MAX_K; ++k) {
            if (p[k] >= threshold)
                return k;
        }

        return MAX_K;
    }

    static void initialize();

private:
    static bool s_initialized;
    static double s_possibility[MAX_N][MAX_N][MAX_N][MAX_N][MAX_K];
};

#endif
