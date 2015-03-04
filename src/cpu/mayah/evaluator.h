#ifndef EVALUATION_FEATURE_COLLECTOR_H_
#define EVALUATION_FEATURE_COLLECTOR_H_

#include <map>
#include <vector>

#include "complement_book.h"
#include "evaluation_feature.h"
#include "opening_book.h"
#include "pattern_book.h"
#include "score_collector.h"

class ColumnPuyoList;
class CoreField;
class GazeResult;
struct PlayerState;
class RefPlan;
struct RensaResult;
class RensaTrackResult;

class EvaluatorBase {
protected:
    EvaluatorBase(const OpeningBook& openingBook,
                  const ComplementBook& complementBook,
                  const PatternBook& patternBook) :
        openingBook_(openingBook),
        complementBook_(complementBook),
        patternBook_(patternBook) {}

    const OpeningBook& openingBook() const { return openingBook_; }
    const ComplementBook& complementBook() const { return complementBook_; }
    const PatternBook& patternBook() const { return patternBook_; }

private:
    const OpeningBook& openingBook_;
    const ComplementBook& complementBook_;
    const PatternBook& patternBook_;
};

class PreEvalResult {
public:
    PreEvalResult() {}

    const std::vector<int>& matchableOpeningIds() const { return matchableOpeningIds_; }
    const std::vector<int>& matchableComplementIds() const { return matchableComplementIds_; }

    std::vector<int>* mutableMatchableOpeningIds() { return &matchableOpeningIds_; }
    std::vector<int>* mutableMatchableComplementIds() { return &matchableComplementIds_; }

private:
    std::vector<int> matchableOpeningIds_;
    std::vector<int> matchableComplementIds_;
};

class PreEvaluator : public EvaluatorBase {
public:
    PreEvaluator(const OpeningBook& openingBook,
                 const ComplementBook& complementBook,
                 const PatternBook& patternBook) :
        EvaluatorBase(openingBook, complementBook, patternBook) {}

    PreEvalResult preEval(const CoreField& currentField);
};

class MidEvalResult {
public:
    void add(EvaluationFeatureKey key, double value)
    {
        collectedFeatures_[key] = value;
    }

    double feature(EvaluationFeatureKey key) const
    {
        auto it = collectedFeatures_.find(key);
        if (it == collectedFeatures_.end())
            return 0;

        return it->second;
    }

    const std::map<EvaluationFeatureKey, double>& collectedFeatures() const { return collectedFeatures_; }

private:
    std::map<EvaluationFeatureKey, double> collectedFeatures_;
};

class MidEvaluator : public EvaluatorBase {
public:
    MidEvaluator(const OpeningBook& openingBook,
                 const ComplementBook& complementBook,
                 const PatternBook& patternBook) :
        EvaluatorBase(openingBook, complementBook, patternBook) {}

    MidEvalResult eval(const RefPlan&, const CoreField& currentField, double score);
};

class EvalResult {
public:
    constexpr EvalResult(double score, int maxVirtualScore) : score_(score), maxVirtualScore_(maxVirtualScore) {}

    double score() const { return score_; }
    int maxVirtualScore() const { return maxVirtualScore_; }

private:
    double score_;
    int maxVirtualScore_;
};

template<typename ScoreCollector>
class RensaEvaluator : public EvaluatorBase {
public:
    // Don't take ownership of |sc|.
    RensaEvaluator(const OpeningBook& openingBook,
                   const ComplementBook& complementBook,
                   const PatternBook& patternBook,
                   ScoreCollector* sc) :
        EvaluatorBase(openingBook, complementBook, patternBook),
        sc_(sc) {}

    void evalRensaScore(double score, double virtualScore);
    void evalRensaChainFeature(const RensaResult&, const PuyoSet&);
    void collectScoreForRensaGarbage(const CoreField& fieldAfterRensa);
    void evalRensaHandWidthFeature(const CoreField&, const RensaTrackResult&);
    void evalFirePointTabooFeature(const RefPlan&, const RensaTrackResult&);
    void evalRensaIgnitionHeightFeature(const RefPlan&, const RensaTrackResult&, bool enemyHasZenkeshi);
    void evalRensaConnectionFeature(const CoreField& fieldAfterDrop);
    void evalRensaRidgeHeight(const CoreField&);
    void evalRensaValleyDepth(const CoreField&);
    void evalRensaFieldUShape(const CoreField&, bool enemyHasZenkeshi);
    void evalComplementationBias(const ColumnPuyoList& keyPuyos, const ColumnPuyoList& firePuyos);
    void evalRensaStrategy(const RefPlan&, const RensaResult&, const ColumnPuyoList& keyPuyos, const ColumnPuyoList& firePuyos,
                           int currentFrameId, const PlayerState& me, const PlayerState& enemy);
private:
    ScoreCollector* sc_;
};

template<typename ScoreCollector>
class Evaluator : public EvaluatorBase {
public:
    // Don't take ownership of |sc|.
    Evaluator(const OpeningBook& openingBook,
              const ComplementBook& complementBook,
              const PatternBook& patternBook,
              ScoreCollector* sc) :
        EvaluatorBase(openingBook, complementBook, patternBook),
        sc_(sc) {}

    void collectScore(const RefPlan&, const CoreField& currentField, int currentFrameId, int maxIteration,
                      const PlayerState& me, const PlayerState& enemy,
                      const PreEvalResult&, const MidEvalResult&, const GazeResult&);

    // ----------------------------------------------------------------------

    bool evalStrategy(const RefPlan&, const CoreField& currentField, int currentFrameId,
                      const PlayerState& me, const PlayerState& enemy, const GazeResult&,
                      const MidEvalResult&);

    void evalBook(const OpeningBook&, const std::vector<int>& matchableBookIds, const RefPlan&);
    void evalFrameFeature(const RefPlan&);
    void evalRestrictedConnectionHorizontalFeature(const CoreField&);
    void evalThirdColumnHeightFeature(const RefPlan&);
    void evalValleyDepth(const CoreField&);
    void evalRidgeHeight(const CoreField&);
    void evalFieldUShape(const CoreField&, bool enemyHasZenkeshi);
    void evalUnreachableSpace(const CoreField&);

    void evalMidEval(const MidEvalResult&);

    void collectScoreForConnection(const CoreField&);
    void evalCountPuyoFeature(const RefPlan& plan);

private:
    ScoreCollector* sc_;
};

#endif
