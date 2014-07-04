#include "ai.h"

#include <algorithm>
#include <tuple>
#include <vector>

#include "core/algorithm/plan.h"
#include "core/algorithm/rensa_info.h"
#include "core/algorithm/rensa_detector.h"
#include "core/field/core_field.h"
#include "core/field/rensa_result.h"

munetoshi::AI::AI() : ::AI("munetoshi") { strategy = GROW; }

void munetoshi::AI::gameWillBegin(const FrameData& frame) {
  UNUSED_VARIABLE(frame);
  strategy = GROW;
}

DropDecision munetoshi::AI::think(int frame_id, const PlainField& field,
                                  const KumipuyoSeq& seq) {
  return think_internal(frame_id, CoreField(field), seq);
}

DropDecision munetoshi::AI::think_internal(int frame_id,
                                           const CoreField& field,
                                           const KumipuyoSeq& seq) {
  UNUSED_VARIABLE(frame_id);

  Decision best_chain_decision;
  Decision best_fire_decision;
  int best_chain_grade = -1;
  int best_fire_grade = -1;
  int previous_chain_grade = evaluate(field);
  auto dicisionMaker = [&](const RefPlan& plan) {
    int fire_grade = plan.rensaResult().score;
    int chain_grade = evaluate(plan.field());

    if (best_fire_grade < fire_grade) {
      best_fire_grade = fire_grade;
      best_fire_decision = plan.decisions().front();
    }

    if (best_chain_grade < chain_grade) {
      best_chain_grade = chain_grade;
      best_chain_decision = plan.decisions().front();
    }
  };

  Plan::iterateAvailablePlans(field, seq, 2, dicisionMaker);
  return best_chain_grade < previous_chain_grade
    || (strategy == FIRE && best_fire_grade > 1000)
             ? DropDecision(best_fire_decision)
             : DropDecision(best_chain_decision);
}

void munetoshi::AI::enemyGrounded(const FrameData& frame) {
  CoreField field(frame.enemyPlayerFrameData().field);
  field.forceDrop();

  if (field.simulate().chains > 1) {
    strategy = FIRE;
  }
}

int munetoshi::AI::evaluate(const CoreField& field) {
  int grade = -1;
  int sum;
  auto adder = [&](std::tuple<int, PuyoColor> t) { sum += std::get<0>(t); };
  std::vector<TrackedPossibleRensaInfo> rensa_info_vect =
      RensaDetector::findPossibleRensasWithTracking(field, 1, RensaDetector::Mode::FLOAT);
  for (auto i = rensa_info_vect.begin(); i != rensa_info_vect.end(); ++i) {
    sum = 0;
    auto key_puyos = i->keyPuyos.list();
    auto fire_puyos = i->firePuyos.list();
    std::for_each(key_puyos.rbegin(), key_puyos.rend(), adder);
    std::for_each(fire_puyos.rbegin(), fire_puyos.rend(), adder);
    grade = std::max(grade, i->rensaResult.chains * 10 - sum * 3
        - std::max(field.height(2) - field.height(1) - 2, 0) * 3
        - std::max(field.height(3) - field.height(2) - 2, 0) * 3
        - std::max(field.height(4) - field.height(3) - 2, 0) * 3
        - std::max(field.height(4) - field.height(5) - 2, 0) * 3
        - std::max(field.height(5) - field.height(6) - 2, 0) * 3);

  }
  return std::max(grade, 0);
}
