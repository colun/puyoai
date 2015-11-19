#include <gflags/gflags.h>
#include <glog/logging.h>

#include <cassert>
#include <random>
#include <tuple>

#include "base/base.h"
#include "core/plan/plan.h"
#include "core/client/ai/ai.h"
#include "core/core_field.h"
#include "core/frame_request.h"
#include "base/time.h"
#include "core/kumipuyo_seq_generator.h"
#include "core/puyo_controller.h"
#include "core/rensa_tracker/rensa_yposition_tracker.h"

unsigned int myRandInt(unsigned int v) {
    static std::mt19937 rnd(1);
    return ((unsigned long long)rnd() * v) >> 32;
}

class MyRensaTracker {
public:
	MyRensaTracker(const CoreField & field) : vanishDropSum_(0) {
		mudaSum_ = 0;
		for(int x=1; x<=6; ++x) {
			int mh = std::max(0, field.height(x) - 3);
			mudaHeight_[x-1] = mh;
			mudaSum_ += mh;
		}
		mudaVanishSum_ = 0;
		mudaVanishLast_ = 0;
		firstVanishHeight_ = -1;
		lowPosBonus_ = 0;
	}

    void trackCoef(int nthChain, int numErasedPuyo, int longBonusCoef, int colorBonusCoef)
    {
        yTracker_.trackCoef(nthChain, numErasedPuyo, longBonusCoef, colorBonusCoef);
    }

    void trackVanish(int nthChain, const FieldBits& vanishedPuyoBits, const FieldBits& vanishedOjamaPuyoBits)
    {
    	int sum = 0;
    	int vanishHeight = 15;
    	int maxVanishHeight = 0;
        vanishedPuyoBits.iterateBitPositions([&](int x, int y) {
        	int oy = yTracker_.originalY(x, y);
        	vanishDropSum_ += oy - y;
        	if(oy<mudaHeight_[x-1]) {
        		--mudaSum_;
        	}
        	++sum;
        	vanishHeight = std::min(vanishHeight, oy);
        	maxVanishHeight = std::max(maxVanishHeight, oy);
        });
        mudaVanishLast_ = std::max(0, (sum-4)*(sum-4));
        mudaVanishSum_ += mudaVanishLast_;
        lowPosBonus_ += (16-maxVanishHeight) * (nthChain+1);
        if(firstVanishHeight_==-1) {
        	firstVanishHeight_ = vanishHeight;
        }

        yTracker_.trackVanish(nthChain, vanishedPuyoBits, vanishedOjamaPuyoBits);
    }

    void trackDrop(FieldBits /*blender*/, FieldBits /*leftOnes*/, FieldBits /*rightOnes*/) {}
#ifdef __BMI2__
    void trackDropBMI2(std::uint64_t /*oldLowBits*/, std::uint64_t /*oldHighBits*/, std::uint64_t /*newLowBits*/, std::uint64_t /*newHighBits*/) {}
#endif

    int getVanishDropSum() const {
    	return vanishDropSum_;
    }
    int getMudaSum() const {
    	return mudaSum_;
    }
    int getMudaVanishSum() const {
    	return mudaVanishSum_;
    }
    int getMudaVanishLast() const {
    	return mudaVanishLast_;
    }
    int getFirstVanishHeight() const {
    	return firstVanishHeight_;
    }
    int getLowPosBonus() const {
    	return lowPosBonus_;
    }
private:
    int vanishDropSum_;
    int mudaSum_;
    int mudaHeight_[6];
    int mudaVanishSum_;
    int mudaVanishLast_;
    int firstVanishHeight_;
    int lowPosBonus_;
    RensaYPositionTracker yTracker_;
};

class ColunAI : public AI {
public:
    ColunAI(int argc, char* argv[]) : AI(argc, argv, "colun") {}
    ~ColunAI() override {}

    DropDecision think(int frameId, const CoreField& f, const KumipuyoSeq& seq,
                       const PlayerState& me, const PlayerState& enemy, bool fast) const override
    {
        long long start_time_ms = currentTimeInMillis();
        static const int maxFrame = 300;
        int frame = maxFrame;
        if(enemy.isRensaOngoing()) {
            frame = enemy.rensaFinishingFrameId() - frameId;
            //fprintf(stderr, "frame: %d\n", frame);
            if(frame<=0) {
                frame = maxFrame;
            }
        }
        UNUSED_VARIABLE(me);
        UNUSED_VARIABLE(fast);

        LOG(INFO) << f.toDebugString() << seq.toString();

        static const Decision DECISIONS[] = {
            Decision(2, 3), Decision(3, 3), Decision(3, 1), Decision(4, 1),
            Decision(5, 1), Decision(1, 2), Decision(2, 2), Decision(3, 2),
            Decision(4, 2), Decision(5, 2), Decision(6, 2), Decision(1, 1),
            Decision(2, 1), Decision(4, 3), Decision(5, 3), Decision(6, 3),
            Decision(1, 0), Decision(2, 0), Decision(3, 0), Decision(4, 0),
            Decision(5, 0), Decision(6, 0),
        };

        int search_turns = 20;
        int time_limit_ms = 600;
        int simCount = 0;
        int counts[22] = {0};
        while(currentTimeInMillis() - start_time_ms < time_limit_ms) {
            ++simCount;
            KumipuyoSeq simSeq = seq;
            if(simSeq.size()<search_turns) {
                simSeq.append(KumipuyoSeqGenerator::generateRandomSequenceWithSeed(search_turns-simSeq.size(), (simCount*1234567891) ^ (frameId*987654321)));
            }
            {
                //simulation
                std::vector<int> bestGenom;
                typedef std::tuple<int, int, int, int, int, int, int, int, int, int, int, int, int, int> TySc;
                TySc bestSc(-1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                for(int tryCount=0; tryCount<200; ++tryCount) {
                    std::vector<int> genom;//(bestGenom.begin(), bestGenom.begin() + myRandInt(bestGenom.size()));
                    CoreField f2 = f;
                    int sc = 0;
                    int maxChain = 0;
                    int mSc = 0;
                    int ff = 0;
                    int sumVanish = 0;
                    int maxChainTurn = 0;
                    UNUSED_VARIABLE(maxChainTurn);
                    int maxChainHeightDiff = 0;
                    int maxChainVanishDropSum = 0;
                    int maxChainMudaSum = 0;
                    int maxChainMudaVanishSum = 0;
                    int maxChainMudaVanishLast = 0;
                    int maxChainBeforeFrame = 0;
                    int maxChainAfterFrame = 0;
                    int maxChainFirstVanishHeight = 0;
                    int maxChainLowPosBonus = 0;
                    int maxChainBalance = 0;
                    bool dead = false;
                    for(int i=0; i<(int)genom.size(); ++i) {
                        auto & de = DECISIONS[genom[i]];
                        int maxHeight = f2.height(1);
                        int minHeight = maxHeight;
                        int goodBalance[] = { 0, 2, 0, -2, -2, 0, 2, 0};
                        int sumHeight = maxHeight - goodBalance[1];
                        for(int j=2; j<=6; ++j) {
                        	int h = f2.height(j);
                        	maxHeight = std::max(maxHeight, h);
                        	minHeight = std::min(minHeight, h);
                        	sumHeight += h - goodBalance[j];
                        }
                        int balance = 0;
                        for(int j=1; j<=6; ++j) {
                        	int h = f2.height(j);
                        	balance += std::abs((h-goodBalance[j])*6-sumHeight);
                        }
                        int dropFrames = f2.framesToDropNext(de);
                        f2.dropKumipuyo(de, simSeq.get(i));
                        MyRensaTracker tracker(f2);
                        const auto & re = f2.simulate(&tracker);
                        if(!f2.isEmpty(3, 12)) {
                            dead = true;
                            break;
                        }
                        if(re.chains) {
                        	++sumVanish;
                        }
                        bool zenkeshi = (f2.isZenkeshi() && i<seq.size());
                        sc += re.score;
                        int beforeFF = ff;
                        ff += dropFrames + re.frames;
                        int newMaxChain = re.chains + (zenkeshi ? 3 : 0) - sumVanish;
                        if(maxChain < newMaxChain) {
                        	maxChain = newMaxChain;
                        	maxChainTurn = i;
                        	maxChainHeightDiff = maxHeight - minHeight;
                        	maxChainVanishDropSum = tracker.getVanishDropSum();
                        	maxChainMudaSum = tracker.getMudaSum();
                        	maxChainMudaVanishSum = tracker.getMudaVanishSum();
                        	maxChainMudaVanishLast = tracker.getMudaVanishLast();
                        	maxChainBeforeFrame = beforeFF;
                        	maxChainAfterFrame = ff;
                        	maxChainFirstVanishHeight = tracker.getFirstVanishHeight();
                        	maxChainLowPosBonus = tracker.getLowPosBonus();
                        	maxChainBalance = balance;
                        }
                        mSc = std::max(mSc, re.score);
                    }
                    if(dead) {
                        continue;
                    }
                    while((int)genom.size()<simSeq.size() && ff<frame) {
                        std::vector<int> candidates;
                        auto & puyo = simSeq.get(genom.size());
                        if(puyo.axis==puyo.child) {
                            for(int v=0; v<11; ++v) {
                                candidates.push_back(v);
                            }
                        }
                        else {
                            for(int v=0; v<22; ++v) {
                                candidates.push_back(v);
                            }
                        }
                        while(true) {
                            if(candidates.empty()) {
                                dead = true;
                                break;
                            }
                            int i = myRandInt(candidates.size());
                            int v = candidates[i];
                            candidates[i] = candidates.back();
                            candidates.pop_back();
                            auto & de = DECISIONS[v];
                            if(!PuyoController::isReachable(f2, de)) {
                                continue;
                            }
                            int maxHeight = f2.height(1);
                            int minHeight = maxHeight;
                            int goodBalance[] = { 0, 2, 0, -2, -2, 0, 2, 0};
                            int sumHeight = maxHeight - goodBalance[1];
                            for(int j=2; j<=6; ++j) {
                            	int h = f2.height(j);
                            	maxHeight = std::max(maxHeight, h);
                            	minHeight = std::min(minHeight, h);
                            	sumHeight += h - goodBalance[j];
                            }
                            int balance = 0;
                            for(int j=1; j<=6; ++j) {
                            	int h = f2.height(j);
                            	balance += std::abs((h-goodBalance[j])*6-sumHeight);
                            }
                            int dropFrames = f2.framesToDropNext(de);
                            if(!f2.dropKumipuyo(de, simSeq.get(genom.size()))) {
                                dead = true;
                                break;
                            }
                            MyRensaTracker tracker(f2);
                            const auto & re = f2.simulate(&tracker);
                            if(!f2.isEmpty(3, 12)) {
                                dead = true;
                                break;
                            }
                            if(re.chains) {
                            	++sumVanish;
                            }
                            sc += re.score;
                            int beforeFF = ff;
                            ff += dropFrames + re.frames;
                            bool zenkeshi = (f2.isZenkeshi() && (int)genom.size()<seq.size());
                            int newMaxChain = re.chains + (zenkeshi ? 3 : 0) - sumVanish;
                            if(maxChain < newMaxChain) {
                            	maxChain = newMaxChain;
                            	maxChainTurn = genom.size();
                            	maxChainHeightDiff = maxHeight - minHeight;
                            	maxChainVanishDropSum = tracker.getVanishDropSum();
                            	maxChainMudaSum = tracker.getMudaSum();
                            	maxChainMudaVanishSum = tracker.getMudaVanishSum();
                            	maxChainMudaVanishLast = tracker.getMudaVanishLast();
                            	maxChainBeforeFrame = beforeFF;
                            	maxChainAfterFrame = ff;
                            	maxChainFirstVanishHeight = tracker.getFirstVanishHeight();
                            	maxChainLowPosBonus = tracker.getLowPosBonus();
                            	maxChainBalance = balance;
                            }
                            mSc = std::max(mSc, re.score);
                            genom.push_back(v);
                            break;
                        }
                        if(dead) {
                            break;
                        }
                    }
                    if(dead) {
                        continue;
                    }
                    TySc sc2(
                    		(maxChainTurn==0 ? maxChain-1 : maxChain) * 1000 - maxChainBalance * 10 + maxChainFirstVanishHeight * 250
                    		, maxChainTurn==0 ? -1000 : 0
							, 0
                    		, 0
							, maxChainLowPosBonus
							, -maxChainMudaVanishSum
							, -maxChainMudaSum
							, -maxChainVanishDropSum
							, -maxChainTurn
							, -maxChainHeightDiff
                    		, -maxChainBeforeFrame
                    		, -maxChainAfterFrame
							, mSc
                    		, 0
							);
                    if(bestSc<sc2) {
                        bestSc = sc2;
                        bestGenom = genom;
                    }
                }
                if(!bestGenom.empty()) {
                    ++counts[bestGenom.front()];
                }
            }
        }
        fprintf(stderr, "simCount: %d\n", simCount);
        int bestAns = 0;
        int bestCnt = 0;
        for(int i=0; i<22; ++i) {
            if(1<=counts[i]) {
                //fprintf(stderr, "%d => %d\n", i, counts[i]);
            }
            if(bestCnt<counts[i]) {
                bestCnt = counts[i];
                bestAns = i;
            }
        }
        return DropDecision(DECISIONS[bestAns]);
    }
};

int main(int argc, char* argv[])
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    google::InstallFailureSignalHandler();

    ColunAI(argc, argv).runLoop();
    return 0;
}
