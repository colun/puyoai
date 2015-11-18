#include <gflags/gflags.h>
#include <glog/logging.h>

#include <cassert>
#include <random>
#include <tuple>
#include <queue>

#include "base/base.h"
#include "core/plan/plan.h"
#include "core/client/ai/ai.h"
#include "core/core_field.h"
#include "core/frame_request.h"
#include "base/time.h"
#include "core/kumipuyo_seq_generator.h"
#include "core/puyo_controller.h"

unsigned int myRandInt(unsigned int v) {
    static std::mt19937 rnd(1);
    return ((unsigned long long)rnd() * v) >> 32;
}

typedef std::tuple<int, int, int, int, int> TySc;
const static TySc zeroSc = TySc(-100, 0, 0, 0, 0);

int nodeCreate = 0;
int nodeDelete = 0;

struct MyNode {
	MyNode * parent;
	int v;
	TySc eval;
	int select;
	int refCount;

	MyNode() : parent(NULL), v(-1), eval(zeroSc), select(-1), refCount(1) {
		++nodeCreate;
	}
	MyNode(MyNode * parent_, int v_, TySc eval_, int select_) : parent(parent_), v(v_), eval(eval_), select(select_), refCount(1) {
		++nodeCreate;
		if(parent_) {
			++parent_->refCount;
		}
	}
	void release() {
		if(refCount==1) {
			if(parent) {
				parent->release();
			}
			++nodeDelete;
			delete this;
		}
		else {
			--refCount;
		}
	}
};

struct MyNodeGreater {
	inline bool operator()(MyNode * a, MyNode * b) {
		return a->eval > b->eval;
	}
};

static const Decision DECISIONS[] = {
    Decision(2, 3), Decision(3, 3), Decision(3, 1), Decision(4, 1),
    Decision(5, 1), Decision(1, 2), Decision(2, 2), Decision(3, 2),
    Decision(4, 2), Decision(5, 2), Decision(6, 2), Decision(1, 1),
    Decision(2, 1), Decision(4, 3), Decision(5, 3), Decision(6, 3),
    Decision(1, 0), Decision(2, 0), Decision(3, 0), Decision(4, 0),
    Decision(5, 0), Decision(6, 0),
};

struct MyField {
	CoreField field;
	int frames;
	int sumScore;
	int maxScore;
	int maxChain;
	MyField(const CoreField & f) : field(f), frames(0), sumScore(0), maxScore(0), maxChain(0) {
	}
	bool put(const Kumipuyo & puyo, int v, bool knownFlag) {
		assert(0<=v && v<22);
        auto & de = DECISIONS[v];
        assert(PuyoController::isReachable(field, de));
        int maxHeight = field.height(1);
        int minHeight = maxHeight;
        for(int j=2; j<=6; ++j) {
        	int h = field.height(j);
        	maxHeight = std::max(maxHeight, h);
        	minHeight = std::min(minHeight, h);
        }
        int dropFrames = field.framesToDropNext(de);
        if(!field.dropKumipuyo(de, puyo)) {
        	return false;
        }
        const auto & re = field.simulate();
        if(!field.isEmpty(3, 12)) {
        	return false;
        }
        sumScore += re.score;
        frames += dropFrames + re.frames;
        bool zenkeshi = (field.isZenkeshi() && knownFlag);
        int newMaxChain = re.chains + (zenkeshi ? 3 : 0);
        if(maxChain < newMaxChain) {
        	maxChain = newMaxChain;
        	//maxChainTurn = genom.size();
        	//maxChainHeightDiff = maxHeight - minHeight;
        }
        maxScore = std::max(maxScore, re.score);
        return true;
	}
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

        const static int beamWidth = 5;
        int search_turns = 20;
        int time_limit_ms = 6000;
        int simCount = 0;
        int counts[22] = {0};
        MyField f1(f);
        while(currentTimeInMillis() - start_time_ms < time_limit_ms) {
            ++simCount;
            KumipuyoSeq simSeq = seq;
            if(simSeq.size()<search_turns) {
                simSeq.append(KumipuyoSeqGenerator::generateRandomSequenceWithSeed(search_turns-simSeq.size(), (simCount*1234567891) ^ (frameId*987654321)));
            }
            {
            	//MC Beam
                TySc bestSc(zeroSc);
                std::vector<MyNode*> stack;
                stack.push_back(new MyNode());
                std::priority_queue<MyNode*, std::vector<MyNode*>, MyNodeGreater> que;
                for(int turn=0; turn<10; ++turn) {
                	assert(que.empty());
                	while(!stack.empty()) {
                		MyNode * parent = stack.back();
                		stack.pop_back();
                        std::vector<MyNode*> parents;
                        {
                        	MyNode * node = parent;
                        	while(node) {
                        		parents.push_back(node);
                        		node = node->parent;
                        	}
                        	std::reverse(parents.begin(), parents.end());
                        }
                		MyField f2(f1);
                        for(int i=1; i<(int)parents.size(); ++i) {
                        	bool ret = f2.put(simSeq.get(i), parents[i]->v, i<seq.size());
                        	UNUSED_VARIABLE(ret);
                        	assert(ret);
                        }
                        auto & puyo = simSeq.get(turn);
                        int ej = (puyo.axis==puyo.child ? 11 : 22);
                        for(int j=0; j<ej; ++j) {
                            if(!PuyoController::isReachable(f2.field, DECISIONS[j])) {
                                continue;
                            }
                        	MyField f3(f2);
                        	if(!f3.put(puyo, j, turn<seq.size())) {
                        		continue;
                        	}
                            TySc bestSc(-100, 0, 0, 0, 0);
                            for(int tryCount=0; tryCount<20; ++tryCount) {
                            	MyField f4(f3);
                            	bool dead = false;
                                for(int pos=turn+1; pos<simSeq.size() && f4.frames<frame; ++pos) {
                                    std::vector<int> candidates;
                                    auto & puyo = simSeq.get(pos);
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
                                        if(!PuyoController::isReachable(f4.field, de)) {
                                            continue;
                                        }
                                        f4.put(puyo, v, pos<seq.size());
                                        break;
                                    }
                                    if(dead) {
                                        break;
                                    }
                                }
                                if(dead) {
                                    continue;
                                }
                                TySc sc2(f4.maxChain<3 ? -10+f4.maxChain : f4.maxChain, 0, 0, f4.maxScore, -f4.sumScore);
                                if(bestSc<sc2) {
                                    bestSc = sc2;
                                }
                        	}
                            if(beamWidth<=que.size()) {
                            	if(bestSc<=que.top()->eval) {
                            		continue;
                            	}
                            	que.top()->release();
                            	que.pop();
                            }
                            que.push(new MyNode(parent, j, bestSc, parent->parent==NULL ? j : parent->select));
                        }
                        parent->release();
                	}
                	while(!stack.empty()) {
                		stack.back()->release();
                		stack.pop_back();
                	}
                	if(que.empty()) {
                		break;
                	}
                	int select = que.top()->select;
                	while(!que.empty()) {
                		if(select!=que.top()->select) {
                			select = -1;
                		}
                		stack.push_back(que.top());
                		que.pop();
                	}
                	if(select!=-1) {
                        ++counts[select];
                        break;
                	}
                }
                while(!stack.empty()) {
                	stack.back()->release();
                	stack.pop_back();
                }
            }
        }
        fprintf(stderr, "simCount: %d\n", simCount);
        int bestAns = -1;
        int bestCnt = 0;
        for(int i=0; i<22; ++i) {
            if(1<=counts[i]) {
                fprintf(stderr, "%d => %d\n", i, counts[i]);
            }
            if(bestCnt<counts[i]) {
                bestCnt = counts[i];
                bestAns = i;
            }
        }
        if(bestAns==-1) {
        	bestAns = 0;
        	fprintf(stderr, "NoAnswer.\n");
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
