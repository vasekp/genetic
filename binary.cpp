#include <iostream>
#include <random>
#include <functional>
#include <thread>
#include <unistd.h> // isatty()

#include "genetic.h"

std::mt19937 Context::rng;

namespace Config {
  const float selectBias = 2.0;
  const float trimBias = 2.0;
  const int popSize = 10000;
  const int popSize2 = 30000;
  const int nGen = 500;

  const float expLengthIni = 30;    // expected length of circuits in 0th generation
  const float expLengthAdd = 10;    // expected length of gates inserted in mutation
  const float pIn = 10;             // penalty for leaving input register modified (= error)
  const float pLength = 1/10000.0;  // penalty for number of gates
  const float pControl = 1/3000.0;  // penalty for control (quadratic in number of C-ing bits)

  const int bIn = 5;
  const int cIn = 1;
  const int nIn = bIn * cIn;
  const int nOut = 3;
  const int nAnc = 0;
  const int nBit = nIn + nOut + nAnc;

  const int nThreads = 4;

  unsigned f(unsigned in) {
    unsigned ins[cIn];
    for(int j = 0; j < cIn; j++) {
      ins[j] = in & ((1 << Config::bIn) - 1);
      in >>= bIn;
    }
    return (ins[0] % 5) & ((1 << Config::nOut) - 1);
  }
}

namespace Colours {
  bool use;

  const char* bold() { return use ? "\033[1m" : ""; }
  const char* warn() { return use ? "\033[1;33m" : ""; }
  const char* error() { return use ? "\033[1;31m" : ""; }
  const char* highlight() { return use ? "\033[1;32m" : ""; }
  const char* reset() { return use ? "\033[0m" : ""; }
};


class Gene {
  int tgt;
  unsigned ctrlEnc;
  unsigned ctrl;

  public:
  Gene(int _target, unsigned _control): tgt(_target), ctrlEnc(_control) {
    ctrl =
      ((ctrlEnc >> tgt) << (tgt+1)) // shift bits left of tgt to the left
        |
      (ctrlEnc & ((1 << tgt) - 1));    // keep bits right of tgt
  }

  int target() const {
    return tgt;
  }

  int control() const {
    return ctrlEnc;
  }

  unsigned apply(unsigned src) const {
    return ((src & ctrl) == ctrl)
      ? (src ^ (1 << tgt))
      : src;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << g.tgt+1;
    auto c = g.ctrl;
    if(c != 0) {
      os << '[';
      for(int i = 1, m = 1; i <= Config::nBit; i++, m <<= 1)
        if(c & m)
          os << i;
      os << ']';
    }
    return os;
  }
};


inline int hamming(register unsigned x, register int mx) {
  register int h = 0;
  for(register int i = 0, m = 1; i < mx; i++, m <<= 1)
    if(x & m)
      h++;
  return h;
}


class Candidate: public ICandidate<float> {
  std::vector<Gene> gt;

  public:
  Candidate() = default;

  Candidate(std::vector<Gene> _gt): gt(_gt) { }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    for(auto it = c.gt.begin(); it != c.gt.end(); it++)
      os << (it == c.gt.begin() ? "" : " ") << *it;
    return os;
  }

  void dump(std::ostream&) const;

  friend class CandidateFactory;

  private:
  float computeFitness() const {
    register unsigned work;
    unsigned cmp;
    int mism = 0;
    for(int in = 0; in < (1 << Config::nIn); in++) {
      work = in;
      for(const Gene &g : gt)
        work = g.apply(work);
      cmp = in | (Config::f(in) << Config::nIn);
      mism += hamming(work ^ cmp, Config::nBit);
      mism += (Config::pIn - 1) * hamming((work & ((1 << Config::nIn) - 1)) ^ in, Config::nBit);
    }
    float penalty = gt.size()*Config::pLength;
    for(const Gene &g : gt) {
      unsigned h = hamming(g.control(), Config::nBit - 1);
      penalty += h*h*Config::pControl;
    }
    return mism + penalty;
  }
};


class CandidateFactory {
  typedef std::function<Candidate()> GetF;

  struct Op {
    Candidate (*fun)(const GetF&);
    float probRel;
    float probCumm;
  };

  static std::vector<Op> ops;

  public:
  static Candidate getNew(const GetF& get) {
    std::uniform_real_distribution<float> rDist(0, 1);
    float x = rDist(Context::rng);
    int which = 0;
    while(ops[which].probCumm < x)
      which++;
    return ops[which].fun(get);
  }

  static void init() {
    ops = std::vector<Op>{
      {mAlterTarget,  1.5},
      {mAlterControl, 1.5},
      {mAddSlice,     2.0},
      {mAddPair,      3.0},
      {mDeleteSlice,  2.0},
      {mSplitSwap2,   3.0},
      {mSplitSwap4,   3.5},
      {mReverseSlice, 3.0},
      {crossover1,    4.0},
      {crossover2,    4.0}
    };

    float cumm = 0;
    for(auto &op : ops) {
      cumm += op.probRel;
      op.probCumm = cumm;
    }
    for(auto &op : ops)
      op.probCumm /= cumm;
  }

  private:
  static Candidate mAlterTarget(const GetF& get) {
    auto p = get();
    auto gm = p.gt;
    if(gm.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, gm.size() - 1);
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    int pos = dPos(Context::rng);
    gm[pos] = Gene(dTgt(Context::rng), gm[pos].control());
    return Candidate(gm);
  }

  static Candidate mAlterControl(const GetF& get) {
    auto p = get();
    auto gm = p.gt;
    if(gm.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, gm.size() - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    int pos = dPos(Context::rng);
    gm[pos] = Gene(gm[pos].target(), dCtrl(Context::rng));
    return Candidate(gm);
  }

  static Candidate mAddSlice(const GetF& get) {
    auto p = get();
    auto gm = p.gt;
    std::uniform_int_distribution<> dPos(0, gm.size());
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    std::uniform_real_distribution<> dProb(0, 1);
    int pos = dPos(Context::rng);
    float probTerm = 1/Config::expLengthAdd;
    do {
      gm.insert(gm.begin() + pos, Gene(dTgt(Context::rng), dCtrl(Context::rng)));
    } while(dProb(Context::rng) > probTerm);
    return Candidate(gm);
  }

  static Candidate mAddPair(const GetF& get) {
    auto p = get();
    auto gm = p.gt;
    std::uniform_int_distribution<> dPos(0, gm.size());
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    int pos1 = dPos(Context::rng),
        pos2 = dPos(Context::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    unsigned tgt = dTgt(Context::rng),
             ctrl = dCtrl(Context::rng);
    gm.insert(gm.begin() + pos2, Gene(tgt, ctrl));
    gm.insert(gm.begin() + pos1, Gene(tgt, ctrl));
    return Candidate(gm);
  }

  static Candidate mDeleteSlice(const GetF& get) {
    auto p = get();
    auto gm = p.gt;
    if(gm.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, gm.size());
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    int pos1 = dPos(Context::rng),
        pos2 = dPos(Context::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    gm.erase(gm.begin() + pos1, gm.begin() + pos2);
    return Candidate(gm);
  }

  static Candidate mSplitSwap2(const GetF& get) {
    auto p = get();
    if(p.gt.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, p.gt.size());
    int pos = dPos(Context::rng);
    std::vector<Gene> gm(p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate(gm);
  }

  static Candidate mSplitSwap4(const GetF& get) {
    auto p = get();
    if(p.gt.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, p.gt.size());
    int pos1 = dPos(Context::rng),
        pos2 = dPos(Context::rng),
        pos3 = dPos(Context::rng);
    if(pos2 < pos1) std::swap(pos1, pos2);
    if(pos3 < pos1) std::swap(pos1, pos3);
    if(pos3 < pos2) std::swap(pos2, pos3);
    std::vector<Gene> gm(p.gt.begin(), p.gt.end() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.begin() + pos3);
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.insert(gm.end(), p.gt.begin() + pos3, p.gt.end());
    return Candidate(gm);
  }

  static Candidate mReverseSlice(const GetF& get) {
    auto p = get();
    int sz = p.gt.size();
    if(sz == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, sz);
    int pos1 = dPos(Context::rng),
        pos2 = dPos(Context::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm(p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.rbegin() + sz - pos2, p.gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(gm);
  }

  static Candidate crossover1(const GetF& get) {
    auto p1 = get(),
         p2 = get();
    auto gt1 = p1.gt,
         gt2 = p2.gt;
    std::uniform_int_distribution<> dPos1(0, gt1.size()), dPos2(0, gt2.size());
    int pos1 = dPos1(Context::rng),
        pos2 = dPos2(Context::rng);
    std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate(gm);
  }

  static Candidate crossover2(const GetF& get) {
    auto p1 = get(),
         p2 = get();
    auto gt1 = p1.gt,
         gt2 = p2.gt;
    std::uniform_int_distribution<> dPos1(0, gt1.size()), dPos2(0, gt2.size());
    int pos1l = dPos1(Context::rng),
        pos1r = dPos1(Context::rng),
        pos2l = dPos2(Context::rng),
        pos2r = dPos2(Context::rng);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate(gm);
  }
};

std::vector<CandidateFactory::Op> CandidateFactory::ops;


void Candidate::dump(std::ostream& os) const {
  register unsigned work;
  for(int in = 0; in < (1 << Config::nIn); in++) {
    work = in;
    for(const Gene &g : gt)
      work = g.apply(work);
    for(int i = 0; i < Config::cIn; i++) {
      for(int j = 0; j < Config::bIn; j++)
        os << ((in & (1 << ((i+1)*Config::bIn - 1 - j)))?1:0);
      os << ' ';
    }
    os << "→ ";
    for(int i = 0; i < Config::cIn; i++) {
      for(int j = 0; j < Config::bIn; j++) {
        bool bW = (work & (1 << ((i+1)*Config::bIn - 1 - j)))?1:0;
        bool bI = (in & (1 << ((i+1)*Config::bIn - 1 - j)))?1:0;
        if(bW != bI)
          os << Colours::error() << bW << Colours::reset();
        else
          os << bW;
      }
      os << ' ';
    }
    work >>= Config::nIn;
    unsigned cmp = Config::f(in);
    for(int j = Config::nOut - 1; j >= 0; j--) {
      bool bW = (work & (1 << j))?1:0;
      bool bC = (cmp & (1 << j))?1:0;
      if(bW != bC)
        os << Colours::warn() << bW << Colours::reset();
      else
        os << bW;
    }
    os << std::endl;
  }
}


int main() {
  Context::rng = std::mt19937((std::random_device())());
  Colours::use = isatty(1);
  CandidateFactory::init();

  Population<Candidate> pop(Config::popSize);
  {
    float probTerm = 1/Config::expLengthIni; // probability of termination; expLength = expected number of genes
    std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    for(int i = 0; i < Config::popSize; i++) {
      std::vector<Gene> gt;
      do {
        gt.push_back(Gene(dTgt(Context::rng), dCtrl(Context::rng)));
      } while(rDist(Context::rng) > probTerm);
      pop.add(Candidate(gt));
    }
  }

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Create popSize2 children and admix parents */
    Population<Candidate> pop2;
    {
      std::vector<Population<Candidate>> pops(Config::nThreads);
      std::vector<std::thread> tasks;

      /* This is to ensure that std::sort won't be called from the threads */
      pop.ensureSorted();
      for(int k = 0; k < Config::nThreads; k++)
        /* Each thread contributes popSize2/nThreads children. This is subject
         * to rounding errors but we don't care that much about exact numbers. */
        tasks.push_back(std::thread([&, k] {
            pops[k].add(Config::popSize2/Config::nThreads, [&] {
              return CandidateFactory::getNew([&] {
                    return pop.rankSelect();
                  });
              });
            pops[k].precomputeFitnesses();
          }));

      for(auto &task : tasks)
        task.join();

      for(auto pp : pops)
        pop2.merge(pp);
    }
    /* Finally merge pop via move semantics, we don't need it anymore. */
    pop2.merge(pop);

    /* Rank-trim down to popSize */
    pop = Population<Candidate>(Config::popSize-1, [&] {
          return pop2.rankSelect(Config::trimBias);
        });

    /* Unconditionally add the best candidate */
    pop.add(pop2.best());

    /* Summarize */
    Candidate best = pop.best();
    Population<Candidate>::Stat stat = pop.stat();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset() <<
      "fitness " << stat.mean << " ± " << stat.stdev << ", "
      "best of pop " << Colours::highlight() << best.fitness() << Colours::reset() <<
      ": " << best << std::endl;
  }

  /* List the best-of-run candidate */
  pop.best().dump(std::cout);
}
