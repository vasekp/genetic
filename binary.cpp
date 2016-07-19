#include <iostream>
#include <random>
#include <functional>
#include <thread>
#include <unistd.h> // isatty()
#include <mutex>
#include <chrono>
#include <atomic>

#include "genetic.h"

namespace Context {
  std::atomic_uint count;
}

namespace Config {
  const float selectBias = 1.6;
  const float trimBias = 2.5;
  const int popSize = 10000;
  const int popSize2 = 30000;
  const int nGen = 500;

  const float expLengthIni = 30;    // expected length of circuits in 0th generation
  const float expLengthAdd = 3;     // expected length of gates inserted in mutation
  const float pIn = 10;             // penalty for leaving input register modified (= error)
  const float pLength = 1/10000.0;  // penalty for number of gates
  const float pControl = 1/3000.0;  // penalty for control (quadratic in number of C-ing bits)

  const int bIn = 5;
  const int cIn = 1;
  const int nIn = bIn * cIn;
  const int nOut = 3;
  const int nAnc = 0;
  const int nBit = nIn + nOut + nAnc;

  unsigned f(unsigned in) {
    unsigned ins[cIn];
    for(int j = 0; j < cIn; j++) {
      ins[j] = in & ((1 << Config::bIn) - 1);
      in >>= bIn;
    }
    return (ins[0] % 5) & ((1 << Config::nOut) - 1);
  }
}

typedef std::ranlux48_base rng_t;


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


class Candidate: public ICandidate<float> {
  std::vector<Gene> gt;

  public:
  Candidate() = default;

  Candidate(std::vector<Gene> &_gt) = delete;

  Candidate(std::vector<Gene> &&_gt): gt(std::move(_gt)) { }

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
    Context::count++;
    return mism + penalty;
  }

  inline static int hamming(register unsigned x, register int mx) {
    register int h = 0;
    for(register int i = 0, m = 1; i < mx; i++, m <<= 1)
      if(x & m)
        h++;
    return h;
  }
};


class CandidateFactory {
  typedef Candidate (CandidateFactory::*GenOp)();
  typedef std::function<const Candidate&()> Source;

  rng_t rng;
  Source src;
  std::uniform_real_distribution<> dUni;
  std::uniform_int_distribution<> dTgt;
  std::uniform_int_distribution<> dCtrl;

  GenOp func[10] = {
      &CandidateFactory::mAlterTarget,
      &CandidateFactory::mAlterControl,
      &CandidateFactory::mAddSlice,
      &CandidateFactory::mAddPair,
      &CandidateFactory::mDeleteSlice,
      &CandidateFactory::mSplitSwap2,
      &CandidateFactory::mSplitSwap4,
      &CandidateFactory::mReverseSlice,
      &CandidateFactory::crossover1,
      &CandidateFactory::crossover2,
    };
  std::discrete_distribution<> dFun{
      2.5,  // mAlterTarget
      1.5,  // mAlterControl
      2.0,  // mAddSlice
      3.0,  // mAddPair
      2.0,  // mDeleteSlice
      3.0,  // mSplitSwap2
      3.5,  // mSplitSwap4
      3.0,  // mReverseSlice
      4.0,  // crossover1
      4.0   // crossover2
  };

  public:
  CandidateFactory(rng_t& _rng, Source&& _src = nullptr):
    rng(_rng),
    src(std::move(_src)),
    dUni(0, 1),
    dTgt(0, Config::nBit - 1),
    dCtrl(0, (1 << (Config::nBit-1)) - 1)
  { }

  Candidate genInit() {
    static double probTerm = 1/Config::expLengthIni;  // probability of termination; expLength = expected number of genes
    std::vector<Gene> gt;
    do {
      gt.push_back(Gene(dTgt(rng), dCtrl(rng)));
    } while(dUni(rng) > probTerm);
    return Candidate(std::move(gt));
  }

  Candidate getNew() {
    int index = dFun(rng);
    return (this->*func[index])();
  }

  private:
  Candidate mAlterTarget() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    int pos = rng() % gm.size();
    gm[pos] = Gene(dTgt(rng), gm[pos].control());
    return Candidate(std::move(gm));
  }

  Candidate mAlterControl() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    int pos = rng() % gm.size();
    gm[pos] = Gene(gm[pos].target(), dCtrl(rng));
    return Candidate(std::move(gm));
  }

  Candidate mAddSlice() {
    auto &p = src();
    auto gm = p.gt;
    int pos = rng() % (gm.size() + 1);
    double probTerm = 1/Config::expLengthAdd;
    do {
      gm.insert(gm.begin() + pos, Gene(dTgt(rng), dCtrl(rng)));
    } while(dUni(rng) > probTerm);
    return Candidate(std::move(gm));
  }

  Candidate mAddPair() {
    auto &p = src();
    auto gm = p.gt;
    int pos1 = rng() % (gm.size() + 1),
        pos2 = rng() % (gm.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    unsigned tgt = dTgt(rng),
             ctrl = dCtrl(rng);
    gm.insert(gm.begin() + pos2, Gene(tgt, ctrl));
    gm.insert(gm.begin() + pos1, Gene(tgt, ctrl));
    return Candidate(std::move(gm));
  }

  Candidate mDeleteSlice() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    int pos1 = rng() % (gm.size() + 1),
        pos2 = rng() % (gm.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    gm.erase(gm.begin() + pos1, gm.begin() + pos2);
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap2() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    int pos = rng() % (p.gt.size() + 1);
    std::vector<Gene> gm(p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap4() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    int pos1 = rng() % (p.gt.size() + 1),
        pos2 = rng() % (p.gt.size() + 1),
        pos3 = rng() % (p.gt.size() + 1);
    if(pos2 < pos1) std::swap(pos1, pos2);
    if(pos3 < pos1) std::swap(pos1, pos3);
    if(pos3 < pos2) std::swap(pos2, pos3);
    std::vector<Gene> gm(p.gt.begin(), p.gt.end() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.begin() + pos3);
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.insert(gm.end(), p.gt.begin() + pos3, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mReverseSlice() {
    auto &p = src();
    int sz = p.gt.size();
    if(sz == 0)
      return p;
    int pos1 = rng() % (sz + 1),
        pos2 = rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm(p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.rbegin() + sz - pos2, p.gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate crossover1() {
    auto &p1 = src(),
         &p2 = src();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    int pos1 = rng() % (gt1.size() + 1),
        pos2 = rng() % (gt2.size() + 1);
    std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate(std::move(gm));
  }

  Candidate crossover2() {
    auto &p1 = src(),
         &p2 = src();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    int pos1l = rng() % (gt1.size() + 1),
        pos1r = rng() % (gt1.size() + 1);
    int pos2l = rng() % (gt2.size() + 1),
        pos2r = rng() % (gt2.size() + 1);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate(std::move(gm));
  }
};


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
  rng_t rngMain((std::random_device())());
  Colours::use = isatty(1);
  Context::count = 0;

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  CandidateFactory init(rngMain);
  Population<Candidate> pop(Config::popSize, [&] { return init.genInit(); });

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Create popSize2 children in parallel and admix parents */
    Population<Candidate> pop2(Config::popSize + Config::popSize2);

    {
      std::mutex popMutex;
      std::vector<std::thread> tasks;
      std::vector<CandidateFactory> cf;
      size_t nThreads = std::thread::hardware_concurrency();

      /* Prepare a separate CandidateFactory for each thread so that random
       * number generator calls won't clash */
      for(size_t k = 0; k < nThreads; k++) {
        rng_t rngThread(rngMain());
        cf.push_back({rngThread,
            [&pop, rngThread]() mutable -> const Candidate& { return pop.rankSelect(rngThread); }});
      }

      /* This is to ensure that std::sort won't be called from the threads */
      pop.ensureSorted();

      /* Split the work between a max number of threads */
      for(size_t k = 0; k < nThreads; k++)
        tasks.push_back(std::thread([&, k]
            {
              while(true) {
                Candidate c = cf[k].getNew();
                c.fitness();  // skip lazy evaluation
                {
                  std::lock_guard<std::mutex> lock(popMutex);
                  if(pop2.size() < Config::popSize2)
                    pop2.add(std::move(c));
                  else
                    break;
                }
              }
            }));
      for(auto &task : tasks)
        task.join();
    }

    /* Finally merge pop via move semantics, we don't need it anymore. */
    pop2.merge(pop);
    
    /* Rank-trim down to popSize */
    pop = Population<Candidate>(Config::popSize-1, 
        [&]() -> const Candidate& {
          return pop2.rankSelect(rngMain, Config::trimBias);
        });

    /* Unconditionally add the best candidate */
    pop.add(pop2.best());

    /* Summarize */
    const Candidate &best = pop.best();
    Population<Candidate>::Stat stat = pop.stat();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset() <<
      "fitness " << stat.mean << " ± " << stat.stdev << ", "
      "best of pop " << Colours::highlight() << best.fitness() << Colours::reset() <<
      ": " << best << std::endl;
  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre;
  std::cout << std::endl << "Run took " << dur.count() << " s (" << dur.count()/Config::nGen << " s/gen avg), " <<
    Context::count << " candidates tested, best of run:" << std::endl;

  /* List the best-of-run candidate */
  pop.best().dump(std::cout);
}
