#include <iostream>
#include <random>
#include <functional>
#include <thread>
#include <unistd.h> // isatty()
#include <mutex>
#include <chrono>
#include <atomic>

#include "genetic.h"


namespace ThreadContext {
  typedef std::ranlux48_base rng_t;
  thread_local rng_t rng;
}


namespace Config {
  const float selectBias = 1.0;
  const float trimBias = 2.5;
  const int popSize = 10000;
  const int popSize2 = 30000;
#ifdef BENCH
  const int nGen = 100;
#else
  const int nGen = 500;
#endif

  const float expLengthIni = 30;    // expected length of circuits in 0th generation
  const float expLengthAdd = 3;     // expected length of gates inserted in mutation
  const float pIn = 10;             // penalty for leaving input register modified (= error)
  const float pLength = 1/10000.0;  // penalty for number of gates
  const float pControl = 1/3000.0;  // penalty for control (quadratic in number of C-ing bits)

  const float heurFactor = 0.15;    // how much prior success of genetic ops should influence future choices

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
    //return (ins[0] + ins[1]) & ((1 << Config::nOut) - 1);
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


class Candidate: public ICandidate<float> {
  std::vector<Gene> gt;
  int origin = -1;
  static std::atomic_ulong count;

  public:

  Candidate() = default;

  /* Suboptimal. This will make the compiler scream at me if I forget. */
  Candidate(std::vector<Gene>& _gt) = delete;

  Candidate(std::vector<Gene>&& _gt): gt(std::move(_gt)) { }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    for(auto it = c.gt.begin(); it != c.gt.end(); it++)
      os << (it == c.gt.begin() ? "" : " ") << *it;
    return os;
  }

  Candidate& setOrigin(int _origin) {
    origin = _origin;
    return *this;
  }

  int getOrigin() const {
    return origin;
  }

  void dump(std::ostream&) const;

  static unsigned long totalCount() {
    return count;
  }

  friend class CandidateFactory;

  private:
  float computeFitness() const {
    register unsigned work;
    unsigned cmp;
    int mism = 0;
    for(int in = 0; in < (1 << Config::nIn); in++) {
      work = in;
      for(const Gene& g : gt)
        work = g.apply(work);
      cmp = in | (Config::f(in) << Config::nIn);
      mism += hamming(work ^ cmp);
      mism += (Config::pIn - 1) * hamming((work & ((1 << Config::nIn) - 1)) ^ in);
    }
    float penalty = gt.size()*Config::pLength;
    for(const Gene& g : gt) {
      unsigned h = hamming(g.control());
      penalty += h*h*Config::pControl;
    }
    count++;
    return mism + penalty;
  }

  inline static uint32_t hamming(register uint32_t x) {
    /* http://stackoverflow.com/a/14555819/1537925 */
    x -= ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
  }
};


class CandidateFactory {
  typedef Candidate (CandidateFactory::*GenOp)();
  typedef std::function<const Candidate&()> Source;

  Source src;
  std::uniform_real_distribution<> dUni;
  std::uniform_int_distribution<> dTgt;
  std::uniform_int_distribution<> dCtrl;

  static const std::vector<std::pair<GenOp, std::string>> func;

  static std::vector<int> weights;
  std::discrete_distribution<> dFun;

  public:
  CandidateFactory(Source&& _src = nullptr):
    src(std::move(_src)),
    dUni(0, 1),
    dTgt(0, Config::nBit - 1),
    dCtrl(0, (1 << (Config::nBit-1)) - 1)
  {
    if(weights.size() == 0) {
      weights = std::vector<int>(func.size(), 1);
      normalizeWeights();
    }
    applyWeights();
  }

  Candidate genInit() {
    const static double probTerm = 1/Config::expLengthIni;  // probability of termination; expLength = expected number of genes
    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do {
      gt.push_back(Gene(dTgt(ThreadContext::rng), dCtrl(ThreadContext::rng)));
    } while(dUni(ThreadContext::rng) > probTerm);
    return Candidate(std::move(gt));
  }

  Candidate getNew() {
    int index = dFun(ThreadContext::rng);
    return (this->*func[index].first)().setOrigin(index);
  }

  static void hit(int ix) {
    if(ix >= 0)
      weights[ix]++;
  }

  static void normalizeWeights() {
    int total = std::accumulate(weights.begin(), weights.end(), 0);
    float factor = 1/Config::heurFactor * (float)func.size()*Config::popSize / total;
    for(int& w : weights)
      w *= factor;
  }

  void applyWeights() {
    dFun = std::discrete_distribution<>(weights.begin(), weights.end());
  }

  static void dumpWeights(std::ostream& os) {
    float total = std::accumulate(weights.begin(), weights.end(), 0);
    int sz = func.size();
    for(int i = 0; i < sz; i++)
      os << func[i].second << ":\t" << weights[i] / total << std::endl;
  }

  private:
  Candidate mAlterTarget() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    int pos = ThreadContext::rng() % p.gt.size();
    gm[pos] = Gene(dTgt(ThreadContext::rng), gm[pos].control());
    return Candidate(std::move(gm));
  }

  Candidate mAlterControl() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    int pos = ThreadContext::rng() % p.gt.size();
    gm[pos] = Gene(gm[pos].target(), dCtrl(ThreadContext::rng));
    return Candidate(std::move(gm));
  }

  Candidate mAddSlice() {
    auto &p = src();
    int pos = ThreadContext::rng() % (p.gt.size() + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dTgt(ThreadContext::rng), dCtrl(ThreadContext::rng));
    } while(dUni(ThreadContext::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() + ins.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mAddPair() {
    auto &p = src();
    int pos1 = ThreadContext::rng() % (p.gt.size() + 1),
        pos2 = ThreadContext::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    unsigned tgt = dTgt(ThreadContext::rng),
             ctrl = dCtrl(ThreadContext::rng);
    gm.reserve(p.gt.size() + 2);
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.emplace_back(tgt, ctrl);
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.emplace_back(tgt, ctrl);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mDeleteSlice() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    int pos1 = ThreadContext::rng() % (p.gt.size() + 1),
        pos2 = ThreadContext::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() - (pos2 - pos1));
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap2() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    int pos = ThreadContext::rng() % (p.gt.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap4() {
    auto &p = src();
    if(p.gt.size() == 0)
      return p;
    int pos1 = ThreadContext::rng() % (p.gt.size() + 1),
        pos2 = ThreadContext::rng() % (p.gt.size() + 1),
        pos3 = ThreadContext::rng() % (p.gt.size() + 1);
    if(pos2 < pos1) std::swap(pos1, pos2);
    if(pos3 < pos1) std::swap(pos1, pos3);
    if(pos3 < pos2) std::swap(pos2, pos3);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
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
    int pos1 = ThreadContext::rng() % (sz + 1),
        pos2 = ThreadContext::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.rbegin() + sz - pos2, p.gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate crossover1() {
    auto &p1 = src(),
         &p2 = src();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    int pos1 = ThreadContext::rng() % (gt1.size() + 1),
        pos2 = ThreadContext::rng() % (gt2.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(pos1 + (gt2.size() - pos2));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate(std::move(gm));
  }

  Candidate crossover2() {
    auto &p1 = src(),
         &p2 = src();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    int pos1l = ThreadContext::rng() % (gt1.size() + 1),
        pos1r = ThreadContext::rng() % (gt1.size() + 1);
    int pos2l = ThreadContext::rng() % (gt2.size() + 1),
        pos2r = ThreadContext::rng() % (gt2.size() + 1);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm;
    gm.reserve(gt1.size() - (pos1r - pos1l) + (pos2r - pos2r));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate(std::move(gm));
  }
};


void Candidate::dump(std::ostream& os) const {
  register unsigned work;
  for(int in = 0; in < (1 << Config::nIn); in++) {
    work = in;
    for(const Gene& g : gt)
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


std::atomic_ulong Candidate::count { 0 };

std::vector<int> CandidateFactory::weights;

const std::vector<std::pair<CandidateFactory::GenOp, std::string>> CandidateFactory::func {
    { &CandidateFactory::mAlterTarget,  "MTarget" },
    { &CandidateFactory::mAlterControl, "MControl" },
    { &CandidateFactory::mAddSlice,     "AddSlice" },
    { &CandidateFactory::mAddPair,      "AddPair" },
    { &CandidateFactory::mDeleteSlice,  "DelSlice" },
    { &CandidateFactory::mSplitSwap2,   "SpltSwp2"  },
    { &CandidateFactory::mSplitSwap4,   "SpltSwp4"  },
    { &CandidateFactory::mReverseSlice, "InvSlice" },
    { &CandidateFactory::crossover1,    "C/Over1" },
    { &CandidateFactory::crossover2,    "C/Over2" }
  };


int main() {
#ifdef BENCH
  ThreadContext::rng = ThreadContext::rng_t(1);
#else
  ThreadContext::rng = ThreadContext::rng_t((std::random_device())());
#endif
  Colours::use = isatty(1);
  CandidateFactory init;

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  Population<Candidate> pop(Config::popSize, [&] { return init.genInit(); });

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Create popSize2 children in parallel and admix parents */
    Population<Candidate> pop2(Config::popSize + Config::popSize2);

    {
      std::mutex popMutex;
#ifndef BENCH
      std::vector<std::thread> tasks;
      /* Split the work between a max number of threads */
      for(size_t k = 0; k < std::thread::hardware_concurrency(); k++) {

        auto seed = ThreadContext::rng();
        /* Let each thread keep adding candidates until the goal is met */
        tasks.push_back(std::thread([&, seed]
            {
              /* Prepare a separate CandidateFactory for each thread so that random
               * number generator calls won't clash */
              ThreadContext::rng = ThreadContext::rng_t(seed);
#endif
              CandidateFactory cf([&]() -> const Candidate& { return pop.rankSelect(ThreadContext::rng, Config::selectBias); });
              while(true) {
                Candidate c = cf.getNew();
                c.fitness();  // skip lazy evaluation
                {
                  std::lock_guard<std::mutex> lock(popMutex);
                  if(pop2.size() < Config::popSize2)
                    pop2.add(std::move(c));
                  else
                    break;
                }
              }
#ifndef BENCH
            }));
      }
      for(auto& task : tasks)
        task.join();
#endif
    }

    /* Finally merge pop via move semantics, we don't need it anymore. */
    pop2.merge(pop);
    
    /* Rank-trim down to popSize */
    pop = Population<Candidate>(Config::popSize-1, 
        [&]() -> const Candidate& {
          const Candidate &c = pop2.rankSelect(ThreadContext::rng, Config::trimBias);
          CandidateFactory::hit(c.getOrigin());
          return c;
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

    /* Make older generations matter less in the choice of gen. op. */
    CandidateFactory::normalizeWeights();
  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre;
  std::cout << std::endl << "Run took " << dur.count() << " s (" << dur.count()/Config::nGen << " s/gen avg), " <<
    Candidate::totalCount() << " candidates tested, best of run:" << std::endl;

  /* List the best-of-run candidate */
  pop.best().dump(std::cout);

  /* Dump the heuristic distribution */
  std::cout << std::endl << "Genetic operator distribution:" << std::endl;
  CandidateFactory::dumpWeights(std::cout);
}
