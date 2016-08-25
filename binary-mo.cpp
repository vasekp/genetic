#include <iostream>
#include <random>
#include <functional>
#include <unistd.h> // isatty()
#include <chrono>
#include <atomic>
#include <iomanip>

#ifdef BENCH
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

#include "genetic"


namespace Config {
  const size_t popSize = 100;
  const size_t popSize2 = 2000;
#ifdef BENCH
  const int nGen = 100;
#else
  const int nGen = 500;
#endif

  const float expLengthIni = 30;      // expected length of circuits in 0th generation
  const float expLengthAdd = 1.5;     // expected length of gates inserted in mutation
  const float pDeleteUniform = 0.10;  // probability of single gate deletion 

  const float heurFactor = 0.15;      // how much prior success of genetic ops should influence future choices

  const float pControl = 0.25;        // how much each bit is likely to be a control bit at gate creation

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
}


class Gene {
  unsigned tgt;     // target qubit
  unsigned ctrl;    // bitmask of control bits: 0 through 2^nBit - 1 with zero at tgt bit
  unsigned ctrlEnc; // 0 through UINT_MAX
  uint32_t hw;      // Hamming weight of ctrl

  public:
  NOINLINE Gene(unsigned _target, unsigned _control): tgt(_target), ctrl(0), ctrlEnc(_control), hw(0) {
    double c = (double)_control / std::numeric_limits<unsigned>::max();
    /* Convert an unsigned between 0 and UINT_MAX to a bit string where the
     * probability of 1 in each position is given by Config::pControl. A value
     * less than 0.5 means that plain NOTs and C-NOTs will be generated more
     * often than CC-NOTs and higher. */
    for(int i = 0; i < Config::nBit - 1; i++) {
      ctrl <<= 1;
      if(c < Config::pControl) {
        ctrl |= 1;
        hw++;
        c /= Config::pControl;
      } else {
        c = (c - Config::pControl)/(1 - Config::pControl);
      }
    }
    /* At this point ctrl has nBit-1 bits. We use this to guarantee that
     * 1<<tgt is left unoccupied. */
    ctrl =
      ((ctrl >> tgt) << (tgt+1))  // shift bits left of tgt to the left
        |
      (ctrl & ((1 << tgt) - 1));  // keep bits right of tgt
  }

  unsigned target() const {
    return tgt;
  }

  unsigned control() const {
    return ctrlEnc;
  }

  unsigned weight() const {
    return hw;
  }

  unsigned apply(unsigned src) const {
    return ((src & ctrl) == ctrl)
      ? (src ^ (1 << tgt))
      : src;
  }

  friend std::ostream& NOINLINE operator<< (std::ostream& os, const Gene& g) {
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


struct Fitness {
  unsigned misInput;
  unsigned misOutput;
  size_t length;
  unsigned controls;

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return a.misInput < b.misInput || (a.misInput == b.misInput && a.misOutput <= b.misOutput && a.length <= b.length && a.controls <= b.controls && !(a == b));
  }

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{' << f.misInput << ',' << f.misOutput << ',' << f.length << ',' << f.controls << '}';
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.misInput == b.misInput && a.misOutput == b.misOutput && a.length == b.length && a.controls == b.controls;
  }
};


class Candidate: public gen::Candidate<Fitness> {
  std::vector<Gene> gt{};
  int origin = -1;
  static std::atomic_ulong count;

  public:

  Candidate() = default;

  /* Note to self: copying is suboptimal. This will make the compiler scream at me if I forget. */
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
  NOINLINE Fitness computeFitness() const {
    register unsigned work;
    unsigned cmp, diff;
    unsigned misIn = 0, misOut = 0;
    for(int in = 0; in < (1 << Config::nIn); in++) {
      work = in;
      for(const Gene& g : gt)
        work = g.apply(work);
      cmp = in | (Config::f(in) << Config::nIn);
      diff = work ^ cmp;
      misIn += hamming(diff & ((1 << Config::nIn) - 1));
      misOut += hamming(diff >> Config::nIn);
    }
    unsigned ctrl = 0;
    for(const Gene& g : gt) {
      unsigned h = g.weight();
      ctrl += h*h;
    }
    count++;
    return {misIn, misOut, gt.size(), ctrl};
  }

  inline static uint32_t hamming(register uint32_t x) {
    /* http://stackoverflow.com/a/14555819/1537925 */
    x -= ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
  }
};


typedef gen::Population<Candidate> Population;


class CandidateFactory {
  typedef Candidate (CandidateFactory::*GenOp)();

  std::uniform_real_distribution<> dUni{0, 1};
  std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
  std::uniform_int_distribution<unsigned> dCtrl{};

  static const std::vector<std::pair<GenOp, std::string>> func;
  static std::vector<unsigned> weights;
  std::discrete_distribution<> dFun{};

  Population& pop;

  public:
  CandidateFactory(Population& _pop): pop(_pop) {
    if(weights.size() == 0) {
      weights = std::vector<unsigned>(func.size(), 1);
      normalizeWeights();
    }
    applyWeights();
  }

  static Candidate genInit() {
    const static double probTerm = 1/Config::expLengthIni;  // probability of termination; expLength = expected number of genes
    static thread_local std::uniform_real_distribution<> dUni{0, 1};
    static thread_local std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    static thread_local std::uniform_int_distribution<unsigned> dCtrl{};
    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do {
      gt.push_back(Gene(dTgt(gen::rng), dCtrl(gen::rng)));
    } while(dUni(gen::rng) > probTerm);
    return Candidate(std::move(gt));
  }

  Candidate NOINLINE getNew() {
    int index = dFun(gen::rng);
    return (this->*func[index].first)().setOrigin(index);
  }

  static void hit(int ix) {
    if(ix >= 0)
      weights[ix]++;
  }

  static void normalizeWeights() {
    unsigned total = std::accumulate(weights.begin(), weights.end(), 0);
    float factor = 1/Config::heurFactor * (float)func.size()*Config::popSize2 / total;
    for(auto& w : weights)
      w *= factor;
  }

  void applyWeights() {
    dFun = std::discrete_distribution<>(weights.begin(), weights.end());
  }

  static void dumpWeights(std::ostream& os) {
    float total = std::accumulate(weights.begin(), weights.end(), 0);
    int sz = func.size();
    /* Find the longest GenOp name */
    typedef decltype(func)::value_type cmp;
    auto max = std::max_element(func.begin(), func.end(),
        [](const cmp& v1, const cmp& v2) {
          return v1.second.length() < v2.second.length();
        });
    auto maxw = max->second.length();
    auto _flags = os.flags(std::ios_base::left | std::ios_base::fixed);
    auto _precision = os.precision(4);
    for(int i = 0; i < sz; i++)
      os << std::setw(maxw+3) << func[i].second + ':' << weights[i] / total << std::endl;
    os.flags(_flags);
    os.precision(_precision);
  }

  private:
  const Candidate& get() {
    return pop.randomSelect();
  }

  Candidate mAlterTarget() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(dTgt(gen::rng), gm[pos].control());
    return Candidate(std::move(gm));
  }

  Candidate mAlterControl() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(gm[pos].target(), gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate(std::move(gm));
  }

  Candidate mAlterSingle() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(dTgt(gen::rng), gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate(std::move(gm));
  }

  Candidate mAddSlice() {
    auto &p = get();
    unsigned pos = gen::rng() % (p.gt.size() + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() + ins.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mAddPairs() {
    auto &p = get();
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins;
    ins.reserve(2*Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() + 2*ins.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.insert(gm.end(), std::make_move_iterator(ins.rbegin()), std::make_move_iterator(ins.rend()));
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mDeleteSlice() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() - (pos2 - pos1));
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mDeleteSliceShort() {
    auto &p = get();
    auto sz = p.gt.size();
    if(sz == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1);
    /* Integer with the same distribution in mAddSlice */
    int len = 1 + floor(log(dUni(gen::rng)) / log(1 - 1/Config::expLengthAdd));
    unsigned pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() - (pos2 - pos1));
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mDeleteUniform() {
    auto &p = get();
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    for(auto& g : p.gt)
      if(dUni(gen::rng) >= Config::pDeleteUniform)
        gm.push_back(g);
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap2() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos = gen::rng() % (p.gt.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate(std::move(gm));
  }

  Candidate mSplitSwap4() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1),
             pos3 = gen::rng() % (p.gt.size() + 1);
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

  Candidate mSplitSwap5() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    std::vector<unsigned> pos;
    for(int i = 0; i < 4; i++)
      pos.push_back(gen::rng() % (p.gt.size() + 1));
    std::sort(pos.begin(), pos.end());
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos[0]);
    gm.insert(gm.end(), p.gt.begin() + pos[2], p.gt.begin() + pos[3]);
    gm.insert(gm.end(), p.gt.begin() + pos[1], p.gt.begin() + pos[2]);
    gm.insert(gm.end(), p.gt.begin() + pos[0], p.gt.begin() + pos[1]);
    gm.insert(gm.end(), p.gt.begin() + pos[3], p.gt.end());
    return Candidate(std::move(gm));
  }

  Candidate mReverseSlice() {
    auto &p = get();
    auto sz = p.gt.size();
    if(sz == 0)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
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
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    unsigned pos1 = gen::rng() % (gt1.size() + 1),
             pos2 = gen::rng() % (gt2.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(pos1 + (gt2.size() - pos2));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate(std::move(gm));
  }

  Candidate crossover2() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
    unsigned pos1l = gen::rng() % (gt1.size() + 1),
             pos1r = gen::rng() % (gt1.size() + 1),
             pos2l = gen::rng() % (gt2.size() + 1),
             pos2r = gen::rng() % (gt2.size() + 1);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm;
    gm.reserve(gt1.size() - (pos1r - pos1l) + (pos2r - pos2r));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate(std::move(gm));
  }

  Candidate threeWay() {
    auto &p1 = get(),
         &p2 = get(),
         &p3 = get();
    std::vector<Gene> gm;
    gm.reserve(p1.gt.size() + p2.gt.size() + p3.gt.size());
    gm.insert(gm.end(), p1.gt.begin(), p1.gt.end());
    gm.insert(gm.end(), p2.gt.begin(), p2.gt.end());
    gm.insert(gm.end(), p3.gt.begin(), p3.gt.end());
    return Candidate(std::move(gm));
  }
};


inline void Candidate::dump(std::ostream& os) const {
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

std::vector<unsigned> CandidateFactory::weights;

const std::vector<std::pair<CandidateFactory::GenOp, std::string>> CandidateFactory::func {
    { &CandidateFactory::mAlterTarget,       "MTarget" },
    { &CandidateFactory::mAlterControl,      "MControl" },
    //{ &CandidateFactory::mAlterSingle,       "MSingle" },
    { &CandidateFactory::mAddSlice,          "AddSlice" },
    { &CandidateFactory::mAddPairs,          "AddPairs" },
    { &CandidateFactory::mDeleteSlice,       "DelSlice" },
    { &CandidateFactory::mDeleteSliceShort,  "DelShort" },
    //{ &CandidateFactory::mDeleteUniform,     "DelUnif" },
    //{ &CandidateFactory::mSplitSwap2,        "SpltSwp2"  },
    { &CandidateFactory::mSplitSwap4,        "SpltSwp4"  },
    { &CandidateFactory::mSplitSwap5,        "SpltSwp5"  },
    //{ &CandidateFactory::mReverseSlice,      "InvSlice" },
    { &CandidateFactory::crossover1,         "C/Over1" },
    { &CandidateFactory::crossover2,         "C/Over2" },
    //{ &CandidateFactory::threeWay,           "3some" }
  };


int main() {
#ifdef BENCH
  gen::rng.seed(1);
#endif
  Colours::use = isatty(1);

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  Population pop(Config::popSize, [&] { return CandidateFactory::genInit(); });

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do popSize */
    Population nondom = pop.front();
    //std::cout << nondom.size();
    nondom.prune([](const Candidate& a, const Candidate& b) -> bool {
        return a.fitness() == b.fitness();
      }, Config::popSize);
    nondom.randomTrim(Config::popSize);
    //std::cout << " → " << nondom.size() << std::endl;
    size_t nd = nondom.size();

    for(auto& c : nondom)
      CandidateFactory::hit(c.getOrigin());

    /* Top up to popSize2 candidates */
    Population pop2(Config::popSize2);
    CandidateFactory cf{pop};

#ifndef SINGLE
#pragma omp parallel for schedule(dynamic)
#endif
    for(size_t k = 0; k < Config::popSize2 - nd; k++) {
      Candidate c{cf.getNew()};
      c.fitness();  // skip lazy evaluation
      pop2.add(std::move(c));
    }

    /* Merge the nondominated subset of the previous population */
    pop2.add(std::move(nondom));
    pop = std::move(pop2);

    /* Summarize */
    /*const Candidate &best = pop.best();
    Population::Stat stat = pop.stat();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset() <<
      "fitness " << stat.mean << " ± " << stat.stdev << ", "
      "best of pop " << Colours::highlight() << best.fitness() << Colours::reset() <<
      ": " << best << std::endl;*/
    nondom = pop.front();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset() <<
      nondom.size() << " nondominated";
    if(nondom.size() > 0) {
      const Candidate& e = nondom.randomSelect();
      std::cout << ", e.g. " << e.fitness() << ' ' << e;
    }
    std::cout << std::endl;

    /* Make older generations matter less in the choice of gen. op. */
    CandidateFactory::normalizeWeights();
  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre;
  std::cout << std::endl <<
    "Run took " << dur.count() << " s (" << dur.count()/Config::nGen << " s/gen avg), " <<
    Candidate::totalCount() << " candidates tested" << std::endl;

  /* List the best-of-run candidate */
  /*Population nondom = pop.front();
  std::cout << nondom.size() << " nondominated";
  if(nondom.size() > 0) {
    const Candidate& e = nondom.randomSelect();
    std::cout << ", e.g. " << e.fitness() << ' ' << e << std::endl;;
    e.dump(std::cout);
  } else
    std::cout << std::endl;*/

  /* Dump the heuristic distribution */
  std::cout << std::endl << "Genetic operator distribution:" << std::endl;
  CandidateFactory::dumpWeights(std::cout);

  /* Delete candidates with duplicate fitnesses */
  Population nondom = pop.front();
  std::cout << std::endl << nondom.size() << " nondominated candidates, ";
  nondom.prune([](const Candidate& a, const Candidate& b) -> bool {
      return a.fitness() == b.fitness();
    });
  std::cout << nondom.size() << " unique fitnesses:" << std::endl;
  for(auto& c : nondom)
    std::cout << c.fitness() << ' ' << c << std::endl;
}
