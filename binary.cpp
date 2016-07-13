#include <iostream>
#include <random>
#include <functional>
#include "unistd.h"

#include "genetic.h"

std::mt19937 Context::rng;

namespace Config {
  const float selectBias = 1;
  const int popSize = 2000;
  const int popSize2 = 5000;
  const int nGen = 500;

  const float expLengthIni = 100;
  const float expLengthAdd = 10;
  const float pIn = 10;
  const float pLength = 1/10000.0;
  const float pControl = 1/3000.0;

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
    //return (ins[0] * ins[1]) & ((1 << Config::nOut) - 1);
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

  unsigned apply(unsigned src) {
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

  Candidate(std::vector<Gene> _gt): gt(_gt) { }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    for(auto it = c.gt.begin(); it != c.gt.end(); it++)
      os << (it == c.gt.begin() ? "" : " ") << *it;
    return os;
  }

  void dump(std::ostream&);

  friend class CandidateFactory;

  private:
  float computeFitness() const {
    register unsigned work;
    unsigned cmp;
    int mism = 0;
    for(int in = 0; in < (1 << Config::nIn); in++) {
      work = in;
      for(Gene g : gt)
        work = g.apply(work);
      cmp = in | (Config::f(in) << Config::nIn);
      mism += hamming(work ^ cmp, Config::nBit);
      mism += (Config::pIn - 1) * hamming((work & ((1 << Config::nIn) - 1)) ^ in, Config::nBit);
    }
    return mism + gt.size()*Config::pLength + Config::pControl*hamming2();
  }

  int hamming2() const {
    int v = 0;
    for(Gene g : gt) {
      unsigned h = hamming(g.control(), Config::nBit - 1);
      v += h*h;
    }
    return v;
  }

  static inline int hamming(register unsigned x, register int mx) {
    register int h = 0;
    for(register int i = 0, m = 1; i < mx; i++, m <<= 1)
      if(x & m)
        h++;
    return h;
  }
};


class CandidateFactory {
  public:
  static Candidate getNew(std::function<Candidate()> get) {
    Candidate (*mutations[])(const Candidate&) = {m1, m2, m3, m4, m5, m6, m7};
    //std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> rDist(0, 7);
    int which = rDist(Context::rng);
    if(which < 7)
      return mutations[which](get());
    else
      return crossover(get(), get());
  }

  private:
  /* Alter target in one gene */
  static Candidate m1(const Candidate& p) {
    auto gm = p.gt;
    if(gm.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, gm.size() - 1);
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    int pos = dPos(Context::rng);
    gm[pos] = Gene(dTgt(Context::rng), gm[pos].control());
    return Candidate(gm);
  }

  /* Alter control in one gene */
  static Candidate m2(const Candidate& p) {
    auto gm = p.gt;
    if(gm.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, gm.size() - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    int pos = dPos(Context::rng);
    gm[pos] = Gene(gm[pos].target(), dCtrl(Context::rng));
    return Candidate(gm);
  }

  /* Add a random operation sequence */
  static Candidate m3(const Candidate& p) {
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

  /* Add an operation pair */
  static Candidate m4(const Candidate& p) {
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

  /* Delete a random slice */
  static Candidate m5(const Candidate& p) {
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

  /* Permute gates (split & swap) */
  static Candidate m6(const Candidate& p) {
    if(p.gt.size() == 0)
      return p;
    std::uniform_int_distribution<> dPos(0, p.gt.size());
    int pos = dPos(Context::rng);
    std::vector<Gene> gm(p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate(gm);
  }

  /* Reverse order of a subsequence of gates */
  static Candidate m7(const Candidate& p) {
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

  /* One-point crossover */
  static Candidate crossover(const Candidate& p1, const Candidate& p2) {
    auto gt1 = p1.gt,
         gt2 = p2.gt;
    std::uniform_int_distribution<> dPos1(0, gt1.size()), dPos2(0, gt2.size());
    int pos1 = dPos1(Context::rng),
        pos2 = dPos2(Context::rng);
    std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate(gm);
  }
};


void Candidate::dump(std::ostream& os) {
  register unsigned work;
  for(int in = 0; in < (1 << Config::nIn); in++) {
    work = in;
    for(Gene g : gt)
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

  Population<Candidate> pop;
  {
    float probTerm = 1/Config::expLengthIni; /* probability of termination; expLength = expected number of genes */
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
    Population<Candidate> pop2 = pop;
    pop2.add(Config::popSize2, [&] {
        return CandidateFactory::getNew([&] {
              return pop.rankSelect();
            });
        });
    pop = Population<Candidate>(Config::popSize, [&] {
          return pop2.rankSelect(2*Config::selectBias);
        });
    Population<Candidate>::Stat stat = pop.stat();

    Candidate best = pop.best();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset() <<
      "fitness " << stat.mean << " ± " << stat.stdev << ", "
      "best of pop " << Colours::highlight() << best.fitness() << Colours::reset() <<
      ": " << best << std::endl;
  }
  pop.best().dump(std::cout);
}
