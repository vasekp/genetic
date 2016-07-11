#include <iostream>
#include <random>

#include "genetic.h"

std::mt19937 Context::rng;

namespace Config {
  const float selectBias = 1.0;
  const int popSize = 500;
  const int popSize2 = 300;
  const int nGen = 100;

  const float expLength = 2;
  const float pLength = 1/20.0;
  const float pControl = 1/100.0;

  const int nIn = 3;
  const int nOut = 1;
  const int nAnc = 0;
  const int nBit = nIn + nOut + nAnc;

  const int rule = 30;  // (a | b) ^ c

  unsigned f(unsigned in) {
    int in1 = in & ((1 << Config::nIn) - 1);  // just for sure
    return (rule & (1 << in1)) ? 1 : 0;
  }
}


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

  Candidate(auto _gt): gt(_gt) { }

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

  /* Add a random operation */
  static Candidate m3(const Candidate& p) {
    auto gm = p.gt;
    std::uniform_int_distribution<> dPos(0, gm.size());
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    int pos = dPos(Context::rng);
    gm.insert(gm.begin() + pos, Gene(dTgt(Context::rng), dCtrl(Context::rng)));
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
    gm.insert(gm.begin() + pos2, Gene(dTgt(Context::rng), dCtrl(Context::rng)));
    gm.insert(gm.begin() + pos1, Gene(dTgt(Context::rng), dCtrl(Context::rng)));
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

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    os << c.fitness() << ':';
    for(Gene g : c.gt)
      os << ' ' << g;
    return os;
  }

  private:
  float computeFitness() const {
    register unsigned work;
    unsigned out;
    int mism = 0;
    for(int in = 0; in < (1 << Config::nIn); in++) {
      work = in;
      for(Gene g : gt)
        work = g.apply(work);
      out = (work >> Config::nIn) & ((1 << Config::nOut) - 1);
      if(out != Config::f(in))
        mism++;
    }
    return mism + gt.size()*Config::pLength + Config::pControl*hamming();
  }

  int hamming() const {
    int v = 0;
    for(Gene g : gt) {
      unsigned ctrl = g.control();
      for(int i = 0, m = 1; i < Config::nBit - 1; i++, m <<= 1)
        if(ctrl & m)
          v++;
    }
    return v;
  }
};


int main() {
  Context::rng = std::mt19937((std::random_device())());

  Population<Candidate> pop;
  {
    float probTerm = 1/Config::expLength; /* probability of termination; expLength = expected number of genes */
    std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> dTgt(0, Config::nBit - 1);
    std::uniform_int_distribution<> dCtrl(0, (1 << (Config::nBit-1)) - 1);
    for(int i = 0; i < Config::popSize; i++) {
      std::vector<Gene> gt;
      while(rDist(Context::rng) > probTerm)
        gt.push_back(Gene(dTgt(Context::rng), dCtrl(Context::rng)));
      pop.add(Candidate(gt));
    }
  }

  for(int gen = 0; gen < Config::nGen; gen++) {
    Population<Candidate> pop2 = pop;
    //std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> rDist(0, 4);
    Candidate (*func[])(const Candidate&) = {Candidate::m1, Candidate::m2, Candidate::m3, Candidate::m4, Candidate::m5};
    for(int b = 0; b < Config::popSize2; b++)
      /*if(rDist(Context::rng) < Config::pMut1)
        pop2.add(Candidate::m1(pop.rankSelect()));
      else
        pop2.add(Candidate::m2(pop.rankSelect()));*/
      pop2.add(func[rDist(Context::rng)](pop.rankSelect()));

    pop = pop2.trim();

    std::cout << "Gen " << gen << ": best of pop " << pop.best() << std::endl;
  }
}
