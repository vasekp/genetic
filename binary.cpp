#include <iostream>
#include <random>

#include "genetic.h"

std::mt19937 Context::rng;

namespace Config {
  const float selectBias = 3;
  const int popSize = 500;
  const int popSize2 = 500;
  const int nGen = 200;

  const float expLengthIni = 100;
  const float expLengthAdd = 5;
  const float pLength = 1/1000.0;
  const float pControl = 1/1000.0;

  const int nIn = 6;
  const int nOut = 3;
  const int nAnc = 3;
  const int nBit = nIn + nOut + nAnc;

  unsigned f(unsigned in) {
    int in1 = in & ((1 << Config::nIn/2) - 1);
    int in2 = (in >> Config::nIn/2) & ((1 << Config::nIn/2) - 1);
    return (in1 + in2) & ((1 << Config::nOut) - 1);
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

  static Candidate m1(const Candidate&);
  static Candidate m2(const Candidate&);
  static Candidate m3(const Candidate&);
  static Candidate m4(const Candidate&);
  static Candidate m5(const Candidate&);
  static Candidate m6(const Candidate&);
  static Candidate crossover(const Candidate&, const Candidate&);

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
      int h = 0;
      for(int i = 0, m = 1; i < Config::nBit - 1; i++, m <<= 1)
        if(ctrl & m)
          h++;
      v += h*h;
    }
    return v;
  }
};



/* Alter target in one gene */
Candidate Candidate::m1(const Candidate& p) {
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
Candidate Candidate::m2(const Candidate& p) {
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
Candidate Candidate::m3(const Candidate& p) {
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
Candidate Candidate::m4(const Candidate& p) {
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
Candidate Candidate::m5(const Candidate& p) {
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
Candidate Candidate::m6(const Candidate& p) {
  if(p.gt.size() == 0)
    return p;
  std::uniform_int_distribution<> dPos(0, p.gt.size());
  int pos = dPos(Context::rng);
  std::vector<Gene> gm(p.gt.begin() + pos, p.gt.end());
  gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
  return Candidate(gm);
}

/* One-point crossover */
Candidate Candidate::crossover(const Candidate& p1, const Candidate& p2) {
  auto gt1 = p1.gt,
       gt2 = p2.gt;
  std::uniform_int_distribution<> dPos1(0, gt1.size()), dPos2(0, gt2.size());
  int pos1 = dPos1(Context::rng),
      pos2 = dPos2(Context::rng);
  std::vector<Gene> gm(gt1.begin(), gt1.begin() + pos1);
  gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
  return Candidate(gm);
}


int main() {
  Context::rng = std::mt19937((std::random_device())());

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
    //std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> rDist(0, 5);
    Candidate (*func[])(const Candidate&) = {Candidate::m1, Candidate::m2, Candidate::m3, Candidate::m4, Candidate::m5, Candidate::m6};
    for(int b = 0; b < Config::popSize2; b++)
      pop2.add(func[rDist(Context::rng)](pop.rankSelect()));

    pop = pop2.trim();

    std::cout << "Gen " << gen << ": best of pop " << pop.best() << std::endl;
  }
}
