#include <iostream>
#include <random>
#include <regex>
#include <fstream>

#include "genetic.h"

std::mt19937 Context::rng;

namespace Context {
  std::vector<std::string> dbAccept;  /* strings to be accepted */
  std::vector<std::string> dbReject;  /* strings to be rejected */
  std::vector<std::regex> qchecks;    /* quick validity checks */

  std::string pool =
    "abcdefghijklmnopqrstuvwxyz"
    ".........................."
    ".........................."
    ".........................."  /* . is 3× more likely than a letter */
    "????????||||||||||        "; /* space is 8× more likely than any single letter */
}

namespace Config {
  const float selectBias = 1.0;
  const float pCrossover = 0.7;
  const int popSize = 500;
  const int popSize2 = 300; /* number of children added to the population */
  const int nGen = 300;     /* maximum number of generations in a run */
  const float penalty = 1/2000.; /* for each character of length */
  const float expLength = 50.0;  /* used in generating the initial population */
}


struct Fitness {
  int fAccept;  /* Number of wrongly rejected strings from dbAccept */
  int fReject;  /* Number of wrongly accepted strings from dbReject */
  size_t fLen;  /* Length of ths regex */
  static size_t tAccept;  /* count of dbAccess elements */
  static size_t tReject;  /* count of dbReject elements */

  operator float() const {
    if(fAccept < 0 || fReject < 0)
      return 10;
    else
      return sqrt(fAccept/(float)tAccept) + pow(fReject/(float)tReject, 0.7) + fLen*Config::penalty;
  }

  friend std::ostream& operator<< (std::ostream &os, const Fitness &f) {
    return os << (float)f << " (" << f.fAccept << "/" << f.tAccept << ", " << f.fReject << "/" << f.tReject << ", L" << f.fLen << ")";
  }
};

/* Provide local definitions */
size_t Fitness::tAccept;
size_t Fitness::tReject;


class Candidate: public ICandidate<Fitness> {
  std::string reS;  /* The string representation of this regex */
  std::regex re;    /* The internal representation */
  bool valid;

  public:
  Candidate(std::string _reS): reS(_reS), valid(false) {
    try {
      /* Fails some quick check ⇒ considered not valid */
      for(auto q : Context::qchecks)
        if(std::regex_search(reS, q))
          return;
      re = std::regex(reS, std::regex::extended | std::regex::optimize);
      valid = true;
    }
    catch(std::regex_error e) {
      /* Thrown by std::regex constructor ⇒ not valid */
    }
  }

  /* Needed so that default copy, move, operator= can be used */
  Candidate() = default;

  std::string getRE() const {
    return reS;
  }

  bool isValid() const {
    return valid;
  }

  /* One-point crossover */
  static Candidate crossover(const Candidate& p1, const Candidate& p2) {
    std::string s1 = p1.getRE();
    std::string s2 = p2.getRE();
    int c1 = std::uniform_int_distribution<>(0,s1.length())(Context::rng);
    int c2 = std::uniform_int_distribution<>(0,s2.length())(Context::rng);
    return Candidate(s1.substr(0, c1) + s2.substr(c2));
  }

  /* One-allele mutation */
  static Candidate mutate(const Candidate& p) {
    std::string str = p.getRE();
    int ix = std::uniform_int_distribution<>(0, str.length())(Context::rng);
    str[ix] = Context::pool[
      std::uniform_int_distribution<>(0, Context::pool.length())(Context::rng)];
    return Candidate(str);
  }

  private:
  Fitness computeFitness() const {
    if(!valid)
      return Fitness{-1, -1};
    int cntA = 0, cntR = 0;
    for(auto str : Context::dbAccept)
      if(!std::regex_match(str, re))
        cntA++;
    for(auto str : Context::dbReject)
      if(std::regex_match(str, re))
        cntR++;
    return Fitness{cntA, cntR, reS.length()};
  }
};


int main() {
  Context::rng = std::mt19937((std::random_device())());

  /* Read the database */
  {
    std::string r;
    std::ifstream file("regex.sz");
    while(getline(file, r)) {
      if(r == "---") break;
      Context::dbAccept.push_back(r);
    }
    while(getline(file, r))
      Context::dbReject.push_back(r);
    Fitness::tAccept = Context::dbAccept.size();
    Fitness::tReject = Context::dbReject.size();
  }

  /* Prepare the quick validity checks */
  {
    std::vector<std::string> qchecksS { "[*+?][*+?]", "[^|]*[*+?][^|]*[*+?][^|]*[*+?]" };
    for(auto s : qchecksS)
      Context::qchecks.push_back(std::regex(s, std::regex::extended | std::regex::optimize));
  }

  /* Initialize the G0 population */
  Population<Candidate> pop;
  {
    float probTerm = 1/Config::expLength; /* probability of termination; expLength = expected length of strings */
    std::uniform_real_distribution<float> rDist(0, 1);
    std::uniform_int_distribution<> iDist(0, Context::pool.length());
    for(int i = 0; i < Config::popSize; i++) {
      std::vector<char> vec;
      while(rDist(Context::rng) > probTerm)
        vec.push_back(Context::pool[iDist(Context::rng)]);
      std::string s(vec.begin(), vec.end());
      pop.add(Candidate(s));
    }
  }

  /* The main loop: (l+m) rule */
  for(int gen = 0; gen < Config::nGen; gen++) {
    Population<Candidate> pop2 = pop;
    std::uniform_real_distribution<float> rDist(0, 1);
    for(int j = 0; j < Config::popSize2; j++)
      /* pCrossover: crossover (adding only one child); 1-pCrossover: mutation */
      if(rDist(Context::rng) < Config::pCrossover)
        pop2.add(Candidate::crossover(pop.rankSelect(), pop.rankSelect()));
      else
        pop2.add(Candidate::mutate(pop.rankSelect()));
    pop = pop2.trim();

    int cnt = 0;
    for(auto c : pop)
      if(c.isValid())
        cnt++;
    auto best = pop.best();

    std::cout << "Gen " << gen << ": " << cnt << " valid, "
      << "best of pop: " << best.fitness()
      << "\t" << best.getRE() << std::endl;
  }
}
