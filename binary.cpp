#include <iostream>
#include <random>

#include "genetic.h"

std::mt19937 Context::rng;

namespace Config {
  const float selectBias = 1.0;
  const float pCrossover = 0.7;
  const int popSize = 50;
  const int popSize2 = 30;
  const int nGen = 50;
}

struct Fitness {
  float x, y;

  operator float() const {
    return std::max(x, y);;
  }

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << (float)f << " (" << f.x << ", " << f.y << ")";
  }
};

class Candidate: public ICandidate<Fitness> {
  float x, y;

  public:
  Candidate(float _x, float _y): x(_x), y(_y) { }

  Candidate() = default;

  static Candidate mutate(const Candidate& p) {
    std::normal_distribution<float> dist(1, 0.1);
    return Candidate{p.x*dist(Context::rng), p.y*dist(Context::rng)};
  }

  static Candidate crossover(const Candidate& p1, const Candidate& p2) {
    return Candidate{p1.x, p2.y};
  }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    return os << "[" << c.x << ", " << c.y << "]";
  }

  private:
  Fitness computeFitness() const {
    return Fitness{x, y};
  }
};

int main() {
  Context::rng = std::mt19937((std::random_device())());

  Population<Candidate> pop;
  {
    std::uniform_real_distribution<float> rDist(0, 1);
    for(int i = 0; i < Config::popSize; i++)
      pop.add(Candidate{rDist(Context::rng), rDist(Context::rng)});
  }

  for(int gen = 0; gen < Config::nGen; gen++) {

    Population<Candidate> pop2 = pop;
    for(int b = 0; b < Config::popSize2; b++)
      if(std::uniform_real_distribution<>(0, 1)(Context::rng) < Config::pCrossover)
        pop2.add(Candidate::crossover(pop.rankSelect(), pop.rankSelect()));
      else
        pop2.add(Candidate::mutate(pop.rankSelect()));

    pop = pop2.trim();

    std::cout << "Gen " << gen << ": best of pop " << pop.best().fitness() << std::endl;
  }
}
