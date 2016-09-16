# Introduction

This is a simple genetic algorithm framework written in C++11 and intending to 
be versatile, parallelized, and simple to use. The encoding, and indeed any 
internal structure of the candidate and fitness, are completely up to the 
particular implementation. Both single- and multi-objective searches are 
supported.


## Prerequisities

The framework makes extensive use of C++11 features and syntax. A fully 
compliant compiler is necessary. GCC 5 or higher is strongly recommended.

OpenMP 2.5 support is required to compile the code without warnings and to use 
the multithreading optimizations. GCC 4 and above are compliant with the 
specification.

Finally, all documentation is provided via [Doxygen][1] markup. Install 
Doxygen to make this accessible via `make doc`.


## Getting started

All the functionality is provided in header files. In order to use the 
framework, start by
```cpp
#include "genetic.hpp"
```
after pointing your compiler to the `include` directory.

Next, define a "candidate base" class providing a **fitness()** function. This 
can return a numeric value or instances of an appropriate structure. In the 
latter case, the return value should provide comparison operators 
corresponding to ordering or dominance, otherwise only very basic 
functionality will be available.

Now you can create a **gen::Population** collecting the candidates:
```cpp
gen::Population<Candidate> pop;
```
This object supports adding, removing candidates, selecting them by various 
strategies, merging populations, some statistics, and more.

Compile with flags `-I` leading to the `include` directory, `-fopenmp` for 
multiprocessor support, and `-O3` for proper performance.

All the functionality is described in the Doxygen documentation embedded in 
the header files. Use `make doc` to extract HTML with hyperlinks into the 
`doc` directory.

Some examples are provided in the `examples` directory. Use `make examples` to 
see them in action.


## Best practices

Move semantics and references are used to speed up the framework internally. 
In order to benefit from the support of the former, use for example
```cpp
pop.add(std::move(pop2));
```
to merge **pop2** into **pop** when the former is no longer needed.

Many functions of **gen::Population** return candidate solutions or their 
groups by reference by default. Use the **auto** keyword to capture these 
easily.

Some variants of **gen::Population** collect metainformation about their 
candidates and need to re-evaluate this if invalidated by changes. It may be 
advantageous to explicitly call for this re-evaluation before using the 
dependent functions. See **gen::NSGAPopulation::precompute** for an example 
scenario of this.

Several function support parallelization via OpenMP tags, and this is on by 
default. Each of such functions accepts a **bool** argument to turn off OpenMP 
if not desirable. This will be the case
- when called from a single worker thread,
- when the sample is too small.
Note that for samples less than about 100 elements, the automatic heuristics 
for splitting the work fairly among worker thread may not work well and the 
overhead imposed by the system can easily outweight any advantage of 
parallelization.

Most of the framework is designed to be thread-safe, using fast user-space 
locking mechanisms to prevent simultaneous writes or reading of an invalid 
memory region. In some cases, however, preventing all possible cases 
automatically would incur significant overhead and has been rejected. The 
following functionality is not guarded by memory locks and must be protected 
by proper program logic:
1. Copying, moving, and assigning of entire **gen::Population** objects
2. Reading from **gen::Population** as a container, namely using 
   iterators, range-based **for**-loops, or
   **gen::BasePopulation::add(const Container&)**.
Destroying a **gen::Population** while read from or written to results in 
a **std::logic_error** being thrown.

[1]: http://www.doxygen.org/index.html
