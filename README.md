# Introduction

This is a simple genetic algorithm framework written in C++11 and intending to 
be generic, thread-safe and simple to use. The encoding, and indeed any 
internal structure of the candidate and fitness, are completely up to the 
particular implementation. At the moment a total ordering of candidates must 
be available, however, limiting the use to single-objective optimization 
schemes.

The header currently declares a virtual ICandidate class template, meant to be 
specialized and derived from, and a Population class template supporting 
adding candidates one by one or using a generating function, merging 
populations, some statistics and more. An algorithm is provided for rank-based 
selection with tunable exponential decay of selection probability in rank but 
other selection methods can be provided.

Two examples are provided: a regular expression-based string classifier and 
an approximate reversible binary quantum circuit search.

# RegEx

This demo searches for a regular expression (regex) separating two sets of 
strings.  For illustratory purposes the [list of cites and towns in Czech 
Republic](https://cs.wikipedia.org/wiki/Seznam_m%C4%9Bst_v_%C4%8Cesku_podle_po%C4%8Dtu_obyvatel) 
sorted by the population was taken, split at an arbitrary point of 50000 
people and normalized (put into lowercase with diacritics removed).

The task is solved using a generational (l+m) scheme with rank-based 
selection, single-point crossover and single-gene random-allele mutation. The 
encoding is direct (the string describing the regular expression). Invalid 
regexes and those considered too complex (by a short series of quick checks) 
are disregarded by being assigned a prohibitive fitness. The algorithm soon 
naturally moderates the population so that little to none invalid candidates 
are generated.

The fitness is based on the number of false negative and false positive 
matches and the length of the regex, but ultimately reduced to one floating 
point number.  Anyway it is internally maintained as a data structure which 
converts to `float` on demand.

## Example solutions

After 300 generations in a (500+300) scheme with a crossover probability 0.7 
(mutation probability 0.3) and penalty 1/2000 per character, the following 
candidate solutions have been obtained:

`....r..|||..b..o|...i.......?...|.....o|...t|b...|.|||j...?...|..l.y.o..|.r.........?...||.|.|o....|t.k||o...?....|..l.y.o|...n||....n||pa.......|..p..?..|.?|..|..g.||.r...|..r.g?.....?...|........?...j...?...|..l.y..|.ar.g?...`

This accepts all 21 towns over 50000 and rejects correctly 538 out of 582 
smaller towns, at a length of 227 characters. Clearly some parts of the regex 
are repetitive or useless. Longer run of the algorithm could eliminate this.

`......?v.|......  .....?v.||.a....p?i..||.|.i...l|....p?i..|.....a?u.|.r..........|..r........|..z?..|.|.l||.|p?i..|..a...|d|d....|..af.|...?.?......a...|..a?..|......b........|....r..`

This also correctly accepts all larger towns but is less efficient in the rest 
(accepting 68 out of 581 as false positives) at the advantage of a shorter 
length (186 characters).


# Binary

This demo searches for a quantum circuit composed entirely of `C^k-NOT` 
operations aiming to compute a given binary function on `(cIn * bIn) + nOut`
qubits, assuming that the output register was initialized to zero.

This is solved using a generational (l+m) scheme with rank-based selection and 
rank-based elitist trimming. The encoding is indirect, preventing invalid 
quantum gates. The genetic operations are plentiful, initially chosen from 
a uniform distribution but gradually biased more towards those having produced 
higher quality solutions:

* single-gene random-allele mutation (target qubit or control set)
* slice addition
* mirror pair slice addition
* slice deletion
* AB → BA permutation
* ABCD → ACBD permutation
* ABCDE → ADCBE permutation
* variable-k-inversion
* single-point crossover
* two-point crossover
* three-way full concatenation

The metaheuristics showed some how each of these operations is influential 
and, in some cases, how it can be enhanced (e.g., single pair addition was 
amended to mirror slice pair addition). It also revealed how 
variable-k-inversion was actually detrimental to performance, polluting the 
population by equivalent candidates of unaffected fitness, so it was removed 
from the selection in recent runs. Most importantly, the analysis showed that 
both crossover operations had little effect on the quality of the sample in 
this case, so the evolution is almost purely mutation-driven. This might be 
understood as that crossovers are not very beneficial in this scenario in 
general, or that good crossover operations are hard to design. For the time 
being, both the operations have been kept anyway.

The fitness function is based on the number of affected input bits (highly 
penalized: the input register should stay unchanged), number of incorrect 
output bits across all possible inputs, length of the circuit, and complexity 
of the control structure (number of control bits per gate), in decreasing 
order of importance. Internally only the float value is maintained.

The function to be compared with can be entered by hardcoding, along with the 
number of input register sections, bit width of one section, and number of 
output qubits. Historically also number of ancilla bits can be specified but 
experience has shown they are not necessary.

This program is more involved than the previous example, featuring parallel 
candidate generation, code optimization where appropriate, fast random number 
generation, and nicer output. However the underlying library is shared.

## Example solutions

Some functions are inherently easier to find. The algorithm converged to 
perfect solutions for addition and multiplication within a few tens of 
generations. For example, for the multiplication problem with two 3-bit 
inputs with a low-3-bit output, this flawless solution was found in 22-th
generation (see below for other parameters):

`8[15] 9[16] 8[24] 9[34] 9[25] 7[14] 9[257]`

Here the number before the square brackets denotes the target qubit of the 
controlled `NOT` operation and the digits enclosed in the brackets are the 
control qubits. Numbers 1–3 and 4–6 are the two sections of the input register 
and 7–9 denotes the output register. We can see that the circuit never acts 
on the input register, which was not enforced in any way in the initialization 
or by the genetic operators used.

Another solution without the need for `C^3-NOT`, found after 25 generations:

`7[14] 8[24] 2[18] 9[25] 2[18] 9[16] 8[15] 9[34]`

Interestingly, this circuit “stores” some information in the second input bit 
and uncomputes this step later on. This, again, was a result purely of the 
evolutionary algorithm.

The evolution works similarly well with an analogous problem of 3-bit 
incomplete addition, where the following zero-error example solutions have 
been found within 100 generations:

`9[3] 8[14] 7[4] 8[2] 9[58] 9[6] 9[127] 7[1] 8[5]`

`8[14] 9[28] 7[4] 8[2] 9[58] 9[3] 8[5] 9[6] 7[1]`

A deliberately harder test function was chosen to be `x mod 5` with one 5-bit 
input register and a 3-bit output. After 500 generations of a (10000+30000) 
scheme candidate solutions with as low as 2 to 5 bit errors out of 96 (`2^5 *
3`) have frequently been found (along with several zero-error solutions).
These are too long to list, except one surprisingly elegant zero-error
exception:

`9[25] 5[12] 9[145] 5[12] 7[1] 8 9[38] 9[6] 8 7[4] 8[2] 8[5] 8[14]`
