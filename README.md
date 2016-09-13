# Introduction

This is a simple genetic algorithm framework written in C++11 and intending to 
be generic, parallelized, and simple to use. The encoding, and indeed any 
internal structure of the candidate and fitness, are completely up to the 
particular implementation. Both single- and multi-objective searches are 
supported.

The header currently declares a virtual `Candidate` class template, meant to be 
specialized and derived from, and a `Population` class template supporting 
adding candidates one by one or using a generating function, merging 
populations, several methods of selection, some statistics, and more.

Some examples are provided in the `examples` directory.
