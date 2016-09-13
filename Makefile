export HEADERS = include/genetic.hpp include/genetic_bits/*.hpp

all:	doc examples

doc:	$(HEADERS) Doxyfile README.md
	doxygen Doxyfile

examples:
	make -C examples

clean:
	rm -r doc
	make -C examples clean

.PHONY: all examples clean
