export HEADERS = include/genetic.hpp include/genetic_bits/*.hpp

examples:
	make -C examples

doc:	$(HEADERS) Doxyfile README.md
	doxygen Doxyfile

clean:
	rm -r doc
	make -C examples clean

.PHONY: all examples clean
