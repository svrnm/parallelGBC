all: test

fresh:
	$(MAKE) clean;
	$(MAKE) all;

include Makefile.rules

run: test
	$(MAKE) -C test run

.PHONY: test
test: library
	$(MAKE) -C test all

library: src
	ar rcs lib/libf4.a src/*.o

.PHONY: src
src: 
	$(MAKE) -C src all

clean:
	$(MAKE) -C src objclean
	$(MAKE) -C test objclean
	rm -f lib/*.a
