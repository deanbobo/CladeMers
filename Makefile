CXX ?= c++
# Use the selected compiler's target (AMD64/x86-64 or ARM64/AArch64).
# Keep default flags portable across CPUs; do not add -march=native.
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lz
PREFIX ?= /usr/local
PYTHON ?= python3

.PHONY: all clean test test-blast install
all: clademers

clademers: clademers.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -pthread $< $(LDFLAGS) $(LDLIBS) -o $@

# Remove only named build artifacts, never clademers_* result directories.
clean:
	$(RM) clademers clademers.o
	$(RM) -r clademers.dSYM

# Python's standard library is only needed to run developer tests.
test: clademers
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v

# Example: make test-blast BLAST_BIN=/path/to/ncbi-blast/bin
# Optional: BLAST_TAXDB=/path/to/extracted/taxdb for the real taxonomy test.
test-blast: clademers
	test -n "$(BLAST_BIN)"
	BLAST_BIN="$(BLAST_BIN)" BLAST_TAXDB="$(BLAST_TAXDB)" $(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v

install: clademers
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 clademers "$(DESTDIR)$(PREFIX)/bin/clademers"
