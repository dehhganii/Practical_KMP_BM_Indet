# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -O3

# Target executables
TARGETS = brute_force kmp_indet bm_indet

all: $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)
