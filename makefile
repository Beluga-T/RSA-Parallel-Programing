# Variables
CXX = g++ 
CXXFLAGS = -fopenmp -pthread
LDFLAGS = -lgmp -lgmpxx
SRC = RSA.cpp
OUT = RSA.out

# Default rule
all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT) $(LDFLAGS)

clean:
	rm -f $(OUT)
	