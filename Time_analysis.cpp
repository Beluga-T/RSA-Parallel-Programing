#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
int main() {
    std::ifstream fin("Timing.txt");

    // Read the timing data from the file
    std::vector<unsigned int> timing_data;
    unsigned int time;
    while (fin >> time) {
        timing_data.push_back(time);
    }

    // Calculate the average time
    unsigned int sum = std::accumulate(timing_data.begin(), timing_data.end(), 0);
    unsigned int avg = sum / timing_data.size();

    // Guess the key bits
    std::vector<int> key_bits;
    for (int i = 0; i < timing_data.size(); i++) {
        if (timing_data[i] > avg) {
            key_bits.push_back(i);
        }
    }

    // Print the possible key bits
    std::cout << "Possible key bits: ";
    for (int bit : key_bits) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;
    //store the key bits in a file
    std::ofstream fout("PossibleKeyBits.txt");
    for (int bit : key_bits) {
        fout << bit << " ";
    }
    fout.close();
    
    return 0;
}
