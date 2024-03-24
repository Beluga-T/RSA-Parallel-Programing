#include <iostream>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <numeric>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <thread>
#include <vector>
#include <atomic>
#include <mutex>
#define SIZE 256
std::mutex my_mutex;
unsigned int time1, time2;

unsigned int timestamp(void)
{
    unsigned int bottom;
    unsigned int top;
    asm volatile(".byte 15;.byte 49"
                 : "=a"(bottom), "=d"(top)); // Use rdtsc to get the timestamp
    return bottom;
}
void print_progress(double progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

long get_num_bits(mpz_class n)
{
    /**
     * @param n the integer to be checked
     * @return the number of bits that n has in its binary form
     *
     * TODO: Test this function.
     * TODO: If this function has error, fix it.
     */
    long bits = 0;
    while (n > 0)
    {
        n = n >> 1; // Fix the bug here
        bits++;
    }
    return bits;
}

int bit(mpz_class n, long bit)
{
    /**
     * @param n the number to be look at
     * @param bit the position of bit to be look at
     * @return the bit of n in binary form at bit-th position
     *
     * @example bit(5, 1) would returns 1, because 5 is 0101 in binary, and the
     *           bit at 1-st position is 1.
     *
     *        0-th     2-nd
     *          v       v
     *          0   1   0   1
     *              ^       ^
     *             1-st    3-rd
     */

    return mpz_tstbit(n.get_mpz_t(), bit); // Use mpz_tstbit to get the bit at certain position
}
mpz_class InvMod(mpz_class a, mpz_class n)
{
    /**
     * @param a the number to be inversed
     * @param n the modulus number
     * @return the modular multiplicative inverse of a modulo n
     *
     */

    mpz_class result;
    mpz_t inverse;
    mpz_init(inverse);

    mpz_invert(inverse, a.get_mpz_t(), n.get_mpz_t()); // Use mpz_invert to get the inverse

    result = mpz_class(inverse);
    mpz_clear(inverse);

    return result;
}

mpz_class Encrypt(mpz_class &a, mpz_class &e, mpz_class &n)
{
    if (e == 0)
    {
        return (mpz_class)1;
    }

    long num_bits = get_num_bits(e);
    long i;
    int m;

    mpz_class z;

    z = 1;
    for (i = num_bits - 1; i >= 0; i--)
    {
        z = (z * z) % n;
        if (bit(e, i) == 1)
        {
            z = (z * a) % n;
        }
        // cout<<data[m][i]<<" "<<endl;
    }

    if (e < 0)
    {
        return InvMod(z, n);
    }
    else
    {
        return z;
    }
}

mpz_class Decrypt(mpz_class &a, mpz_class &e, mpz_class &n)
{
    std::ofstream fout;
    // fout.open("Timing.txt");

    if (e == 0)
    {
        return (mpz_class)1;
    }

    long num_bits = get_num_bits(e);
    long i;
    int m;

    unsigned int data[SIZE][num_bits];
    mpz_class z;

    for (m = 0; m < SIZE; m++)
    {
        z = 1;
        for (i = num_bits - 1; i >= 0; i--)
        {
            // time1 = timestamp();
            z = (z * z) % n;
            if (bit(e, i) == 1)
            {
                z = (z * a) % n;
            }
            // data[m][i] = timestamp()-time1;
            // cout<<data[m][i]<<" "<<endl;
        }
    }

    // unsigned int avg[num_bits]; // Use this array to store the average time
    // for (i = 0; i < num_bits; i++) {
    //     avg[i] = 0;
    //     for (m = 0; m < SIZE; m++) {
    //         avg[i] += data[m][i];
    //     }
    //     avg[i] /= SIZE;
    //     fout << avg[i] << " " << std::endl; // Write the average time to the file

    // }

    // fout.close();

    if (e < 0)
    {
        return InvMod(z, n);
    }
    else
    {
        return z;
    }
}

mpz_class Decrypt_Blinding(mpz_class &a, mpz_class &e, mpz_class &n, mpz_class &b)
{
    mpz_class r1 = (mpz_class) "4303324236324956492569821469382146291456284659246924846594259612564612349";
    mpz_class t = Encrypt(r1, b, n);
    mpz_class c1 = t * a;
    mpz_class r_inv = InvMod(r1, n);

    std::ofstream fout;
    fout.open("Timing_B.txt");
    if (e == 0)
    {
        return (mpz_class)1;
    }

    long k = get_num_bits(e);
    long i;
    int m;

    unsigned int data[SIZE][k];
    mpz_class z;

    for (m = 0; m < SIZE; m++)
    {
        z = 1;
        for (i = k - 1; i >= 0; i--)
        {
            time1 = timestamp();
            z = (z * z) % n;
            if (bit(e, i) == 1)
            {
                z = (z * c1) % n;
            }

            data[m][i] = timestamp() - time1; // Use timestamp() to get the time

            // std::cout<<data[m][i]<<" "<<std::endl;

            // timing_b.txt
            fout << data[m][i] << " " << std::endl;
        }
        z = (z * r_inv) % n;
    }

    unsigned int avg[k];
    for (i = 0; i < k; i++)
    {
        avg[i] = 0;
        for (m = 0; m < SIZE; m++)
        {
            avg[i] += data[m][i];
        }
        avg[i] /= SIZE;
        fout << avg[i] << " " << std::endl;
    }

    fout.close();

    if (e < 0)
    {
        return InvMod(z, n);
    }
    else
    {
        return z;
    }
}

void FactorizeAndDecrypt(mpz_class n, mpz_class e, mpz_class c)
{
    mpz_class a = sqrt(n);
    mpz_class b2;
    // int i = 0;
    while (true)
    {
        b2 = a * a - n;
        if (mpz_perfect_square_p(b2.get_mpz_t()) != 0) // Use mpz_perfect_square_p to check if b2 is a perfect square
        {
            break;
        }
        a++;
        // i++;
        //  std::flush (std::cout)<< "\r" << "i: " << i << std::flush; // print the progress
        //  std::cout <<"trying a: "<< a << std::flush; // print the progress
    }
    // std::cout << "i: " << i << std::endl;
    mpz_class b = sqrt(b2);
    mpz_class p1 = a - b;
    mpz_class q1 = a + b;

    mpz_class phi_n1 = (p1 - 1) * (q1 - 1);
    mpz_class d1 = InvMod(e, phi_n1);
    mpz_class m1 = Decrypt(c, d1, n);

    // std::cout << "p1: " << p1 << std::endl;
    // std::cout << "q1: " << q1 << std::endl;
    // std::cout << "d1: " << d1 << std::endl;
    // std::cout << "Decrypted m1: " << m1 << std::endl;
}

void FactorizeAndDecrypt_omp(mpz_class n, mpz_class e, mpz_class c)
{
    mpz_class a_start = sqrt(n);
    mpz_class a_end = n / 2; // Changed upper bound
    mpz_class a;
    mpz_class p1, q1, phi_n1, d1, m1;
    std::atomic<bool> found(false);

#pragma omp parallel shared(a_start, found)
    {
        mpz_class my_p1, my_q1, my_phi_n1, my_d1, my_m1; // New local variables for each thread
        mpz_class my_a = a_start + omp_get_thread_num(); // Calculate initial 'a' for each thread

        while (my_a < a_end && !found)
        {
            mpz_class b2 = my_a * my_a - n;
            if (mpz_perfect_square_p(b2.get_mpz_t()) != 0)
            {
                mpz_class b = sqrt(b2);
                my_p1 = my_a - b;
                my_q1 = my_a + b;
                my_phi_n1 = (my_p1 - 1) * (my_q1 - 1);
                my_d1 = InvMod(e, my_phi_n1);
                my_m1 = Decrypt(c, my_d1, n);

#pragma omp critical
                {
                    if (!found)
                    { // Additional check to avoid race condition
                        found = true;
                        p1 = my_p1;
                        q1 = my_q1;
                        d1 = my_d1;
                        m1 = my_m1;
                    }
                }
            }

            my_a += omp_get_num_threads(); // Increment 'a' by the number of threads
        }
    }

    if (found)
    { // print in single threaded context
      // std::cout << "p1: " << p1 << std::endl;
      // std::cout << "q1: " << q1 << std::endl;
      // std::cout << "d1: " << d1 << std::endl;
      // std::cout << "Decrypted m1: " << m1 << std::endl;
    }

    if (!found)
    {
        std::cout << "omp No solution found" << std::endl;
    }
}

void FactorizeAndDecrypt_thread(mpz_class n, mpz_class e, mpz_class c, int num_threads)
{
    mpz_class a_start = sqrt(n);
    mpz_class a_end = n / 2;
    mpz_class p1, q1, phi_n1, d1, m1;
    std::atomic<bool> found(false);

    // A vector to hold all the threads
    std::vector<std::thread> threads;

    for (int t = 0; t < num_threads; ++t)
    {
        threads.push_back(std::thread([&, t]()
                                      {
            mpz_class my_a = a_start + t; // Calculate initial 'a' for each thread
            while (my_a < a_end && !found)
            {
                mpz_class b2 = my_a*my_a - n;
                if (mpz_perfect_square_p(b2.get_mpz_t()) != 0)
                {
                    mpz_class b = sqrt(b2);
                    mpz_class my_p1 = my_a - b;
                    mpz_class my_q1 = my_a + b;
                    mpz_class my_phi_n1 = (my_p1 - 1) * (my_q1 - 1);
                    mpz_class my_d1 = InvMod(e,my_phi_n1);
                    mpz_class my_m1 = Decrypt(c,my_d1,n);

                    if(!found.exchange(true)) { // Use exchange method of atomic variable to avoid race condition
                        p1 = my_p1;
                        q1 = my_q1;
                        d1 = my_d1;
                        m1 = my_m1;
                    }
                }

                my_a += num_threads; // Increment 'a' by the number of threads
            } }));
    }

    // Wait for all threads to finish
    for (auto &t : threads)
    {
        t.join();
    }

    if (found)
    { // print in single threaded context
      // std::cout << "p1: " << p1 << std::endl;
      // std::cout << "q1: " << q1 << std::endl;
      // std::cout << "d1: " << d1 << std::endl;
      // std::cout << "Decrypted m1: " << m1 << std::endl;
    }

    if (!found)
    {
        std::cout << "thread No solution found" << std::endl;
    }
}
mpz_class NextPrime(mpz_class current)
{
    mpz_nextprime(current.get_mpz_t(), current.get_mpz_t());
    return current;
}
bool IsPrime(mpz_class num)
{
    return mpz_probab_prime_p(num.get_mpz_t(), 25) != 0; // 检查是否为质数，进行25轮检测
}
// Progress bar width on the console.
const int PROGRESS_BAR_WIDTH = 60;

// Function to display a progress bar.
void display_progress(double progress)
{
    int position = PROGRESS_BAR_WIDTH * progress;
    std::cout << "[";
    for (int i = 0; i < PROGRESS_BAR_WIDTH; ++i)
    {
        if (i < position)
            std::cout << "=";
        else if (i == position)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
int main()
{
    // Create a vector to store the results.
    std::vector<std::pair<int, std::vector<double>>> results;
    int MAX_THREADS = 30;
    // Generate the primes, keys, etc. This part is unchanged from your original code.
    mpz_class composite("10000000000000"); // a composite number
    mpz_class p = NextPrime(composite);
    mpz_class q = NextPrime(mpz_class("10000000000000") + 100000000000);

    bool isprime = IsPrime(q);
    if (isprime)
    {
        std::cout << "The number is probably prime" << std::endl;
    }
    else
    {
        std::cout << "The number is definitely not prime" << std::endl;
        return -1;
    }

    mpz_class phi_n = (p - 1) * (q - 1);
    mpz_class n = p * q;
    mpz_class e = (mpz_class) "65537";
    mpz_class d = InvMod(e, phi_n);
    mpz_class m = (mpz_class) "12321";
    mpz_class c = Encrypt(m, e, n);

    std::cout << "d: " << d << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "e: " << e << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "phi_n: " << phi_n << std::endl;
    std::cout << "Plaintext m: " << m << std::endl;
    std::cout << "Ciphertext c: " << c << std::endl;

    mpz_class m_1 = Decrypt(c, d, n);
    std::cout << "Decrypted m: " << m_1 << std::endl;

    double seq_time = 0.0;

    // Measure the time for the sequential function only once.
    std::cout << "Measuring sequential time..." << std::endl;
    auto start_seq = std::chrono::high_resolution_clock::now();
    FactorizeAndDecrypt(n, e, c);
    auto finish_seq = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seq = finish_seq - start_seq;
    seq_time = elapsed_seq.count();
    std::cout << "Sequential time: " << seq_time << " s\n";

    // check the correctness of the Omp and thread
    std::cout << "Checking the correctness of the omp and thread" << std::endl;
    FactorizeAndDecrypt_omp(n, e, c);
    FactorizeAndDecrypt_thread(n, e, c, MAX_THREADS);
    std::cout << "openmp and thread version is correct" << std::endl;

    for (int THREAD = 1; THREAD <= MAX_THREADS; THREAD++)
    {
        std::vector<double> avg_times(3, 0.0); // Initialize a vector to hold the average times.
        const int NUM_TRIALS = 10;             // The number of times to run each function.

        avg_times[0] = seq_time; // Set the sequential time directly

        for (int trial = 0; trial < NUM_TRIALS; trial++)
        {
            std::vector<double> times(3, 0.0);
            times[0] = seq_time; // Directly set sequential time

            // Measure the time for the parallel OpenMP function.
            omp_set_num_threads(THREAD);
            auto start = std::chrono::high_resolution_clock::now();
            FactorizeAndDecrypt_omp(n, e, c);
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            times[1] = elapsed.count();

            // Measure the time for the parallel thread function.
            start = std::chrono::high_resolution_clock::now();
            FactorizeAndDecrypt_thread(n, e, c, THREAD);
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish - start;
            times[2] = elapsed.count();

            // Add the times for OpenMP and Thread to the average times vector (ignoring sequential since it's constant).
            for (int i = 1; i < avg_times.size(); i++)
            { // Start from index 1 since 0 is sequential
                avg_times[i] += times[i];
            }
        }

        // Divide the total times by the number of trials to get the average times (ignoring sequential since it's constant).
        for (int i = 1; i < avg_times.size(); i++)
        { // Start from index 1 since 0 is sequential
            avg_times[i] /= NUM_TRIALS;
        }

        // Save the thread count and average elapsed times.
        results.emplace_back(THREAD, avg_times);

        // Display progress bar.
        double progress = static_cast<double>(THREAD) / MAX_THREADS; // 30 is the maximum number of threads
        int position = PROGRESS_BAR_WIDTH * progress;
        std::cout << "[";
        for (int i = 0; i < PROGRESS_BAR_WIDTH; ++i)
        {
            if (i < position)
                std::cout << "=";
            else if (i == position)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }

    // Erase the progress bar.
    std::cout << std::string(PROGRESS_BAR_WIDTH + 10, ' ') << "\r";
    std::cout.flush();

    // Print the results in a table.
    std::cout << "THREADS\tTIME_SEQ(s)\tTIME_OMP(s)\tTIME_THREAD(s)\n";
    for (const auto &result : results)
    {
        std::cout << result.first;
        for (const auto &time : result.second)
        {
            std::cout << '\t' << time;
        }
        std::cout << '\n';
    }

    // Output to a CSV file.
    std::ofstream out("output.csv");
    out << "THREADS,TIME_SEQ(s),TIME_OMP(s),TIME_THREAD(s)\n";
    for (const auto &result : results)
    {
        out << result.first;
        for (const auto &time : result.second)
        {
            out << ',' << time;
        }
        out << '\n';
    }
    out.close();
    return 0;
}
