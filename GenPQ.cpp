#include <gmpxx.h>
#include <iostream>
#include <cmath>

bool IsPrime(mpz_class num) {
    return mpz_probab_prime_p(num.get_mpz_t(), 25) != 0; // 检查是否为质数，进行25轮检测
}

mpz_class NextPrime(mpz_class current) {
    mpz_nextprime(current.get_mpz_t(), current.get_mpz_t());
    return current;
}

int main() {
    mpz_class iterations;
    std::cout << "Enter the desired growth times for a: ";
    std::cin >> iterations;

    mpz_class base = sqrt(iterations) + iterations;
    mpz_class p = base, q = base;
    
    // 找到p和q，p和q是离base最近的两个质数，且p < base < q
    while(!IsPrime(p)) {
        p--;
    }
    q = NextPrime(base);

    std::cout << "p: " << p << std::endl;
    std::cout << "q: " << q << std::endl;
    
    return 0;
}
