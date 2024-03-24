#include <fstream>
#include <gmpxx.h>
#include <string>
#include <iostream>
mpz_class NextPrime(mpz_class current) {
    mpz_nextprime(current.get_mpz_t(), current.get_mpz_t());
    return current;
}

bool IsPrime(mpz_class num) {
    return mpz_probab_prime_p(num.get_mpz_t(), 25) != 0; // 检查是否为质数，进行25轮检测
}

void printProgressBar(int percentage) {
    std::cout << "[";
    int pos = percentage / 5;
    for (int i = 0; i < 20; i++) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << percentage << " %\r";
    std::cout.flush();
}

mpz_class GetLastPrimeFromTXT() {
    std::ifstream file("primes.txt");
    std::string line;
    mpz_class lastPrime = 0;
    bool isEmpty = file.peek() == std::ifstream::traits_type::eof();

    if (isEmpty) {
        return mpz_class("10000000000000000000000000000000000000000000000000"); // 默认值，你可以自己设置
    }

    while (std::getline(file, line)) {
        if (line.empty())
            continue;
        try {
            lastPrime = mpz_class(line.substr(3));
        } catch (const std::invalid_argument& e) {
            // 非法输入，继续读取下一行
            continue;
        }
    }

    return lastPrime;
}


int main() {
    std::ofstream file("primes.txt", std::ios::app); // 以追加模式打开文件
    mpz_class current = GetLastPrimeFromTXT(); // 从文件中获取上次的最后一个质数
    int iterations = 1000000; // 按需调整迭代次数
    mpz_class constantToAdd = mpz_class("1000000000"); // 按需调整常数

    for (int i = 0; i < iterations; i++) {
        mpz_class p = NextPrime(current); // 生成质数p

        // 验证p是否为质数，如果不是，继续下一轮迭代
        if (!IsPrime(p)) {
            current = p;
            continue;
        }

        // 直接在p上添加一个常数，然后生成下一个质数
        mpz_class q = NextPrime(p + constantToAdd);

        // 验证q是否为质数，如果不是，继续下一轮迭代
        if (!IsPrime(q)) {
            current = q;
            continue;
        }

        file << "p: " << p << "\n";
        file << "q: " << q << "\n";

        current = q; // 将q设为下一个迭代的初始值

        printProgressBar(i * 100 / iterations); // 打印进度条
    }

    printProgressBar(100); // 打印完成的进度条
    std::cout << std::endl;

    file.close();
    return 0;
}
