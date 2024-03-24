
# RSA_Timing Project

## Description

The RSA_Timing project is designed to analyze the performance and timing of RSA cryptographic operations. This project includes a variety of tools and scripts, including prime number generators, RSA algorithm implementations, and timing analysis scripts to measure the execution time across different phases of RSA operations.

## Features

- **Prime Generation**: Utilizes efficient algorithms to generate prime numbers required for RSA encryption.
- **RSA Algorithm Implementation**: A complete implementation of RSA encryption and decryption operations.
- **Performance Analysis**: Scripts for analyzing and visualizing the performance of RSA operations under various conditions.
- **Multi-threading Support**: Analyzes the impact of using multiple CPU cores on RSA operation times.

## Prerequisites

Before you can run the RSA_Timing project, you'll need the following installed on your system:

- Python 3.x
- C++ Compiler (GCC or Clang)
- Make
- Git (for cloning the repository)

## Installation

Follow these steps to get your development environment set up:

1. **Clone the repository**

   ```bash
   git clone https://github.com/yourusername/RSA_Timing.git
   cd RSA_Timing
   ```

2. **Compile the C++ components**

   ```bash
   make all
   ```

3. **Install Python dependencies** (Optional)

   If your project relies on Python scripts and has external dependencies, list them in a `requirements.txt` file and install them using pip:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

Describe how to use the project, including any scripts and commands. For example:

To run the prime number generator:

```bash
./PrimeGenerator.out
```

To perform RSA encryption and decryption:

```bash
./RSA.out
```

To analyze RSA operation timings:

```bash
python plot_rsa.py
```


## License

Distributed under the MIT License. See `LICENSE` for more information.

## Contact

Siwei Tan - tswlovexyb@gmail.com


