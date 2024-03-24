import matplotlib.pyplot as plt
import os
import pandas as pd

pwd = os.path.expanduser('~/2023科研/RSA_Timing/')

df = pd.read_csv(pwd + 'output.csv')

plt.figure(figsize=(10,6))
plt.plot(df['THREADS'], df['TIME_SEQ(s)'], label='Sequential')
plt.plot(df['THREADS'], df['TIME_OMP(s)'], label='OpenMP')
plt.plot(df['THREADS'], df['TIME_THREAD(s)'], label='Thread')
plt.xlabel('Number of Threads')
plt.ylabel('Time (s)')
plt.legend()
plt.title('Time vs. Number of Threads')
plt.grid(True)
plt.savefig(pwd + 'plot_time_vs_threads.png')
plt.show()

# Calculate speedups
df['SPEEDUP_OMP'] = df['TIME_SEQ(s)'] / df['TIME_OMP(s)']
df['SPEEDUP_THREAD'] = df['TIME_SEQ(s)'] / df['TIME_THREAD(s)']

# Determine the sequential time from OpenMP with 1 thread
sequential_time_omp_1thread = df[df['THREADS'] == 1]['TIME_OMP(s)'].values[0]

# Create an ideal speedup line based on this value
df['IDEAL_SPEEDUP'] = sequential_time_omp_1thread / df['TIME_OMP(s)']

# Plotting Speedup vs. Number of Threads
plt.figure(figsize=(10,6))
plt.plot(df['THREADS'], df['SPEEDUP_OMP'], label='Speedup (OpenMP)')
plt.plot(df['THREADS'], df['SPEEDUP_THREAD'], label='Speedup (Thread)')
plt.plot(df['THREADS'], df['THREADS'], linestyle='--', color='grey', label='Ideal Speedup')  # Plotting the ideal speedup line
plt.xlabel('Number of Threads')
plt.ylabel('Speedup')
plt.legend()
plt.title('Speedup vs. Number of Threads')
plt.grid(True)
plt.savefig(pwd + 'plot_speedup.png')
plt.show()