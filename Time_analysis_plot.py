import matplotlib.pyplot as plt
import numpy as np

# 读取文件
with open('/home/siwtan/2023科研/RSA_Timing/Timing.txt', 'r') as f:
    times = [int(line.strip()) for line in f]

# 计算平均时间
avg_time = np.mean(times)

# 找出超过平均时间的位
key_bits = [i for i, t in enumerate(times) if t > avg_time]

# 打印可能的私钥位
print('Possible key bits:', key_bits)

# 绘制时间图
plt.plot(times, label='Execution times')
plt.axhline(y=avg_time, color='r', linestyle='-', label='Average time')
plt.legend()
plt.show()
