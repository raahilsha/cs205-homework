import matplotlib.pyplot as plt
import math
import numpy as np

x = range(1, 512)
y1 = [i - 1 for i in x]
y2 = [math.ceil(math.log(i, 2)) for i in x]

plt.plot(x, y1, '.b', label='Lone Cashier')
plt.plot(x, y2, '.g', label='Infinite Employees')
plt.yscale('log')
plt.xlabel('Number of Bags')
plt.ylabel('Time Taken (s)')
plt.title('Problem 5 - Part 4')
plt.show()