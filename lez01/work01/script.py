import matplotlib
matplotlib.use('Agg') # set the backend before importing pyplot

import matplotlib.pyplot as plt # etc. etc.
import numpy as np
import math

x = np.arange(2, len(open('averages_PI.dat').readlines()) + 2)
y = np.loadtxt('averages_PI.dat') 
error_bars = np.loadtxt('errors_PI.dat')

plt.figure(figsize=(12, 9))
plt.errorbar(x, y, yerr=error_bars, fmt='o', color='goldenrod', ecolor='khaki', markevery=20)
mean_pi=y[-1]
plt.axhline(y=mean_pi, color='darkgoldenrod', label='<$\pi$> = '+str(mean_pi)+'+-'+str(error_bars[-1]), ls='dashed' ) #ls=linestyle
plt.axhline(y=math.pi, color='darkgoldenrod', label='$\pi$ value from Python = '+str(math.pi), ls='dashdot' ) #ls=linestyle
plt.xlabel('#blocks')
plt.ylabel('approximate $\pi$')
plt.grid(True)
plt.title('Buffons$\'$ experiment simulation')
plt.legend()
plt.show()
plt.savefig('plot_buff.png')  # Save the plot to a file
