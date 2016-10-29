import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,10,0.1)

# The lines to plot
y1 = 4 - 2*x
y2 = 3 - 0.5*x
y3 = 1 -x

# The upper edge of polygon (min of lines y1 & y2)
y4 = np.minimum(y1, y2)

# Set y-limit, making neg y-values not show in plot
plt.ylim(0, 5)

# Plotting of lines
plt.plot(x, y1,
         x, y2,
         x, y3)

# Filling between line y3 and line y4
plt.fill_between(x, y3, y4, color='grey', alpha='0.5')
plt.show()
