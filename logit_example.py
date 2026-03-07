import numpy as np
import matplotlib.pyplot as plt
# vibe code slop
x = np.linspace(0.002,0.998,1000)
logit = np.log(x/(1-x))

f = logit
# f = f/np.trapezoid(f,x)

plt.plot(x,f)
plt.title("PDF proportional to |logit(x)|")
plt.xlabel("x")
plt.ylabel("density")
plt.show()