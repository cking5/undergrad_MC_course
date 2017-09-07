import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12,6))

def lj(r, epsilon=1.0, sigma=1.0):
    if r > 0.0:
        return 4*epsilon*(sigma**12/r**12-sigma**6/r**6)
    else:
        return None 
    
def flj(r, epsilon=1.0, sigma=1.0):
    if r > 0.0:
        return -4*epsilon*((-12*(sigma)**12)/r**13 + (6*(sigma)**6)/r**7)
    else:
        return None
    
x = np.arange(0.1, 3.0, 0.01)
z = np.arange(0, 2**(1/6), 0.1)
w = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
a = np.arange(-1, 0.1, 0.1)
b = [2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6), 2**(1/6)]

vlj = np.vectorize(lj)
vflj= np.vectorize(flj)

plt.subplot(121)
plt.plot(x, vlj(x), 'b-')
plt.plot(z, w, 'k--')
plt.plot(b, a, 'k--')
plt.title('(a)', fontsize=20, fontweight='bold')
plt.xlabel('$r_{ij}$', fontsize=18, fontweight='bold')
plt.text(0.1, 0.2, r'$\phi_{ij}$', fontsize=18, fontweight='bold')
plt.annotate(r'$\sigma$', xy=(1.0, 0.0), xytext=(1.5, 0.8), fontsize=14,
             fontweight='bold', arrowprops=dict(facecolor='black', shrink=0.05))
plt.text(0.4, -0.9, r'$\epsilon$', fontsize=18, fontweight='bold')
plt.text(1.15, -0.1, '$r_{min}$', fontsize=18, fontweight='bold')
plt.xlim(0.5, 3.0)
plt.ylim(-1.5, 1.0)

ax = fig.add_subplot(121)
ax.plot(x, vlj(x))
ax.spines['left'].set_position(('axes', 0))
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position(('axes', 0.6))
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.subplot(122)
plt.plot(x, vflj(x), 'b-')
plt.title('(b)', fontsize=20, fontweight='bold')
plt.text(0.7, 0.3, '$F_{ij}$', fontsize=18, fontweight='bold')
plt.text(1.57, -0.3, '$r_{ij}$', fontsize=18, fontweight='bold')
plt.annotate('$r_{min}$', xy=(1.12, 0.0), xytext=(1.5, 0.5), fontsize=18,
             fontweight='bold', arrowprops=dict(facecolor='black', shrink=0.05))
plt.xlim(1.0, 3.0)
plt.ylim(-3, 1.0)

ax = fig.add_subplot(122)
ax.plot(x, vflj(x))
ax.spines['left'].set_position(('axes', 0))
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position(('axes', 0.75))
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.show()