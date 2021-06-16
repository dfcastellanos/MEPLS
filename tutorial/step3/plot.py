import numpy as np
import matplotlib.pyplot as plt
import sys

# read the data
data = np.genfromtxt( sys.argv[1] ).T
strain = data[0]
stress = data[1]

# strain in % and stress in MPa
strain *= 100
stress *= 1000

plt.plot(strain, stress, linewidth=1.)

# tune plot details
plt.xlabel(r'$\varepsilon_{\rm xy}$ (%)', fontsize=15)
plt.ylabel(r'$\Sigma_{\rm xy}$ (MPa)', fontsize=15)
plt.tick_params(labelsize=13)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.tight_layout()  
plt.show()