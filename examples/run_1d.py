import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src/molmc")
#sys.path.append("../src/molmc")

from molmc import MonteCarlo1D
import numpy as np

# define Monte-Carlo simulation and parameters for Argon
simulator = MonteCarlo1D(epsilon=0.2378, sigma=3.41, r_cut=30)
n_steps = 500
max_it = 500
criteria = 0.001 / 100

# run simulation varying number of maximum ierations
result_var = np.array(
    [simulator.calculate_work(n_it, criteria) for n_it in range(10, n_steps + 1)]
)

# run simulation with fixed number of maximum iterations
result_fix = np.array(
    [simulator.calculate_work(max_it, criteria) for _ in range(n_steps)]
)

# initialize plot
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# plot fisrt simulation results
ax[0].scatter(
    range(10, n_steps + 1),
    result_var[:, 1],
    color="black",
    s=10,
    label="Calculated work",
)
ax[0].plot(
    [0, n_steps],
    [-simulator.epsilon, -simulator.epsilon],
    color="red",
    ls="--",
    lw=2,
    label="Exact work",
)
ax[0].set_xlim([0, n_steps * 1.01])
ax[0].set_xlabel("Number of Iterations")
ax[0].set_ylabel(r"Mean work (kcal/mol)")
ax[0].set_title("Varying maximum iterations")
ax[0].legend()

# plot second simulations results
ax[1].hist(result_fix[np.where(result_fix[:, 0] != max_it)][:, 0], ec="black", bins=20)
ax[1].set_xlabel("Number of iterations")
ax[1].set_ylabel("Occurance")

# convergence rate of the simulation
conv_rate = sum(n_it < max_it for n_it in result_fix[:, 0]) / n_steps
ax[1].set_title(
    f"Fixed maximum iterations (< {max_it})\nConvergence rate = {np.around(conv_rate, 5)}"
)

plt.tight_layout()
plt.show()
# plt.savefig("result_1d.png", dpi=300)
