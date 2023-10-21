import matplotlib.pyplot as plt
from monte_carlo import MonteCarlo1D
import numpy as np

# define Monte-Carlo simulation and parameters for Argon
simulator_1d = MonteCarlo1D(epsilon=0.2378, sigma=3.41, r_cut=30)
n_steps = 1000
max_it = 1000
criteria = 0.001 / 100

# run simulation varying number of maximum ierations
result_var = np.array(
    [simulator_1d.calculate_work(n_it, criteria) for n_it in range(10, n_steps + 1)]
)

# run simulation with fixed number of maximum iterations
result_fix = np.array(
    [simulator_1d.calculate_work(max_it, criteria) for _ in range(n_steps)]
)

# initialize plot
fig, ax = plt.subplots(1, 2, figsize=(20, 10))

# plot fisrt simulation results
ax[0].scatter(range(10, n_steps + 1), result_var[:, 1], color="black")
ax[0].plot([0, 1010], [-0.2378, -0.2378], color="red", ls="--", lw=2)
ax[0].set_xlim([0, 1010])
ax[0].set_xlabel("Number of Iterations")
ax[0].set_ylabel(r"Mean work (kcal/mol)")
ax[0].set_title("Varying maximum itetarions")

# plot second simulations results
ax[1].hist(result_fix[np.where(result_fix[:, 0] != 0)][:, 0], ec="black", bins=20)
ax[1].set_xlabel("Number of steps")
ax[1].set_ylabel("Occurance")

# convergence rate of the simulation
conv_rate = sum(n_it < max_it for n_it in result_fix[:, 2]) / n_steps
ax[0].set_title(
    f"Fixed maximum itetarions\nConvergence rate = {np.around(conv_rate, 5)}"
)

plt.show()
