import matplotlib.pyplot as plt
from molmc import MonteCarlo3D
import numpy as np

# define Monte-Carlo parameters
n_steps = 1000
n_particles = 100

# define simulation for liquid Argon
simulator = MonteCarlo3D(box_lenght=16.819, temp=90, epsilon=0.2378, sigma=3.41)

# run simulator
energy = np.array([simulator.calculate_energy(n_particles) for _ in range(n_steps)])

# write results
conv = sum(((energy > -2) & (energy < 0)) / n_steps)

# plot results
fig, ax = plt.subplots(
    figsize=(14, 6))
log_bins = np.logspace(np.log10(min(energy)), np.log10(max(energy)), 20)
ax.hist(energy, bins=[-0.2, 0] + list(log_bins), ec="black", color="tab:blue")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Mean energy (kcal/mol)")
ax.set_ylabel("Occurance")
ax.set_title(f"Argon box mean energy per particle (log scale)\nrate of acceptance = {conv}")

plt.tight_layout()
plt.show()
