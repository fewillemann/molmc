import matplotlib.pyplot as plt
import numpy as np
from monte_carlo import MonteCarlo3D

# define Monte-Carlo parameters
n_steps = 100_000
n_particles = 100

# run simulation
simul = MonteCarlo3D(box_lenght=16.819, temperature=90)
result = [simul.mc_step(n_particles) for _ in range(n_steps)]
energy = [step[1] for step in result]

# write results
print("rate of acceptance:", sum([step[0] for step in result]) / n_steps)

fig, ax = plt.subplots(figsize=(20, 10))
log_bins = np.logspace(np.log10(min(energy)), np.log10(max(energy)), 20)
ax.hist(energy, bins=[-0.2, 0] + list(log_bins), ec="black", color="tab:blue")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Mean energy (kcal/mol)")
ax.set_ylabel("Occurance")
ax.set_title("Argon box mean energy per particle (log scale)")

plt.show()
