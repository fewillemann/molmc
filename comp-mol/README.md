# Monte Carlo integration methods

These codes were written to estimate the energy of two model systems: 
1. Interaction potential of the Argon dimer (one dimensional)
2. Mean interaction energy of liquid Argon (three dimensional)

> The pair-wise potential was approximated as the **Lennard-Jones potential** in both cases.

## The Argon dimer

Using random number generation, the work **W** required to separate two argon atoms (from the equilibrium position to an "infinitely large" distance of $\text{30\AA}$, where the potential would be zero) was estimated through the Monte Carlo method of one-dimensional integration, in which a probability density is used to calculate the average of the integrand in a uniform distribution within the interval of integration in question. Therefore, we have:

$$W = \int_{x_{1}}^{x_{2}} F(r) \ dr = (x_{2} - x_{1}) \langle F(\xi) \rangle_{\tau} = (x_{2} - x_{1}) \frac{1}{\tau} \sum_{i}^{\tau} F(\xi_{i})$$

Where $\xi$ is a random number between $x_1$ (**equilibrium position**) and $x_2$ (**"infinity"**). Then, the Lennard-Jones force is:

$$U(r) = 4\varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right) ^6  \right] \rightarrow  F(r) = -\frac{\partial U(r)}{\partial r} = 24\varepsilon \left[ \left( \frac{\sigma^{12}}{r^{13}} \right) - \left( \frac{\sigma^6}{r^7} \right)  \right]$$

Where $\varepsilon$ is an energy characteristic parameter and $\sigma$ is a distance characteristic parameter. For argon, $\varepsilon = 0.2378 \text{ kcal/mol}$ and $\sigma = 3.41 \text{ \AA}$. Using the above equations, the work was calculated with a precision of 0.001\% compared to the exact work $W_{exact} = -\varepsilon$.

## Liquid Argon

Through the Monte Carlo method, the average interaction energy in the argon liquid was estimated. The density of this liquid is typically 1.394 g/cm$^3$, which is equivalent to 100 atoms in a box of \text{L = 16.8192 \AA} side. For this, the x, y and z coordinates of each atom were generated randomly (using a uniform distribution), creating a configuration in which the potential between each atom is given by the Lennard-Jones potential U(r) of the equation \ref{potencial_LJ} . Thus, the average energy will be the sum of all potentials divided by the number of atoms plus a kinetic part.

$$E_{media} = \frac{U_{total}}{N} + \frac{3}{2}k_{B}T$$

Where $k_{B}$ is the Boltzmann constant, T is the temperature, and the total potential is given by:

$$U_{total} = \sum_{1}^{N - 1} \sum_{j > i}^{N} U_{ij}(r_{ij})$$

As we want to simulate a real liquid, we use the **image method** (periodic boundary conditions) to correct the potential of very distant atoms, which consists of reducing the separation distance by half of L when it is greater than half of L, thus considering that As one atom moves further away from another (as if it were leaving the box), another closer atom will be used in calculating the potential (as if another atom were entering the box, keeping their total number constant).

The objective is to verify if the result will be close to that expected by the Boltzmann distribution for liquid argon, that is, if the average energy is between -2 kcal/mol and 0 kcal/mol, which is what is observed experimentally.
