"""This module defines two molecular Monte Carlo algorithms. Written as a programming exercise.
Author: Felipe Reibnitz Willemann"""
from typing import Tuple

import numpy


class MonteCarlo1D:
    def __init__(self, epsilon: float, sigma: float, r_cut: float) -> None:
        """Initialize.

        Args:
            epsilon (float): Lennard-Jones energy at equilibrium in kcal/mol.
            sigma (float); Lennard-Jones radius parameter in Angstrom.
            r_cut (float): equilibrium radius in Angstrom.
        """
        # interaction parameters
        self.epsilon = epsilon  # 0.2378 kcal/mol
        self.sigma = sigma  # 3.41 angstrom
        self.r_cut = r_cut
        self.r_eq = self.sigma * (2 ** (1 / 6))

        # random number generator
        self.random = numpy.random.uniform

    # ----------------------------------------------------------------------------------
    # public methods
    def calculate_work(self, max_it: int, criteria: float) -> Tuple[int, float]:
        """Define Monte Carlo work calculation.

        Args:
            max_it (int): maximum number of Monte Carlo iterations.
            criteria (float): energy convergence criteria.

        Returns:
            n_it, work (Tuple[int, float]): number of iterations calculated and total work in kcal/mol.
        """
        total_force = 0.0

        for n_it in range(1, max_it + 1):
            # calculate force sum and the work function
            total_force += self._lj_force(self.random(self.r_eq, self.r_cut))
            work = (self.r_cut - self.r_eq) * total_force / n_it

            # check for convergence
            if abs(work / self.epsilon + 1) < criteria:
                break

        return n_it, work  # / self.epsilon + 1

    # ----------------------------------------------------------------------------------
    # private methods
    def _lj_force(self, distance: float) -> float:
        """Define pair-wise force from Lennard-Jones potential.

        Args:
            distance (float): separation distance between atoms in Angstrom.

        Returns:
            value of the force (float).
        """
        par = (self.sigma / distance) ** 6

        return 24 * self.epsilon * par * (2 * par - 1) / distance


class MonteCarlo3D:
    def __init__(
        self, box_lenght: float, temp: float, epsilon: float, sigma: float
    ) -> None:
        """Initialize.

        Args:
            box_lenght (float): lenght of the simulation box in Angstroms.
            temp (float); temperature of the system in Kelvin.
            epsilon (float): Lennard-Jones energy at equilibrium in kcal/mol.
            sigma (float); Lennard-Jones radius parameter in Angstrom.
        """
        # physical and interaction parametrs
        self.L = box_lenght
        self.T = temp
        self.epsilon = epsilon
        self.sigma = sigma

        k_boltzmann = 0.001987  # kcal/mol/K
        self.kinetic = 3 * k_boltzmann * self.T / 2  # kcal/mol

        # random number generator and L2 norm
        self.random = numpy.random.uniform
        self.norm = numpy.linalg.norm

    # ----------------------------------------------------------------------------------
    # public methods
    def calculate_energy(self, n_particles: int) -> Tuple[bool, float]:
        """Define Monte Carlo mean energy calculation.

        Args:
            n_particles (int): number of particles in the system.

        Returns:
            mean energy (float): mean energy per particle of the system.
        """
        # create the gas box configuration (coordinate matrix)
        coords = self.random(-self.L / 2, self.L / 2, (n_particles, 3))
        total_potential = 0.0

        # iterate over atoms pairs
        for i in range(n_particles - 1):
            for j in range(i + 1, n_particles):
                # calculate effective distances considering periodic boundary conditions
                dist = coords[i] - coords[j]
                efective_dist = dist - numpy.round(dist * 2 / self.L) * self.L / 2

                # add to acumulated potential
                total_potential += self._lj_potential(self.norm(efective_dist))

        return total_potential / n_particles + self.kinetic

    # ----------------------------------------------------------------------------------
    # private methods
    def _lj_potential(self, distance: float) -> float:
        """Define pair-wise force from Lennard-Jones potential.

        Args:
            distance (float): separation distance between atoms in Angstrom.

        Returns:
            value of the potential (float).
        """
        par = (self.sigma / distance) ** 6

        return 4 * self.epsilon * par * (par - 1)
