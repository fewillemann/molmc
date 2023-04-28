# exercício desafio 1 - PGF5216
# nome: Felipe Reibnitz Willemann

import numpy

# define pair-wise force from Lennard-Jones potential 
def Forca_LJ(epsilon, sigma, distancia):

	potencia = (sigma/distancia)**6
	return 24*epsilon*potencia*(2*potencia - 1)/distancia

# define physical variables and parameters
epsilon = 0.2378                 # kcal/mol
sigma = 3.41                     # angstrom
r_equilibrio = sigma*(2**(1/6))  # angstrom
r_corte = 30                     # angstrom
delta_r = r_corte - r_equilibrio # angstrom

# define Monte-Carlo parameters
max_passos = 1000000
numero_de_MC = 50000
criterio_convergencia = 0.001/100

# strings used to track necessary MC steps unitl convergence
str_passos, str_trabalhos = '', ''

# initiate MC loop
for i in range (numero_de_MC):

	# initialte total force, steps index and convergence variable
	soma_Forca = 0
	i = 1
	converge = False

	# initiate loop to approximate the work function
	while (i < max_passos) and (not converge):

		# generate random position
		rand = numpy.random.uniform(r_equilibrio, r_corte)
		
		# calculate force sum and the work function
		soma_Forca += Forca_LJ(epsilon, sigma, rand)
		Trab_MC = delta_r * soma_Forca / k

		# check for convergence
		if (abs(abs(Trab_MC) - epsilon) <= criterio_convergencia):
			converge = True
			str_passos += str(k) + '\n'
			str_trabalhos += str(Trab_MC) + '\n'

		else: i += 1

	if (not converge):
		str_passos += str(max_passos) + '\n'

# write results
passos, trabalho = open('saida_passos.txt', 'w'), open('saida_trabalhos.txt', 'w')
passos.write(str_passos)
passos.close()
trabalho.write(str_trabalhos)
trabalho.close()
