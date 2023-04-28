# exercicio desafio 2 - PGF5216
# nome: Felipe Reibnitz Willemann

import numpy

# define Lennard-Jones potential
def Potencial_LJ(epsilon, sigma, dist_quad):

	potencia = dist_quad*dist_quad*dist_quad/sigma/sigma/sigma
	return 4*epsilon*potencia*(potencia - 1)

# define physical variables and constants
epsilon = 0.2378                   # kcal/mol
sigma = 3.41                       # angstrom
L = 16.8192                        # angstrom
T = 90                             # K
k_boltzmann = 0.001987             # kcal/mol/K
termo_cinetico = 3*k_boltzmann*T/2 # kcal/mol

# define monte-carlo integration parameters
str_energias = ''        # string  mean energy
passos_MC = 1000000
numero_particulas = 100
prox_boltzmann = 0

# initiate monte-carlo iteration
for k in range(passos_MC):
	# create the gas box configuration (coordinate matrix) 
	coordenadas = [[numpy.random.uniform(-L/2, L/2), numpy.random.uniform(-L/2, L/2), numpy.random.uniform(-L/2, L/2)] for i in range(numero_particulas)]

	# initialte total potential energy
	potencial_acumulado = 0

	# run through atoms pairs
	for i in range(numero_particulas - 1):
		for j in range(i + 1, numero_particulas):
			# calculate distances
			distancia_x = abs(coordenadas[i][0] - coordenadas[j][0])
			distancia_y = abs(coordenadas[i][1] - coordenadas[j][1])
			distancia_z = abs(coordenadas[i][2] - coordenadas[j][2])
			
			# calculate effective distances (considering periodic boundary conditions)
			dist_efetiva_x = distancia_x - round(distancia_x/L)*L/2
			dist_efetiva_y = distancia_y - round(distancia_y/L)*L/2
			dist_efetiva_z = distancia_z - round(distancia_z/L)*L/2

			# calculate total distance
			distancia_quadrado = dist_efetiva_x**2 + dist_efetiva_y**2 + dist_efetiva_z**2
			
			# add potential energy of the pair to the total
			potencial_acumulado += Potencial_LJ(epsilon, sigma, distancia_quadrado)
				
	# calculate mean energy and add to the ensamble
	energia_media = potencial_acumulado/N + termo_cinetico
	str_energias += str(energia_media) + '\n'

	# check if energy corresponds to the boltmann distributions prediction
	if energia_media > -2 and energia_media < 0:
		prox_boltzmann += 1

# write results
porcentagem = 100*prox_boltzmann/passos_MC
energias = open('saida_energias.txt', 'w')
energias.write(str_energias)
energias.write(str(porcentagem))
energias.close()
