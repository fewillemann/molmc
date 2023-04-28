#Exercicio desafio 2 - PGF5216
#Nome: Felipe Reibnitz Willemann
#Numero USP: 10740078

import numpy

def Potencial_LJ(epsilon, sigma, dist_quad):

	potencia = dist_quad*dist_quad*dist_quad/sigma/sigma/sigma
	return 4*epsilon*potencia*(potencia - 1)

L = 16.8192
sigma = 3.41
#em angstron

T = 90
#em kelvin

k_boltzmann = 0.001985
termo_cinetico = 3*k_boltzmann*T/2
epsilon = 0.2378
#em kcal/mol

str_energias = ''
passos_MC = 1000000
numero_particulas = 100
prox_boltzmann = 0
#parametros

for k in range(passos_MC):

	coordenadas = [[numpy.random.uniform(-L/2, L/2), numpy.random.uniform(-L/2, L/2), numpy.random.uniform(-L/2, L/2)] for i in range(numero_particulas)]
	#gera a caixa contendo os atomos de argonio

	N = len(coordenadas)
	potencial_acumulado = 0

	for i in range(N - 1):
		for j in range(i + 1, N):
			distancia_x = abs(coordenadas[i][0] - coordenadas[j][0])
			distancia_y = abs(coordenadas[i][1] - coordenadas[j][1])
			distancia_z = abs(coordenadas[i][2] - coordenadas[j][2])
			dist_efetiva_x = distancia_x - round(distancia_x/L)*L/2
			dist_efetiva_y = distancia_y - round(distancia_y/L)*L/2
			dist_efetiva_z = distancia_z - round(distancia_z/L)*L/2
			#metodo das imagens

			distancia_quadrado = dist_efetiva_x**2 + dist_efetiva_y**2 + dist_efetiva_z**2
			potencial_acumulado += Potencial_LJ(epsilon, sigma, distancia_quadrado)
			#soma o potencial do atomo i com todos os outros atomos
				
	energia_media = potencial_acumulado/N + termo_cinetico
	str_energias += str(energia_media) + '\n'

	if energia_media > -2 and energia_media < 0:
		prox_boltzmann += 1

porcentagem = 100*prox_boltzmann/passos_MC
energias = open('saida_energias.txt', 'w')
energias.write(str_energias)
energias.write(str(porcentagem))
energias.close()