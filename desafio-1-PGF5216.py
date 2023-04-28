#Exercício desafio 1 - PGF5216
#Nome: Felipe Reibnitz Willemann
#Numero USP: 10740078

import numpy

def Forca_LJ(epsilon, sigma, distancia):

	potencia = (sigma/distancia)**6
	return 24*epsilon*potencia*(2*potencia - 1)/distancia

epsilon = 0.2378
#em kcal/mol

sigma = 3.41
#em angstron

r_equilibrio = sigma*(2**(1/6))
r_corte = 30
delta_r = r_corte - r_equilibrio
#em angstron

Trab_exato = -epsilon
max_passos = 1000000
numero_de_MC = 50000
criterio_convergencia = 0.001/100
#parametros

str_passos, str_trabalhos = '', ''
#string que ira conter os passos necessarios para convergir

for i in range (numero_de_MC):

	soma_Forca = 0
	k = 1
	converge = False

	while (k < max_passos) and (not converge):

		rand = numpy.random.uniform(r_equilibrio, r_corte)
		#gera numero aleatorio

		soma_Forca += Forca_LJ(epsilon, sigma, rand)
		Trab_MC = delta_r*soma_Forca/k

		if (abs(Trab_exato + abs(Trab_MC)) <= criterio_convergencia):
			converge = True
			str_passos += str(k) + '\n'
			str_trabalhos += str(Trab_MC) + '\n'
			#salva quantos passos demorou para convergir e recomeca o MC

		else: k += 1

	if (not converge):
		str_passos += str(max_passos) + '\n'

passos, trabalho = open('saida_passos.txt', 'w'), open('saida_trabalhos.txt', 'w')
passos.write(str_passos)
passos.close()
trabalho.write(str_trabalhos)
trabalho.close()
#guarda os passos em um arquivo de texto