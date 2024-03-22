import time #Biblioteca para medir el tiempo
import numpy as np
from sympy import symbols, sympify, diff, lambdify
def chebyshev_halley(valor_inicial,funcion, tolerancia, iterMax):
    start4 = time.time() #Inicio del cronómetro
    variable_simbolica = symbols("x")#Definir la variable simbólica 'x'.
    funcion_simbolica = sympify(funcion)#Convertir la cadena de texto funcion en una expresión simbólica.
    funcion_derivada = diff(funcion_simbolica, variable_simbolica)#Calcular la derivada de la función.
    funcion_derivada_2 = diff(funcion_derivada, variable_simbolica)#Calcular la segunda derivada de la función.
    funcion_original_evaluada = lambdify(variable_simbolica, funcion_simbolica, "numpy")#Convertir la función original a una función numérica.
    funcion_derivada_evaluada = lambdify(variable_simbolica, funcion_derivada, "numpy")#Convertir la función derivada a una función numérica.
    funcion_derivada_2_evaluada = lambdify(variable_simbolica, funcion_derivada_2, "numpy")#Convertir la segunda derivada a una función numérica.

    for i in range(iterMax):
        # Fórmula de Chebyshev-Halley
        num = funcion_original_evaluada(valor_inicial) * funcion_derivada_2_evaluada(valor_inicial)
        deno = funcion_derivada_evaluada(valor_inicial)**2
        aprox_cero = valor_inicial - (1 + (1/2)*((num / deno)/1-(num / deno)))*(funcion_original_evaluada(valor_inicial)/funcion_derivada_evaluada(valor_inicial))
        error = max(abs(funcion_original_evaluada(aprox_cero)),abs(aprox_cero - valor_inicial))  # Se calcula el error, con el valor que haga a la función cero.

        # Verificar convergencia
        if error < tolerancia:
            end4 = time.time()
            print("---Método iterativo: chebyshev_halley---")
            print("El tiempo de ejecución fue de:", end4 - start4, "segundos")
            print("Aproximación del cero:", aprox_cero)
            print("Error obtenido:", error)
            print("Cantidad de iteraciones:", i)
            break

        valor_inicial = aprox_cero

chebyshev_halley(2,"416 - 101 * (exp((580 + 182 * x) / 1575) - 1) * ((580 + 182 * x) / 1575) - x", 1e-8, 1000)

