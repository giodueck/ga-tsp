# ga-tsp
TSP resuelto usando un algoritmo genetico para el curso de Sistemas Paralelos y Distribuidos

# Uso
El script `compile.sh` compila el codigo, para pasar parametros como `-g` o `-O3` se pueden pasar como parametros al script, ej. `./compile.sh -O3`, luego `./ga-tsp -h` para ayuda de uso. Para pasar multiples opciones usar la notacion de comillas: `./compile.sh '-O3 -Werror'`

El codigo y el programa estan en ingles, ya que suelo programar e investigar en ingles.

# OpenMP y Pthreads
Estos dos metodos son los que se usan para la paralelizacion, y la seleccion se hace al compilar con o sin la bandera de compilador `-fopenmp`.

El programa tiene el mismo resultado bajo las dos metodologias, pero usar OpenMP es bastante mas simple que las funciones de pthreads. En cambio, el programa es ligeramente mas rapido cuando es compilado con pthreads.

# Estrategia de paralelizacion
El algoritmo genetico entero se ejecuta en varias instancias semi-independientes, esto se llama el modelo de islas. Cada cierto numero de generaciones, las poblaciones de las islas son cruzadas para intercambiar estrategias efectivas y mantener una buena diversidad genetica.

# Features
El algoritmo genetico tiene parametros que pueden ser especificados al ejecutar el programa, y funciona con archivos [TSPLIB](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/).

Los parametros son:
- Numero de generaciones
- Tamaño de poblacion
- Probabilidad de mutacion
- Tamaño de torneo
- Porcentaje de elitismo
- Cantidad de islas paralelas (poblacion dividida entre las islas)
- Frecuencia de cruce entre islas
- Frecuencia de impresion de estadisticas en la consola
- Seed para el PRNG

# Creditos
TSPLIB es un proyecto de la universidad de Heidelberg con una libreria de problemas de prueba y una documentacion del tipo de archivo que los representa.

Problemas incluidos son de la universidad de Waterloo en su pagina [TSP Test Data](https://www.math.uwaterloo.ca/tsp/data/).
