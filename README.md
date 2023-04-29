# ga-tsp
TSP resuelto usando un algoritmo genetico para el curso de Sistemas Paralelos y Distribuidos

# Uso
El script `compile.sh` compila el codigo, para pasar parametros como `-g` o `-O3` se pueden pasar como parametros al script, ej. `./compile.sh -O3`, luego `./ga-tsp -h` para ayuda de uso.

El codigo y el programa estan en ingles, ya que suelo programar e investigar en ingles.

# Estrategia de paralelizacion
El algoritmo genetico entero se ejecuta en varias instancias semi-independientes, esto se llama el modelo de islas. Cada cierto numero de generaciones, las poblaciones de las islas son cruzadas para intercambiar estrategias efectivas y mantener una buena diversidad genetica.

# Features
El algoritmo genetico tiene parametros que pueden ser especificados al ejecutar el programa, y funciona con archivos [TSPLIB](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/).

Los parametros son:
- Numero de generaciones
- Tamaño de poblacion
- Probabilidad de mutacion
- Estrategia de seleccion (torneo o truncaciion)
- Tamaño de torneo
- Porcentaje de truncacion
- Porcentaje de cruce (para truncacion)
- Porcentaje de elitismo (para truncacion)
- Cantidad de islas paralelas (poblacion dividida entre las islas)
- Frecuencia de cruce entre islas
- Frecuencia de impresion de estadisticas en la consola
- Seed para el PRNG
