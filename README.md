# Monte-Carlo-Simulations-For-Path-Integral-Calculations

En este proyecto se calculan, mediante un algoritmo de Metropolis (MCMC), las llamadas integrales de camino de Feynman. Matemáticamente, estas son integrales funcionales; es decir, integrales sobre un espacio de funciones en lugar de espacios finito-dimensionales como $\mathbb{R}^N$ o $\mathbb{C}^N$.

El archivo `Proyecto_Integral_de_Caminos.ipynb` contiene la información matemática y física necesaria para una comprensión profunda del método. En este documento se desarrolla principalmente la arquitectura algorítmica empleada para lograr este cometido.

## Tecnologías Usadas
* **Python**: Lenguaje principal de simulación.
* **NumPy**: Para el manejo eficiente de arrays, trayectorias discretizadas y generación de números aleatorios.
* **Matplotlib / FuncAnimation**: Para la visualización de resultados, densidades de probabilidad y renderizado de animaciones de las trayectorias de Markov.

---

## Fundamento Físico y Matemático

En física cuántica, podemos expresar un propagador $\langle x_f, t_f \mid x_i, t_i \rangle$ de la siguiente manera:

$$\langle x_f, t_f \mid x_i, t_i \rangle = \int_{x_i}^{x_f} \mathcal{D}[x(t)] \exp \left( \frac{i}{\hbar} \int_{x_i}^{x_f} dt \, L_{\text{clásico}}(x,\dot{x}) \right)$$

donde la medida funcional $\mathcal{D}[x(t)]$ se define como:

$$\mathcal{D}[x(t)] = \lim_{(n, \epsilon) \to (\infty, 0)} \left( \frac{m}{2 \pi i \hbar \epsilon} \right)^{\frac{n}{2}} \int dx_{n-1} \int dx_{n-2} \dots \int dx_1$$

y $L_{\text{clásico}}$ es el lagrangiano del sistema, definido como:

$$L_{\text{clásico}} = \frac{1}{2}mv^2 - V(x),$$

siendo $m$ la masa de la partícula, $v$ su velocidad y $V(x)$ el potencial al que está sometida.

Físicamente, este enfoque es valioso porque permite abordar problemas complejos de manera sistemática. Sin embargo, dado que estas integrales solo tienen soluciones analíticas exactas para casos muy específicos, es indispensable recurrir a métodos de aproximación numérica. (En el *notebook* adjunto, se demuestra cómo una aproximación ingenua sobre $\mathbb{R}^N$ resulta computacionalmente ineficiente).

## Implementación del Algoritmo (Metropolis MCMC)

Este programa calcula las integrales de camino en **tiempo imaginario** ($t \to -i \tau$). Esta transformación (rotación de Wick) convierte la integral oscilatoria cuántica en una integral con un peso exponencial real, equivalente a un ensamble estadístico de la mecánica clásica. Este método es ideal para obtener el estado base $E_0$.

Al aplicar la transformación a tiempo imaginario, la acción se vuelve euclidiana y la integral se puede aproximar usando Monte Carlo:

$$\int_{x_i}^{x_f} \mathcal{D}[x(\tau)] \exp \left( - \frac{S_E}{\hbar} \right) \approx \frac{1}{M} \sum_{j=0}^{M} \exp\left(- \sum_{i=0}^{n-1} \left[ \frac{m}{2} \left( \frac{y_{i+1}^j - y_i^j}{\Delta \tau} \right)^2 + V(y_i^j ) \right] \Delta \tau \right)$$

donde $n-1$ es el número de discretizaciones temporales y $M$ es el número de trayectorias aleatorias generadas.

Para evitar el gasto computacional de trayectorias aleatorias con acción altísima (que no aportan físicamente), aplicamos el **Algoritmo de Metropolis**. Este genera una cadena de Markov que muestrea preferentemente las trayectorias físicamente relevantes. 

Si proponemos un cambio local a la trayectoria (un paseo gaussiano simétrico), la probabilidad de aceptar la nueva configuración $X$ desde la actual $Y$ es:

$$A(X \mid Y) = \min\left(1, \exp\left(-(S_{\text{nuevo}} - S_{\text{viejo}})\right) \right)$$

El algoritmo genera un número aleatorio $u \in [0, 1)$ y se acepta la nueva trayectoria si $u < A$.

### Estructura del Código

El script `algoritmoMCMC.py` consta de tres funciones principales:

1. **`accion(path, dtau, potencial, m)`**: Calcula la acción euclidiana discreta $S_E$ de una trayectoria `path` dada.
2. **`metropolis(xi, xf, n, pasos_m, m, delta, potencial, dtau)`**: Partiendo de una trayectoria inicial, selecciona un punto intermedio aleatorio $j$ y le aplica un desplazamiento gaussiano $y_j = x_j + \mathcal{N}(0, \delta)$. Se evalúa la condición de aceptación $u < \exp(-(S_{\text{nuevo}} - S_{\text{viejo}}))$ para reemplazar el punto. Este "barrido" se repite `pasos_m` veces.
3. **`propagador_MCMC(m, xi, xf, tau, n, M, pasos_m, delta, potencial)`**: Orquesta la simulación general. Genera $M$ trayectorias termalizadas independientes mediante la función anterior y las promedia para aproximar la integral funcional.

---

## Casos de Estudio y Resultados

Para validar la robustez del algoritmo, se evaluaron tres sistemas físicos distintos, logrando reproducir los comportamientos esperados para el estado base:

1. **Oscilador Armónico Cuántico:** Dado 

$$
V(x) = \frac{1}{2}m v^2 + \frac{1}{2}\omega^2 x^2
$$

Para $m=1, x_i=0, x_f=0, \tau \in (0.2, 2), n=500, M=500, \text{pasosm}=30, \delta=0.3$ se obtuvo lo siguiente:

> [!IMPORTANT]
> **Resultados del Estado Base ($E_0$)**
> - **$E_0$ numérico:** 0.4834
> - **$E_0$ exacto:** 0.5 
> - **Error relativo:** 3.32%

Esto se puede replicar ejecutando el archivo [oscilador_armonico.py](oscilador_armonico.py).

2. **Pozo de Potencial Infinito:** Dado 

$$
V(x) = 
\begin{cases}
0 \ \text{si} \ |x| < 1 \\
\infty \ \text{en otro caso}
\end{cases}
$$

Para $m=1, x_i=0.5, x_f=0.5, \tau \in (0.5, 2), n=10000, M=10000, \text{pasosm}=500, \delta=0.1$ se obtuvo lo siguiente:

> [!IMPORTANT]
> **Resultados del Estado Base ($E_0$)**
> - **$E_0$ numérico:** 4.7284
> - **$E_0$ exacto:** 4.9348
> - **Error relativo:** 4.18%

Esto se puede replicar ejecutando el archivo [pozo_infinito.py](pozo_infinito.py).

3. **Potencial 'Calabera':** Dado

$$
V(x) = 10 \ \sin ( \pi x)^2 + x^4
$$

Para $m=1, x_i=0.5, x_f=0.5, \tau \in (0.5, 2), n=5000, M=5000, \text{pasosm}=200, \delta=0.1$ se obtuvo lo siguiente:

> [!IMPORTANT]
> **Resultados del Estado Base ($E_0$)**
> - **$E_0$ numérico por MCMC:** 12.7694
> - **$E_0$ numérico por diferencias finitas:** 12.3592
> - **Error relativo:** 3.32%

Esto se puede replicar ejecutando el archivo [pozo_cal.py](pozo_cal.py).

### Animación de trayectorias
A continuación, se muestran las trayectorias generadas para el caso del oscilador armónico:

[Animación Metropolis](animacion_trayectorias.gif)

## Estructura del Repositorio y Uso

Este proyecto ofrece dos formas de explorar las simulaciones, dependiendo del nivel de detalle que busques:

1. **El Notebook Interactivo (`Proyecto_Integral_de_Caminos.ipynb`)**: 
   Es el documento principal y el punto de partida recomendado. Contiene toda la fundamentación teórica (matemática y física), la construcción paso a paso del algoritmo de Metropolis y la ejecución detallada de todos los casos de estudio (oscilador armónico, pozo infinito, etc.) junto con sus visualizaciones. Ideal para entender el "porqué" del proyecto.

2. **Los Scripts de Python (`.py`)**: 
   Si prefieres evaluar el código fuente o ejecutar las simulaciones directamente desde la terminal, el repositorio incluye scripts independientes para cada caso de estudio. Cada archivo `.py` es autocontenido: incluye las funciones del algoritmo MCMC y la configuración específica del potencial físico a simular. 

**Para ejecutar localmente:**
Solo necesitas tener instalados `numpy` y `matplotlib`. Puedes clonar el repositorio y correr cualquiera de los scripts individuales en tu terminal, o abrir el *notebook* para correr las celdas de forma interactiva.
