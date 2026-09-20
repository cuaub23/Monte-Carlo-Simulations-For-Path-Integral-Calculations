# Monte-Carlo-Simulations-For-Path-Integral-Calculations
En este proyecto se calcula mediante un algoritmo metropolis las llamadas integrales de caminos de Feynman que matematicamente son integrales funcionales, es decir, integrales sobre un espacio de fuciones en lugar de espacios como $\mathbb{R}^N$ o $\mathbb{C}^N$.

El notebook {notebook}() cuenta con la informacion matematica y fisica necesaria para el entendimiento. Por otro lado, aqui solo se desarrolla la idea algoritmica para lograr este cometido.

En fisica cuantica podemos expresar un propagador $\langle x_f, t_f | x_i, t_i \rangle$ de la siguiente manera 

$$
\braket{x_f,t_f| x_i,t_i} = \int_{x_i}^{x_f} \mathcal{D}[x(t)] \exp \left( \frac{i}{\hbar} \int_{x_i}^{x_f} dt \, L_{clásico}(x,\dot{x}) \right) 
$$

donde se define la medida funcional $\mathcal{D}[x(t)]$:

$$
\mathcal{D}[x(t)] = \lim_{(n, \epsilon) \to (\infty, 0)} \left( \frac{m}{2 \pi i \hbar \epsilon} \right)^{\frac{n}{2}} \int dx_{n-1} \int dx_{n-2} ... \int dx_1 
$$

y $L_{clásico}$ se define como 

$$
L_{clásico} = \frac{1}{2}mv^2 - V(x),
$$

siendo $m$ la masa de la particula, $v$ su velocidad y $V(x)$ una funcion suave que decae en el infinito. 

## Importancia del proyecto
Fisicamente es importante ya que nos permite abordar problemas de cierta indole de manera mas sencilla. Sin embargo, solo se ha calculado analiticamente para ciertos casos por lo que es necesario abordar la aproximacion numerica.

Dentro del notebook se aproxima primeramente como integrales sobre $\mathbb{R}^N$ cuando $N \to \infty$ pero esto no resulta ya que son demasiados calculos para una mala aproximacion como se muestra ahi. 

## Algoritmo 
Primeramente es de mencionar que este programa calcula las integrales de camino en tiempos imaginarios $t \to -i \tau$. Este caso sigue siendo relevante ya que nos permite obtener el espectro de energias en particular es util para obtener el estado base $E_0$. 
Para el calculo de estas integrale haremos uso de las aproximaciones Monte Carlo. Sea $(X_i)_{i=1}^N$ un conjunto de variables aleatorias, independientes e idénticamente distribuidas, entonces para una función $g(\mathbb{X})$ se cumple  

$$  
\begin{align}
G_N = \int g(\mathbb{X}) \, d\mathbb{X} \\
G_N \approx \frac{1}{N} \sum_{i=1}^{N} g(X_i)  
\end{align}
$$  

Esto nos permite aproximar la integral de caminos de la siguiente manera:

$$  
 \int_{x_i}^{x_f}  \mathcal{D}[x(t)] \exp \left( \frac{i}{\hbar} \int_{x_i}^{x_f} dt \, L_{clásico}(x,\dot{x}) \right)  \approx \frac{1}{M} \sum_{j=0}^{M} \exp\left(- \sum_{i=0}^{n-1} \left[ \frac{m}{2} \left( \frac{y_{i+1}^j - y_i^j}{\Delta \tau} \right)^2 + V(y_i^j ) \right] \Delta \tau \right)  
$$  

En donde $n-1$ es el número de variables y $M$ es el número de trayectorias aleatorias, independientes e idénticamente distribuidas.
