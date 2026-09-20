def accion(path, dtau, potencial, m): #Aproximacion numerica de la accion 
    a = np.diff(path)**2
    int_K = 0.5*m*np.sum(a)/dtau   #Aproximacion de la integral de la energia cinetica
    int_potencial = np.sum([potencial(punto, m) for punto in path[:-1]])*dtau  #Aproximacion de la integral del potencial

    return int_potencial + int_K

def metropolis(xi, xf, n, pasos_m, m, delta, potencial, dtau):
    path = np.linspace(xi, xf, n+1)

    for _ in range(pasos_m):
        i = np.random.randint(1, n)  #Elige un numero aleatorio entre 1 y n
        
        y = path[i]   #Que sea i y no el numero en el rango del for hace mas estocastico el metodo 
        x = y + np.random.normal(0, delta) # x + delta*epsilon

        x_prev, x_next = path[i-1], path[i+1] 

        dy1, dy2 = y - x_next, y - x_prev
        dx1, dx2 = x - x_next, x - x_prev

        #Ahora calculamos la accion alrededor de estos puntos
        S_y = 0.5*m*(dy1**2 + dy2**2)/dtau + potencial(y, m)*dtau
        S_x = 0.5*m*(dx1**2 + dx2**2)/dtau + potencial(x, m)*dtau
        #Aqui se esta calculando la accion de una manera diferente que en la funcion por eso no se utiliza

        if np.random.rand() < np.exp(S_y - S_x) : # Aqui se decide si se cambia y -> x con la condicion de metropolis
            path[i] = x
        
    return path

def propagador_MCMC(m, xi, xf, tau, n, M, pasos_m, delta, potencial):
    '''
    m : masa 
    xi : posicion inicial 
    xf : posicion final 
    tau : diferencia entre los tiempos imganarios 
    n : numeros de puntos en una trayectoria
    M : numero de trayectorias 
    pasos_m : numero de pasos en metropolis
    delta : Se usa para modificar un punto xi a otro xi' = xi + delta para generar la trayectoria de metropolis
    '''
    dtau = tau/n    
    weights = [] #Lista en la que se pondran los sumandos de la aproximacion de Monte Carlo
    paths = []
    for _ in range(M): 
        path = metropolis(xi, xf, n, pasos_m, m, delta, potencial, dtau)  #Hacemos una trayectoria aleatoria de metropolis
        S_val = accion(path, dtau, potencial, m)
        paths.append(path)
        
        if np.isfinite(S_val):   #Aqui verifica que la accion sea finita, sirve para los casos donde hay inf involucrados
            weights.append(np.exp(-S_val))
        
    if weights:    
        K_tau = np.mean(weights) #Esta ya es la aproximacion de Monte Carlo 
        return K_tau, paths
    print('Error: construccion invalida')
    return np.nan
