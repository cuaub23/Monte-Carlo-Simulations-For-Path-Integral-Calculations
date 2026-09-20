from algoritmoMCMC import accion, metropolis, propagador_MCMC
def pozo_inf(x, m, L=1.0):
    if x <= 0.0 or x >= L:
        return np.inf 
    return 0.0

tau_values = np.linspace(0.5, 2.0, 10)
log_K = []
for tau in tau_values:
    K, paths = propagador_MCMC(1.0, 0.5, 0.5, tau, 10000, 10000, 500, 0.1, pozo_inf)
    log_K.append(np.log(K))

slope, intercept, *_ = linregress(tau_values, log_K)
E0_est = -slope
E0_exact = (np.pi**2) / 2

# Gráfico del resultado
plt.figure(figsize=(8, 5))
plt.plot(tau_values, log_K, 'o-', label='log K_E (Metropolis)')
plt.plot(tau_values, intercept - E0_est * np.array(tau_values), '--',
         label=f'Ajuste: E0 ≈ {E0_est:.3f}')
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\log K_E$")
plt.title("Estimación de $E_0$ para pozo cuadrado infinito (1D)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Resultados
print(f"E0 estimado: {E0_est:.4f}")
print(f"E0 exacto:   {E0_exact:.4f}")
print(f"Error relativo: {abs(E0_est - E0_exact) / E0_exact * 100:.2f}%")
