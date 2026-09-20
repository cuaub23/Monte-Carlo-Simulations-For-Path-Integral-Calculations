def potencial_oscilador_armonico(x, m, omega=1):
    return 0.5*m*(omega**2)*(x**2)

tau_values = np.linspace(0.2, 2.0, 10)
log_K = []
for tau in tau_values:
    K, paths = propagador_MCMC(1.0, 0.0, 0.0, tau, 500, 200, 30, 0.3 , potencial_oscilador_armonico)
    log_K.append(np.log(K))

slope, intercept, *_ = linregress(tau_values, log_K)
E0_est = -slope
E0_exact = 0.5

# Graficando el ajuste
plt.figure(figsize=(8, 5))

plt.plot(tau_values, log_K, 'o-', label='log K_E (Metropolis)')
plt.plot(tau_values, intercept - E0_est * tau_values, '--', label=f'Ajuste: E0 ≈ {E0_est:.3f}')

plt.xlabel(r"$\tau$")
plt.ylabel(r"$\log K_E$")
plt.title("Estimación de $E_0$ con Monte Carlo por caminos (Metropolis)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Resultados
print(f"E0 estimado: {E0_est:.4f}")
print(f"E0 exacto:   {E0_exact:.4f}")
print(f"Error relativo: {abs(E0_est - E0_exact) / E0_exact * 100:.2f}%")
