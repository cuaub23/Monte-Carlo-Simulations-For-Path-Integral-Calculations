from algoritmoMCMC import accion, metropolis, propagador_MCMC

def pozo_calabera(x, m, L=1.0, alpha=1.0, V0 = 10.0):
    return V0 * np.sin(np.pi * x / L)**2 + alpha*x**4

x = np.linspace(-3.3, 3.3 , 200)
plt.plot(x, pozo_calabera(x, 1))
plt.title(r"Potencial $V(x) = V_0 sin( \pi x /L)^2 + \alpha x^4$")
plt.grid()
plt.show()

tau_values = np.linspace(0.5, 2.0, 10)
log_K = []
for tau in tau_values:
    K, paths = propagador_MCMC(1.0, 0.5, 0.5, tau, 5000, 5000, 200, 0.1, pozo_calabera)
    log_K.append(np.log(K))

slope, intercept, *_ = linregress(tau_values, log_K)
E0_est = -slope
E0_df = 12.359163564778655

# Gráfica
plt.figure(figsize=(8, 5))
plt.plot(tau_values, log_K, 'o-', label='log K_E')
plt.plot(tau_values, intercept - E0_est * np.array(tau_values), '--',
         label=f'Ajuste: E0 ≈ {E0_est:.3f}')
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\log K_E$")
plt.title(r"Estimación de $E_0$ para $V(x) = V_0 \sin^2(\pi x/L) + \alpha x^4$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

print(f"E0 estimado por MCMC: {E0_est:.4f}")
print(f"E0 estimado por metodo de diferencias finitas : {E0_df:.4f}")
print(f"Error relativo: {abs(E0_est - E0_df) / E0_df * 100:.2f}%")
