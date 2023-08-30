import matplotlib.animation as anm
import matplotlib.pyplot as plt
import numpy as np

xhalfrange, xpoints, dt, steps, thresh = 37, 1024, 0.1, 10000, 1e-12
dx, nstates = 2 * xhalfrange / (xpoints - 1), 12
legend, plotlim = True, 1e-8

x = np.linspace(-xhalfrange, xhalfrange, xpoints)
k = 2 * np.fft.fftfreq(len(x), dx) * np.pi

psi0, V = np.exp(-(x - 1)**2), 0.5 * x**2

def energy(wfn):
    Ek = 0.5 * np.conj(wfn) * np.fft.ifft(k**2 * np.fft.fft(wfn))
    Er = np.conj(wfn) * V * wfn; return np.sum(Ek + Er).real * dx

if __name__ == "__main__":
    R, K, states = np.exp(-0.5 * V * dt), np.exp(-0.5 * k**2 * dt), []

    for i in range(nstates):
        psi = [0j + psi0 / np.sqrt(np.sum(np.abs(psi0)**2) * dx)]

        for _ in range(steps):
            psi.append(R * np.fft.ifft(K * np.fft.fft(R * psi[-1])));

            for j in range(i):
                psi[-1] -= np.sum(np.conj(states[j][-1]) * psi[-1]) * dx * states[j][-1]
                # psi[0] -= np.sum(np.conj(states[0][-1] * np.exp(0.5j * np.pi)) * psi[0]) * dx * states[0][-1] * np.exp(0.5j * np.pi)


            psi[-1] /= np.sqrt(np.sum(np.abs(psi[-1])**2) * dx)

            # if np.sum(np.abs(psi[-1] - psi[-2])**2) < thresh or np.abs(energy(psi[-1]) - energy(psi[-2])) < thresh: break

        states.append(psi); print("E_{}:".format(i), energy(psi[-1]))

    D = [[energy(psi) + np.abs(psi)**2 for psi in state] for state in states]

    xmin = np.min([np.min([np.min(x[np.abs(psi)**2 > plotlim]) for psi in Si]) for Si in states])
    xmax = np.max([np.max([np.max(x[np.abs(psi)**2 > plotlim]) for psi in Si]) for Si in states])
    [fig, ax], xrange = plt.subplots(), np.max([np.abs(xmin), np.abs(xmax)]); plt.tight_layout()

    ax.set_xlim(-xrange, xrange);
    ax.set_ylim(0, 1.3 * np.max(D[-1][-1]))
    ax.plot(x, V); psiplots = [ax.plot(x, Di[0])[0] for Di in D]

    def update(j):
        labels = ["$\\psi_{{{},{}}}$".format(i, j if j < len(D[i]) else len(D[i]) - 1) for i in range(len(psiplots))]
        for i in range(len(psiplots)): psiplots[i].set_ydata(D[i][j if j < len(D[i]) else -1])
        if legend: ax.legend(psiplots, labels, loc="upper right")

    ani = anm.FuncAnimation(fig, update, frames=np.max([len(state) for state in states]), interval=10); plt.show(); plt.close("all")
