import numpy as np

if __name__ == "__main__":
    # load one-electron integrals
    T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")

    # load two electrons integrals
    J = np.loadtxt("J.mat").reshape(2 * T.shape)

    print(J)
