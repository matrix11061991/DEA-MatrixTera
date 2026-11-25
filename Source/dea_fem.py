"""
DEA de physique - Résolution par éléments finis (1D, P1) des transferts
couplés chaleur-humidité (formulation de type Luikov / Palomo Del Barrio).

Partie A : support de toiture (chaleur seule), validé contre la série de Fourier.
Partie B : matériau poreux, système couplé non linéaire T-omega
           (éléments P1, Euler implicite, itérations de Picard).
Les lois de coefficients de la partie B sont ILLUSTRATIVES : à remplacer par
celles de l'article.
"""
import numpy as np

# ---------------------------------------------------------------- maillage
def mesh(L, n_el):
    z = np.linspace(0.0, L, n_el + 1)
    return z, np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)])

GP = np.array([-1.0, 1.0]) / np.sqrt(3.0)      # points de Gauss (2 pts)
GW = np.array([1.0, 1.0])

def shape(xi):
    return np.array([(1 - xi) / 2, (1 + xi) / 2]), np.array([-0.5, 0.5])

# ============================================================== PARTIE A
def fem_heat(L=1.0, n_el=100, a=1.0, T0=100.0, Tg=0.0, Td=0.0, dt=1e-3, t_end=0.1):
    """rho_C dT/dt = lambda d2T/dz2, Dirichlet aux deux bords. rho_C=1, lambda=a."""
    z, el = mesh(L, n_el); N = len(z)
    M = np.zeros((N, N)); K = np.zeros((N, N))
    for e in el:
        h = z[e[1]] - z[e[0]]
        Me = h / 6 * np.array([[2, 1], [1, 2]])
        Ke = a / h * np.array([[1, -1], [-1, 1]])
        M[np.ix_(e, e)] += Me; K[np.ix_(e, e)] += Ke
    A = M / dt + K
    T = np.full(N, T0); T[0], T[-1] = Tg, Td
    free = np.arange(1, N - 1)
    for _ in range(int(round(t_end / dt))):
        b = M @ T / dt
        # relèvement des valeurs imposées
        Tn = T.copy(); Tn[0], Tn[-1] = Tg, Td
        rhs = (b - A @ np.where(np.isin(np.arange(N), free), 0, Tn))[free]
        T[free] = np.linalg.solve(A[np.ix_(free, free)], rhs)
        T[0], T[-1] = Tg, Td
    return z, T

def exact(z, t, T0=100.0, L=1.0, a=1.0, n=400):
    s = np.zeros_like(z)
    for k in range(1, n, 2):
        s += 4*T0/(k*np.pi)*np.sin(k*np.pi*z/L)*np.exp(-a*(k*np.pi/L)**2*t)
    return s

# ============================================================== PARTIE B
class Coeffs:
    """Lois illustratives (unités adimensionnées)."""
    rho_C0, cw = 1.0, 2.0          # (rhoC)_g = rho_C0 + cw*omega
    lam0, lam1 = 1.0, 0.5          # lambda_g = lam0 + lam1*omega
    Lat, DVT   = 1.0, 0.05         # terme Lambda*D_VT
    Dvw        = 0.02              # couplage chaleur <- gradient d'humidité
    Dw0, beta  = 0.2, 2.0          # D_w = Dw0*exp(beta*omega)
    DT         = 0.01              # thermomigration
    K0         = 0.05              # K = K0*omega^3 (gravité/capillarité)
    @staticmethod
    def C(w):   return Coeffs.rho_C0 + Coeffs.cw * w
    @staticmethod
    def A(w):   return Coeffs.lam0 + Coeffs.lam1 * w + Coeffs.Lat * Coeffs.DVT
    @staticmethod
    def B(w):   return Coeffs.Dvw + 0 * w
    @staticmethod
    def Dw(w):  return Coeffs.Dw0 * np.exp(Coeffs.beta * w)
    @staticmethod
    def DT_(w): return Coeffs.DT + 0 * w
    @staticmethod
    def Kh(w):  return Coeffs.K0 * w**3

def fem_coupled(L=1.0, n_el=100, dt=1e-3, t_end=0.2,
                T_init=20.0, w_init=0.1, T_s=60.0, w_s=0.5,
                phi=lambda z, t: 0.0 * z, tol=1e-10, max_picard=30):
    """
    Inconnues u = (T, omega), P1, 2 ddl/noeud (ordre [T_0..T_N, w_0..w_N]).
    Bord z=0   : T = T_s, omega = w_s (Dirichlet)
    Bord z=L   : flux nul de chaleur et d'eau (Neumann homogène)
    Forme faible (v fonction test, IPP) :
      ∫ C T_t v + ∫ (A T_z + B w_z) v' = 0
      ∫ w_t v   + ∫ (Dw w_z + DT T_z) v' = ∫ K v' + ∫ phi v
    """
    z, el = mesh(L, n_el); N = len(z)
    T = np.full(N, T_init); w = np.full(N, w_init)
    T[0], w[0] = T_s, w_s
    free = np.concatenate([np.arange(1, N), N + np.arange(1, N)])
    t, hist = 0.0, []
    for _ in range(int(round(t_end / dt))):
        t += dt
        Tn, wn = T.copy(), w.copy()               # état au pas précédent
        Tk, wk = T.copy(), w.copy()               # itéré de Picard
        for it in range(max_picard):
            S = np.zeros((2 * N, 2 * N)); F = np.zeros(2 * N)
            for e in el:
                h = z[e[1]] - z[e[0]]
                ix = np.concatenate([e, N + e])
                Se = np.zeros((4, 4)); Fe = np.zeros(4)
                for xi, gw in zip(GP, GW):
                    Nn, dNxi = shape(xi); dN = dNxi * 2 / h; J = h / 2 * gw
                    wq = Nn @ wk[e]                # omega au point de Gauss
                    Tn_q, wn_q = Nn @ Tn[e], Nn @ wn[e]
                    zq = Nn @ z[e]
                    # --- équation de chaleur
                    Se[:2, :2] += J * (Coeffs.C(wq) / dt * np.outer(Nn, Nn)
                                       + Coeffs.A(wq) * np.outer(dN, dN))
                    Se[:2, 2:] += J * Coeffs.B(wq) * np.outer(dN, dN)
                    Fe[:2]     += J * Coeffs.C(wq) / dt * Tn_q * Nn
                    # --- équation d'eau
                    Se[2:, 2:] += J * (np.outer(Nn, Nn) / dt
                                       + Coeffs.Dw(wq) * np.outer(dN, dN))
                    Se[2:, :2] += J * Coeffs.DT_(wq) * np.outer(dN, dN)
                    Fe[2:]     += J * (wn_q / dt * Nn + Coeffs.Kh(wq) * dN
                                       + phi(zq, t) * Nn)
                S[np.ix_(ix, ix)] += Se; F[ix] += Fe
            u = np.concatenate([Tk, wk]); u[0], u[N] = T_s, w_s
            g = np.zeros(2 * N); g[0], g[N] = T_s, w_s
            rhs = (F - S @ g)[free]
            sol = np.linalg.solve(S[np.ix_(free, free)], rhs)
            unew = g.copy(); unew[free] = sol
            err = np.max(np.abs(unew - np.concatenate([Tk, wk])))
            Tk, wk = unew[:N], unew[N:]
            if err < tol: break
        T, w = Tk, wk
        hist.append(it + 1)
    return z, T, w, hist

if __name__ == "__main__":
    # ---- Partie A : validation
    for n_el in (25, 50, 100, 200):
        z, T = fem_heat(n_el=n_el, dt=1e-3)
        print(f"A  n_el={n_el:4d}  erreur max = {np.max(np.abs(T-exact(z,0.1))):.4f}")
    # ---- Partie B : convergence en maillage et nb d'itérations de Picard
    ref = fem_coupled(n_el=400, dt=2.5e-4)
    for n_el, dt in ((25, 1e-3), (50, 1e-3), (100, 1e-3)):
        z, T, w, h = fem_coupled(n_el=n_el, dt=dt)
        zr, Tr, wr, _ = ref
        eT = np.max(np.abs(T - np.interp(z, zr, Tr)))
        ew = np.max(np.abs(w - np.interp(z, zr, wr)))
        print(f"B  n_el={n_el:4d}  err T={eT:.4f}  err w={ew:.5f}  Picard moy={np.mean(h):.1f}")
    z, T, w, h = fem_coupled(n_el=100)
    print("B  T(z) :", np.round(T[::20], 2)); print("B  w(z) :", np.round(w[::20], 4))
