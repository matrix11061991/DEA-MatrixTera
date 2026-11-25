# DEA de physique [2013 - 2015]

Ce dépôt contient les codes de calcul et l'outil mis en place pendant la réalisation de mon DEA (Diplôme d'Études Approfondies, équivalent d'un master recherche).

Le travail porte sur la **résolution numérique d'un système d'équations aux dérivées partielles (EDP)** avec **conditions initiales** et **conditions aux limites**. Le sujet s'appuie sur l'article d'Elena Palomo Del Barrio relatif aux transferts couplés de chaleur et d'humidité dans les parois de bâtiment (support de toiture et matériau poreux).

## Sommaire

1. [Contexte physique](#1-contexte-physique)
2. [Le problème à résoudre](#2-le-problème-à-résoudre)
3. [Conditions initiales et conditions aux limites](#3-conditions-initiales-et-conditions-aux-limites)
4. [Méthodes de résolution](#4-méthodes-de-résolution)
5. [Exemple 1 : solution analytique par séparation des variables](#5-exemple-1--solution-analytique-par-séparation-des-variables)
6. [Exemple 2 : différences finies (schéma implicite)](#6-exemple-2--différences-finies-schéma-implicite)
7. [Validation et critères de qualité](#7-validation-et-critères-de-qualité)
8. [Notes sur les codes initiaux](#8-notes-sur-les-codes-initiaux)
9. [Références](#9-références)

---

## 1. Contexte physique

Une toiture est composée de plusieurs couches. On distingue ici :

- **Le support de toiture** (indice `s`), matériau solide où seule la chaleur est transportée par conduction.
- **Le matériau poreux** (indice `g`), où la chaleur **et** l'humidité circulent en même temps et s'influencent mutuellement (l'eau qui s'évapore absorbe de la chaleur, un gradient de température déplace la vapeur, etc.).

L'objectif est de prédire l'évolution dans le temps `t` et selon l'épaisseur `z` de :

- la température `T(z, t)` ;
- la teneur en eau `ω(z, t)` (pour le matériau poreux).

## 2. Le problème à résoudre

### 2.1 Support de toiture : équation de la chaleur

$$(\rho C)_s \frac{\partial T_s(z,t)}{\partial t} = \lambda_s \frac{\partial^2 T_s(z,t)}{\partial z^2}$$

| Symbole | Signification | Unité usuelle |
|---|---|---|
| `T_s(z,t)` | température dans le support | K ou °C |
| `(ρC)_s` | capacité thermique volumique (masse volumique × chaleur massique) | J·m⁻³·K⁻¹ |
| `λ_s` | conductivité thermique | W·m⁻¹·K⁻¹ |
| `z` | coordonnée dans l'épaisseur | m |
| `t` | temps | s |

C'est une équation de diffusion (parabolique). La grandeur `a = λ_s / (ρC)_s` est la **diffusivité thermique** (m²·s⁻¹) : plus elle est grande, plus la chaleur se propage vite.

> Remarque : dans la version image du README d'origine, le dénominateur est écrit `∂²x²` ; la variable d'espace est bien `z`, donc `∂z²`.

### 2.2 Matériau poreux : système couplé chaleur-humidité

$$
\begin{cases}
(\rho C)_g(\omega,t)\,\dfrac{\partial T}{\partial t} = \dfrac{\partial}{\partial z}\left[\big(\lambda_g(\omega,t) + \Lambda\, D_{VT}\big)\dfrac{\partial T}{\partial z} + D_{vw}(\omega,t)\dfrac{\partial \omega}{\partial z}\right] \\[3mm]
\dfrac{\partial \omega}{\partial t} = \dfrac{\partial}{\partial z}\left[D_w(\omega,t)\dfrac{\partial \omega}{\partial z} + D_T(\omega,t)\dfrac{\partial T}{\partial z}\right] - \dfrac{\partial K}{\partial z} + \varphi(z,t)
\end{cases}
$$

**Première équation : bilan d'énergie.**

| Terme | Rôle physique |
|---|---|
| `(ρC)_g(ω,t) ∂T/∂t` | stockage de chaleur ; dépend de la teneur en eau (l'eau a une forte capacité calorifique) |
| `λ_g(ω,t) ∂T/∂z` | conduction de la chaleur ; la conductivité augmente avec l'humidité |
| `Λ D_VT ∂T/∂z` | chaleur latente transportée par la vapeur sous l'effet du gradient de température (`Λ` : chaleur latente de vaporisation) |
| `D_vw ∂ω/∂z` | chaleur latente transportée par la vapeur sous l'effet du gradient d'humidité |

**Deuxième équation : bilan de masse d'eau.**

| Terme | Rôle physique |
|---|---|
| `∂ω/∂t` | variation de la teneur en eau |
| `D_w(ω,t) ∂ω/∂z` | diffusion de l'eau due au gradient d'humidité (loi de type Fick) |
| `D_T(ω,t) ∂T/∂z` | thermomigration : l'eau se déplace sous l'effet d'un gradient de température (effet Soret) |
| `∂K/∂z` | transport par gravité / capillarité (`K` : conductivité hydraulique) |
| `φ(z,t)` | terme source ou puits (condensation, évaporation, pluie, etc.) |

**Pourquoi c'est difficile.** Les coefficients `λ_g`, `(ρC)_g`, `D_w`, `D_T`, `D_vw` dépendent de la solution `ω` elle-même : le système est **non linéaire** et **fortement couplé**. Il n'a pas de solution analytique générale ; on le résout numériquement.

> À vérifier avec l'article : dans l'équation de chaleur, le terme en `D_vw` apparaît avec `∂ω/∂t` dans la version image d'origine. Dans la formulation classique de Luikov, ce terme porte le gradient d'espace `∂ω/∂z` (comme écrit ci-dessus, par cohérence dimensionnelle avec l'opérateur `∂/∂z` qui l'entoure). À confirmer avec la forme exacte de l'article.

## 3. Conditions initiales et conditions aux limites

Une EDP seule n'a pas de solution unique : il faut préciser l'état de départ et ce qui se passe aux bords du domaine `z ∈ [0, L]`.

**Conditions initiales (à t = 0)**

- `T(z, 0) = T₀` (profil uniforme ou connu)
- `ω(z, 0) = ω₀`

**Conditions aux limites** (trois types classiques)

| Type | Forme | Signification |
|---|---|---|
| Dirichlet | `T(0,t) = T_imposée` | température imposée en surface |
| Neumann | `-λ ∂T/∂z = q` | flux de chaleur imposé (flux nul = paroi adiabatique) |
| Robin (convectif) | `-λ ∂T/∂z = h (T - T_∞)` | échange avec l'air ambiant, `h` : coefficient d'échange |

Pour l'humidité, on a des conditions analogues : teneur en eau imposée, flux d'eau imposé (pluie, évaporation), ou échange de vapeur avec l'air.

## 4. Méthodes de résolution

| Méthode | Principe | Avantages | Limites |
|---|---|---|---|
| **Séparation des variables** | on cherche `T(z,t) = Z(z)·θ(t)` ; solution en série de Fourier | solution exacte, référence de validation | uniquement problèmes linéaires à coefficients constants |
| **Différences finies explicites** | `T^{n+1}_i = T^n_i + r (T^n_{i+1} - 2T^n_i + T^n_{i-1})` | très simple à coder | **stable seulement si** `r = a·Δt/Δz² ≤ 1/2` |
| **Différences finies implicites** | on résout un système linéaire tridiagonal à chaque pas | **inconditionnellement stable**, adapté aux problèmes raides | un système à résoudre par pas de temps |
| **Crank-Nicolson** | moyenne explicite/implicite | précision d'ordre 2 en temps | oscillations possibles sur condition initiale discontinue |

Pour le système couplé non linéaire, on retient en pratique le schéma **implicite**, avec **itérations de Picard** (on recalcule les coefficients avec la solution la plus récente jusqu'à convergence à chaque pas de temps).

## 5. Exemple 1 : solution analytique par séparation des variables

Cas de référence : équation de la chaleur sur `[0, L]`, température initiale uniforme `T₀`, bords maintenus à 0.

En posant `T(z,t) = Z(z)θ(t)`, on obtient `Z'' + k²Z = 0` et `θ' = -a k² θ`. Les conditions aux bords imposent `k = nπ/L`, d'où :

$$T(z,t) = \sum_{n\ \mathrm{impair}} \frac{4T_0}{n\pi}\,\sin\!\left(\frac{n\pi z}{L}\right)\exp\!\left(-a\left(\frac{n\pi}{L}\right)^2 t\right)$$

```python
import numpy as np

def solution_analytique(z, t, T0=100.0, L=1.0, a=1.0, n_termes=200):
    """Série de Fourier pour T(0,t)=T(L,t)=0 et T(z,0)=T0."""
    T = np.zeros_like(z)
    for n in range(1, n_termes, 2):          # seuls les n impairs contribuent
        k = n * np.pi / L
        T += 4 * T0 / (n * np.pi) * np.sin(k * z) * np.exp(-a * k**2 * t)
    return T
```

Cette solution sert de **référence** pour valider tout schéma numérique.

## 6. Exemple 2 : différences finies (schéma implicite)

```python
import numpy as np
from scipy.linalg import solve_banded

# --- Paramètres physiques ---
rho_C, lambda_s, L = 1.0, 1.0, 1.0
a = lambda_s / rho_C                  # diffusivité thermique

# --- Conditions initiale et aux limites ---
T0, T_gauche, T_droite = 100.0, 0.0, 0.0

# --- Discrétisation ---
N, dt, t_fin = 201, 1e-3, 0.1
z = np.linspace(0, L, N)
dz = z[1] - z[0]
r = a * dt / dz**2                    # nombre de Fourier de maille

# --- État initial ---
T = np.full(N, T0)
T[0], T[-1] = T_gauche, T_droite

# --- Matrice tridiagonale (identique à chaque pas, car coefficients constants) ---
ab = np.zeros((3, N))
ab[0, 2:] = -r                        # diagonale supérieure
ab[1, 1:-1] = 1 + 2 * r               # diagonale principale
ab[2, :-2] = -r                       # diagonale inférieure
ab[1, 0] = ab[1, -1] = 1              # lignes des conditions de Dirichlet
ab[0, 1] = 0
ab[2, -2] = 0

# --- Boucle en temps (nombre de pas fini) ---
t = 0.0
while t < t_fin - 1e-12:
    b = T.copy()
    b[0], b[-1] = T_gauche, T_droite
    T = solve_banded((1, 1), ab, b)
    t += dt
```

**Vérification effectuée :** avec ces paramètres, l'écart maximal avec la solution analytique à `t = 0,1` est d'environ 0,23 pour un profil de 100 (soit environ 0,2 %). À `z = L/2`, on obtient 47,68 en numérique contre 47,45 en analytique. L'écart provient surtout de la discontinuité de la condition initiale aux bords et de l'ordre 1 du schéma en temps.

## 7. Validation et critères de qualité

- **Comparaison avec la solution analytique** (section 5) sur le cas linéaire.
- **Étude de convergence en maillage** : diviser `Δz` et `Δt` et vérifier que l'erreur diminue (ordre 2 en espace, ordre 1 en temps pour l'implicite).
- **Conservation** : vérifier le bilan d'énergie et de masse d'eau sur le domaine.
- **Critère d'arrêt** : soit un temps final `t_fin`, soit une convergence vers le régime permanent, `max|T^{n+1} - T^n| < ε`.
- **Stabilité** : pour tout schéma explicite, contrôler `r ≤ 1/2`.

## 8. Notes sur les codes initiaux

Les deux exemples de code présents dans la première version de ce dépôt contenaient plusieurs défauts. Ils sont corrigés par les sections 5 et 6 :

1. **Le premier exemple n'utilise pas la séparation des variables.** C'est un schéma explicite en différences finies avec `np.gradient`.
2. **Mauvais rapport des constantes.** L'équation donne `∂T/∂t = λ/(ρC) · ∂²T/∂z²`, alors que le code utilisait `ρC/λ`.
3. **Instabilité.** Avec `Δz ≈ 10⁻³` et `Δt = 10⁻³`, on a `r ≈ 1000`, alors que le schéma explicite exige `r ≤ 0,5`. La solution diverge.
4. **Boucle `while True` sans fin.** La condition `T[-1] < T_inf` ne peut pas être vraie car la température ne passe pas sous la valeur imposée. Il faut une durée maximale ou un critère de convergence.
5. **Conditions aux limites absentes** du premier exemple.
6. **Deuxième exemple : double comptage.** Le terme `alpha*(...)` et le terme `dT_dt` représentent la même quantité et étaient additionnés ; de plus `alpha` était calculé avec le rapport inversé et `dT_dt` n'était pas multiplié par `dt`.
7. **Équation modélisée.** Seule l'équation de chaleur du support est traitée ; le système couplé `T-ω` du matériau poreux reste à implémenter (itérations de Picard, coefficients non linéaires).

## 9. Références

- E. Palomo Del Barrio, article source du modèle de transferts couplés chaleur-humidité (à compléter : titre, revue, année).
- A. V. Luikov, *Heat and Mass Transfer in Capillary-Porous Bodies*, Pergamon, 1966.
- J. R. Philip et D. A. de Vries, « Moisture movement in porous materials under temperature gradients », *Trans. AGU*, 1957.
