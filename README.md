# DEA in physics [2013 - 2015]
This deposit contains the calculation codes and the tool that I set up during the realization of my DEA(master's degree).
This focuses on solving a system of partial differential equations with **initial** conditions and **boundary** conditions.
## The problem to be solve
<img
  src="https://latex.codecogs.com/svg.image?(\rho&space;C)_{s}\frac{\partial&space;T_{s}(z,t)}{\partial&space;t}&space;=&space;\lambda&space;_{s}\frac{\partial&space;^{2}T_{s}(z,t)}{\partial&space;^{2}x^{2}}"
/>

<img
  src="https://latex.codecogs.com/svg.image?\left\{\begin{matrix}(\rho&space;C)_{g}(\omega&space;,t)\frac{\partial&space;T(z,t)}{\partial&space;t}&space;=&space;\frac{\partial}{\partial&space;z}\left&space;[&space;(\lambda&space;_{g}(\omega&space;,t)&plus;\wedge&space;D_{VT})\frac{\partial&space;T(z,t)}{\partial&space;z}&space;&plus;&space;D_{vw}(\omega&space;,t)\frac{\partial&space;\omega&space;(z,t)}{\partial&space;t}&space;\right&space;]&space;&&space;&space;\\\frac{\partial&space;\omega&space;(z,t)}{\partial&space;t}&space;=&space;\frac{\partial&space;}{\partial&space;z}\left&space;[&space;D_{w}(\omega&space;,t)\frac{\partial&space;\omega(z,t)&space;}{\partial&space;z}&space;&plus;&space;D_{T}(\omega&space;,t)\frac{\partial&space;T(z,t)}{\partial&space;z}&space;\right&space;]-\frac{\partial&space;K(z,t)}{\partial&space;z}&space;&plus;&space;\varphi&space;(z,t)\end{matrix}\right."
/>

## Example of solving the equation governing the roof support
This is an example of Python code that can be used to solve the roof support equation using the separation of variables method:
```python
import numpy as np

# Constantes de l'équation
rho_C = 1.0
lambda_s = 1.0

# Conditions initiales
T_0 = 100.0

# Conditions aux limites
T_inf = 0.0

# Discrétisation du domaine
N = 1000
z = np.linspace(0, 1, N)

# Pas de temps
dt = 0.001

# Solution initiale
T = T_0*np.ones(N)

# Boucle principale de résolution
while True:
    # Calcul de la dérivée de T en z
    dT_dz = np.gradient(T, z)
    
    # Calcul de la dérivée seconde de T en z
    d2T_dz2 = np.gradient(dT_dz, z)
    
    # Calcul de la dérivée de T en t
    dT_dt = (rho_C/lambda_s)*d2T_dz2
    
    # Mise à jour de la solution
    T = T + dT_dt*dt
    
    # Condition de sortie de la boucle
    if T[-1] < T_inf:
        break
```
This code solves the differential equation using an infinite loop, which is broken when the temperature reaches the value of the boundary condition T_inf. The main loop of the solution consists of calculating the derivatives of T with respect to z and t, and then updating the solution using these derivatives.

This is an example of Python code that uses the finite differences method to solve the differential equation governing the roof support:
```python
import numpy as np

# Constantes de l'équation
rho_C = 1.0
lambda_s = 1.0

# Conditions initiales
T_0 = 100.0

# Conditions aux limites
T_inf = 0.0

# Discrétisation du domaine
N = 1000
z = np.linspace(0, 1, N)
dz = z[1] - z[0]

# Pas de temps
dt = 0.001

# Solution initiale
T = T_0*np.ones(N)

# Coefficients de la méthode des différences finies
alpha = rho_C*dt/(lambda_s*dz**2)

# Boucle principale de résolution
while True:
    # Calcul de la dérivée seconde de T en z
    d2T_dz2 = (T[2:] - 2*T[1:-1] + T[:-2])/(dz**2)
    
    # Calcul de la dérivée de T en t
    dT_dt = (rho_C/lambda_s)*d2T_dz2
    
    # Mise à jour de la solution
    T[1:-1] = T[1:-1] + alpha*(T[2:] - 2*T[1:-1] + T[:-2]) + dT_dt
    
    # Conditions aux limites
    T[0] = T_0
    T[-1] = T_inf
    
    # Condition de sortie de la boucle
    if T[-1] < T_inf:
        break
```
This code solves the differential equation using the finite differences method to approximate the derivatives of T with respect to z and t. The boundary conditions are applied by updating the values of T at the edges of the domain (z = 0 and z = 1). The main loop of the solution consists of updating the solution using the values of T at the interior points of the domain, and then applying the boundary conditions at the edge points.
