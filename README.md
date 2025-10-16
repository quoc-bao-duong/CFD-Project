# CFD-Project
## Advection-Diffusion Equation 2D

$$
\rho \frac{\partial \phi}{\partial t} + \rho u \frac{\partial \phi}{\partial x} + \rho v \frac{\partial \phi}{\partial y}
= \Gamma \frac{\partial^2 \phi}{\partial x^2} + \Gamma \frac{\partial^2 \phi}{\partial y^2}
$$

**Velocity Field:**  
u = x, v = -y 

**Boundary Conditions:**

$$
\left( \frac{\partial \phi}{\partial y} \right)_{y=0} = 0, \quad
\phi_{y=1} = 0, \quad
\phi_{x=0} = 1 - y, \quad
\left( \frac{\partial \phi}{\partial y} \right)_{x=1} = 0
$$

---

# ⚙️ Approach

Finite Volume Method using:
- **First-order Upwind**
- **Central Scheme**
- **Forward Euler Explicit**

Implementation: **MATLAB**, Visualization: **TECPLOT**

---

# 🎯 Goal

Solve the 2D Advection-Diffusion Equation using **Finite Volume Method (FVM)**  
and **Visualize the Results**.
