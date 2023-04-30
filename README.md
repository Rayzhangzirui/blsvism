# Binary Level Set Method for Variational Implicit Solvation Model

The repository contain codes for the binary level set method for the Variational Implicit Solvaiton Model (VISM).

### Shape optimization
The goal is to find an interface $\Gamma = \partial \Omega \subset R^3$  that minimizes some energy functional $G[\Gamma]$, which consists of the surface area and volumne integral of $\Omega^c$ region.

$$
G[\Gamma] = \gamma \text{Area}(\Gamma) + \int_{\Omega^c} U(x) dx
$$

In VISM, $\gamma$ is the surface tension and $U(x)$ includes the Lennard-Jones potential and electrostatic potential.

### Binary Level Set Method
The interface is represented by a binary level set function: -1 inside and +1 outside. The surface area is approximated by the "heat content": convolution of indicator function with a compact radially symmetric kernel $K_\epsilon$ of radius $\epsilon$:

$$
\text{Area}(\Gamma) = C \int \int 1_{\Omega}(x) K_{\epsilon}(|x-y|)  1_{\Omega^c}(y)dxdy + O(\epsilon^2)
$$

After discretization of $G$ using midpoint rule in a grid of size $h$, we arrived at a discrete energy with $O(h)$ error

### Optimization scheme
The energy is minimized by iteratively "flipping the pixel" in the steepest descent direction of the energy, using an efficient, customized heap data structure.

### Usage

To run numerical examples in the paper:

kerneltest: convergence of surface area estimation.

atomtest: convergence of 1 atom.

atom2test: 2 atom example.

vismtest: biomolecule example.

vismfft: comparison with threshold-dynamics with fft.



### Reference
Zhang, Z., Cheng, L.-T., 2021. Binary Level Set Method for Variational Implicit Solvation Model. https://arxiv.org/abs/2110.12815