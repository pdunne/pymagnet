# Forces and Torques

The total force acting on a magnetic body with magnetisation $\mathbf{M}$ due to an external field $\mathbf{B'}$ can be divided into the sum of a volume and surface integral:

$$
\mathbf{F} = -\int_V (\nabla \cdot \mathbf{M}) \mathbf{B'} \: dv + \int_S (\mathbf{M} \cdot \mathbf{\hat{n}}) \mathbf{B'} \: ds
$$

which can be rewritten as

$$
\mathbf{F} = \int_V \rho_m \mathbf{B'} \: dv + \int_S \sigma_m \mathbf{B'} \: ds
$$

where $\rho_m = -\nabla \cdot \mathbf{M}$ and $\sigma_m =  \mathbf{M} \cdot \mathbf{\hat{n}}$ are the volume and surface charge densities.

For a uniformly magnetised body, $\rho_m = -\nabla \cdot \mathbf{M} = 0$, thus the total magnetic force acting on a magnet due to an external field $\mathbf{B'}$ can be deduced using the surface charge alone:

$$
\mathbf{F} =  \int_S \sigma_m \mathbf{B'} \: ds
$$

Similarly, the torque on a uniformly magnetised body can be determined using

$$
\mathbf{T} = \int_S \sigma_m (\mathbf{r} \times \mathbf{B'}) \: ds
$$

where $\mathbf{r}$ is the vector from the point about which the torque is computed, usually the centre of mass.
