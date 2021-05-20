# Magnetic Field 1D Methods

!!! note
    - SI Units with the Sommerfeld convetion are used for this discussion:
    $\mathbf{B} = \mu_0 \left( \mathbf{H} + \mathbf{M}  \right)$[^1]
    - However, the Kennelly Convetion is used for the creation of each magnet object
    in the library:
    $\mathbf{B} = \mu_0\mathbf{H} + \mathbf{J}$[^1]
    - In free space $\mathbf{B} = \mu_0 \mathbf{H}$
    - In magnetised bodies, the demagnetising field $\mathbf{H_d} = - N \mathbf{M}$,
    where $N$ is the demagnetising factor.

## Cylinder

The magnetic field directly above the centre of a cylinder is[^2]:

<figure>
  <img src="../../img/cylinder.png" />
  <figcaption>Magnetic cylinder schematic</figcaption>
</figure>

$$
B_z = \frac{\mu_0 M_r}{2} \left[ \frac{z+L}{\sqrt{(z+L)^2 + R^2} } - \frac{z}{\sqrt{z^2 + R^2}} \right]
$$

## Cuboid

While for a cuboid, this equation is:

<figure>
  <img src="../../img/cuboid.png" />
  <figcaption>Magnetic cuboid</figcaption>
</figure>

$$
B_z = \frac{\mu_0 M_r}{2}  {\left[ \tan^{-1}{\left(
\frac{(z+L)\sqrt{a^2 + b^2 + (z+L)^2} }{ab}
\right)} - \tan^{-1}{\left( \frac{z\sqrt{a^2 + b^2 + z^2} }{ab}
\right)}
\right]}
$$

[^1]: J. M. D. Coey, Magnetism and Magnetic Materials (Cambridge University Press, 2010).
[^2]: E. P. Furlani, Permanent Magnet and Electromechanical Devices (Academic Press, San Diego, 2001).
