# 2D Magnets

!!! note
    - SI Units with the Sommerfeld convetion are used for this discussion:
    $\mathbf{B} = \mu_0 \left( \mathbf{H} + \mathbf{M}  \right)$[^1]
    - However, the Kennelly Convetion is used for the creation of each magnet object
    in the library:
    $\mathbf{B} = \mu_0\mathbf{H} + \mathbf{J}$[^1]
    - In free space $\mathbf{B} = \mu_0 \mathbf{H}$
    - In magnetised bodies, the demagnetising field $\mathbf{H_d} = - N \mathbf{M}$,
    where $N$ is the demagnetising factor.

For infinitely long objects, the field problems can be approximated in 2D, where the magnet field consists of:

$$
\mathbf{B} = B_x \mathbf{\hat{x}} + B_y\mathbf{\hat{y}}
$$

## Rectangles

<figure>
    <img src="../../img/2d_rectangle.png" width="300" />
    <figcaption>2D Magnet Rectangle</figcaption>
</figure>

The magnetic field due to rectangle magnetised in $y$ is[^2]:

$$
B_x = \frac{\mu_0 M_r}{4\pi} \left[\ln {\left(
\frac{{\left(x+a\right)}^2 + {\left(y-b\right)}^2}{{\left(x+a\right)}^2
+{\left(y+b\right)}^2}
\right)}
-\ln{\left(
\frac{{\left(x-a\right)}^2+{\left(y-b\right)}^2}{ {\left(x-a\right)}^2 +
{\left(y+b\right)}^2}
\right)}\right]
$$

$$
B_y = \frac{\mu_0 M_r}{2\pi}
\left[{\tan}^{-1}{\left( \frac{2b \left(x+a\right)}{y^2-b^2+{\left(x+a\right)}^2}
\right)}
- {\tan}^{-1}{\left(\frac{2b\left(x-a\right)}{y^2-b^2+{\left(x-a\right)}^2}\right)}\right]
$$

or magnetised in $x$ is:

TODO:

$$
B_x = \frac{\mu_0 M_r}{2\pi} []
$$

$$
B_y = \frac{\mu_0 M_r}{4\pi}[]
$$

## Biaxial Rods (Circle)

<figure>
    <img src="../../img/2d_circle.png" width="300" />
</figure>

A long bipolar rod of radius $a$ can be approximated as circular source. The magnetic
stray field is most conveniently written in polar coordinates, as[^2]:

$$
\mathbf{B} = \frac{\mu_0 M_r}{2} \left( \frac{a^2}{r^2}\right) \left[
    \cos(\phi) \mathbf{\hat{r}} + \sin(\phi) \mathbf{\hat{\phi}}
     \right]
$$

[^1]: J. M. D. Coey, Magnetism and Magnetic Materials (Cambridge University Press, 2010).
[^2]: E. P. Furlani, Permanent Magnet and Electromechanical Devices (Academic Press, San Diego, 2001).
