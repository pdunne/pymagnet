# Transforms

## Coordinate System Transforms

For certain magnetic sources, like a sphere or solenoid, their field equations are
more conveniently written in non-cartesian coordinates. The resulting vector fields
then need to be transformed back into cartesian coordinates of the form
$\mathbf{B} = B_x \mathbf{\hat{x}} + B_y \mathbf{\hat{y}} + B_z \mathbf{\hat{z}}$ or $\mathbf{B} = B_x \mathbf{\hat{x}} + B_y \mathbf{\hat{y}}$.

### Scalar Transforms

#### 2D: Cartesian ↔ Polar

**Cartesian to Polar:**

$$
\begin{align}
\rho &= \sqrt{x^2 + y^2} \\
\phi &= \arctan(y/x)
\end{align}
$$

**Polar to Cartesian:**

$$
\begin{align}
x &= \rho \cos \phi \\
y &= \rho \sin \phi
\end{align}
$$

#### 3D: Cartesian ↔ Spherical

**Cartesian to Spherical:**

$$
\begin{align}
r &= \sqrt{x^2 + y^2 + z^2} \\
\theta &= \arccos(z/r) \\
\phi &= \arctan(y/x)
\end{align}
$$

**Spherical to Cartesian:**

$$
\begin{align}
x &= r \sin\theta \cos\phi \\
y &= r \sin\theta \sin\phi \\
z &= r \cos\theta
\end{align}
$$

---

### Vector Transforms

#### Polar to Cartesian

$\mathbf{B} = B_r \mathbf{\hat{r}} + B_\phi \mathbf{\hat{\phi}}$

$$
\begin{align}
B_x &= B_r \cos \phi - B_\phi \sin \phi \\
B_y &= B_r \sin \phi + B_\phi \cos \phi
\end{align}
$$

#### Cylindrical to Cartesian

$\mathbf{B} = B_r \mathbf{\hat{r}} + B_\phi \mathbf{\hat{\phi}} + B_h \mathbf{\hat{h}}$

$$
\begin{align}
B_x &= B_r \cos \phi - B_\phi \sin \phi \\
B_y &= B_r \sin \phi + B_\phi \cos \phi \\
B_z &= B_h
\end{align}
$$

For a solenoid or cylinder $B_\phi = 0$, thus:

$$
\begin{align}
B_x &= B_r \cos \phi \\
B_y &= B_r \sin \phi \\
B_z &= B_h
\end{align}
$$

#### Spherical to Cartesian

$\mathbf{B} = B_r \mathbf{\hat{r}} + B_\theta \mathbf{\hat{\theta}} + B_\phi \mathbf{\hat{\phi}}$

$$
\begin{align}
B_x &= B_r \sin\theta \cos\phi + B_\theta \cos\theta\cos\phi - B_\phi \sin\phi \\
B_y &= B_r \sin\theta\sin\phi + B_\theta \cos\theta\sin\phi + B_\phi \cos\phi  \\
B_z &= B_r \cos\theta - B_\theta \sin\theta
\end{align}
$$

---

## API Reference

The conversion functions are available in `pymagnet.utils`:

```python
from pymagnet.utils import (
    cart2pol, pol2cart, vector_pol2cart,
    cart2sph, sph2cart, vector_sph2cart
)
```

### Scalar Conversion Functions

| Function | Description |
|----------|-------------|
| `cart2pol(x, y)` | Cartesian to polar: returns `(rho, phi)` |
| `pol2cart(rho, phi)` | Polar to cartesian: returns `(x, y)` |
| `cart2sph(x, y, z)` | Cartesian to spherical: returns `(r, theta, phi)` |
| `sph2cart(r, theta, phi)` | Spherical to cartesian: returns `(x, y, z)` |

### Vector Conversion Functions

| Function | Description |
|----------|-------------|
| `vector_pol2cart(Brho, Bphi, phi)` | Polar vector to cartesian: returns `(Bx, By)` |
| `vector_sph2cart(Br, Btheta, Bphi, theta, phi)` | Spherical vector to cartesian: returns `(Bx, By, Bz)` |

### Unit Conversion Functions

| Function | Description |
|----------|-------------|
| `get_unit_value_meter(unit)` | Get SI prefix factor for length units (e.g., "mm" → 1e-3) |
| `get_unit_value_tesla(unit)` | Get SI prefix factor for field units (e.g., "mT" → 1e-3) |

**Example:**

```python
from pymagnet.utils import get_unit_value_meter

factor = get_unit_value_meter("mm")
print(f"1 mm = {factor} m")  # Output: 1 mm = 0.001 m
```

---

## Polygon Geometry

### Centroid

The centroid of a non-self-intersecting closed polygon defined by $n$ vertices,
is the point $(C_x, C_y)$ where

$$
C_x = \frac{1}{6A}\sum_{i=0}^{n-1}(x_i+x_{i+1})(x_i\ y_{i+1} - x_{i+1}\ y_i)
$$

and

$$
C_y = \frac{1}{6A}\sum_{i=0}^{n-1}(y_i+y_{i+1})(x_i\ y_{i+1} - x_{i+1}\ y_i)
$$

### Signed Area

The signed area $A$ of a polygon:

$$
A = \frac{1}{2}\sum_{i=0}^{n-1} (x_i\ y_{i+1} - x_{i+1}\ y_i)
$$

Also written as:

$$
A = \frac{1}{2} \left(
    \begin{vmatrix}
        x_1 & x_2 \\
        y_1 & y_2  \\
    \end{vmatrix}
   + \begin{vmatrix}
        x_2 & x_3 \\
        y_2 & y_3  \\
    \end{vmatrix}
    + \cdots
    + \begin{vmatrix}
        x_n & x_1 \\
        y_n & y_1  \\
    \end{vmatrix}
   \right)
$$

For a convex polygon, if $A > 0$ the vertices are listed in counter-clockwise order, and clockwise if $A < 0$.
