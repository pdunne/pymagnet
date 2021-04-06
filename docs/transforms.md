# Transforms

## Coordinate System Transforms

### Scalar Transforms

### Vector Transforms

For certain magnetic sources, like a sphere or solenoid, their field equations are
more conveniently written in non-cartesian coordinates. The resulting vector fields
then need to be transformed back into cartesian coordinates of the form
$\mathbf{B} = B_x \mathbf{\hat{x}} + B_y \mathbf{\hat{y}} + B_z \mathbf{\hat{z}}$ or $\mathbf{B} = B_x \mathbf{\hat{x}} + B_y \mathbf{\hat{y}}$.

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


## Misc

Signed area of a polygon:

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
