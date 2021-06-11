# Methods

## Magnetic Fields

The overall approach is to

1. Instantiate a set of magnets
2. Generate an array of points to be calculated
3. Loop over each magnet, calcuate the field at each point and sum this to the
total field.
4. Draw the resulting data as a line, contour, slice, or volume plot.

```mermaid
flowchart
    1.Create_Magnets-->2.Create_Points
    subgraph 3.Calculate_B
        direction TB
        Transform_Local_Coords -->
        Calc_B -->
        Transform_Global_Coords -->
        Sum_Fields
    end
    subgraph Globals
        direction TB
        Registry --> Magnet1,Magnet2...
        Points --> x,y,z
        Field --> x,y,z,n
    end
    subgraph 4.Plot
        direction TB
        Draw_Field --> Draw_Magnets
    end
    1.Create_Magnets -->Registry
    2.Create_Points -->Points
    Registry -->|Loop over magnets| 3.Calculate_B
    Sum_Fields -->Field
    Points --> Transform_Local_Coords
    3.Calculate_B --> 4.Plot
    Globals --> 4.Plot
```

## Forces and Torques

1. Find the faces of a magnet with $|\mathbf{M}\cdot \mathbf{\hat{n}}| > 0$

2. Generate a grid of points on each of those faces

3. Calculate the magnetic field at these points due to all other magnets

4. Calculate the mean force and torque on each face

5. Sum the forces and torques on each face

## Classes

At the top of the hierarchy is the Registry class which records a set of `Weakref`
references to instances of each class, which is used for the `Magnet` and `Polyhedron`
child classes.

```mermaid
classDiagram
Registry <|-- Magnet
Registry <|-- Polyhedron
class Registry{
    +set instances
    -list _class_instances
    +print_instances()
    +get_instances() set
    +get_num_instances() int
    +reset()
}
```

### Magnet Classes

```mermaid
classDiagram
    Registry <|-- Magnet
    Magnet <|-- Magnet_2D
    Magnet_2D <|-- Rectangle
    Rectangle <|-- Square
    Magnet_2D <|-- Circle
    Magnet <|-- Magnet_3D
    Magnet_3D <|-- Prism
    Prism <|-- Cube
    Magnet_3D <|-- Cylinder
    Magnet_3D <|-- Sphere
    Magnet_3D <|-- Mesh
    class Registry{
        +set instances
        +get_instances() 
        +reset()
    }
    class Magnet{
        +rest_magnets()
        +list_magnets()
    }
    class Magnet_2D{
        +get_field()
    }

    class Rectangle{
    }

    class Circle{
    }
    class Magnet_3D{
        +get_field()
        +get_force_torque()
    }
    class Prism{
    }
    class Cylinder{
    }
    class Sphere{
    }
    class Mesh{
    }
```

### Mesh Class

The Mesh magnet class performs

### Quaternion Class

This is a convenience class for performing rotations of points/vectors about arbitrary axes.

## Plot Methods

### Draw Magnets on Plot
