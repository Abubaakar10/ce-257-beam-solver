# CE 257 Capstone Project: Technical Implementation Report

## 1. Introduction
This project provides a structural analysis application that computes the internal forces of an arbitrary continuous beam using the Three-Moment Equation. The software allows users to draw the beam system, apply loads, and visualize the Bending Moment Diagram (BMD) and Shear Force Diagram (SFD) through an interactive Graphical User Interface (GUI).

## 2. Mathematical Derivation: The Three-Moment Equation
The Three-Moment Equation expresses the relationship between the internal bending moments at three consecutive supports of a continuous beam. It is derived from the principle of consistent deformations, which states that the slope of the deflection curve must be continuous across any intermediate support.

For any two adjacent spans, say $L_1$ (between supports 1 and 2) and $L_2$ (between supports 2 and 3), let $M_1$, $M_2$, and $M_3$ be the internal moments at the supports. The slope $\theta$ at support 2 from the left span must equal the slope from the right span:
$$\theta_{2, left} = \theta_{2, right}$$

Using the Moment-Area Theorem or Conjugate Beam Method, the slope at the support due to applied loads on simply supported spans is balanced by the slope caused by the unknown support moments. This balance yields the general form:

$$M_{i-1}L_{i} + 2M_{i}(L_{i} + L_{i+1}) + M_{i+1}L_{i+1} = -6 \left( \frac{A_i \bar{x}_i}{L_i} + \frac{A_{i+1} \bar{x}_{i+1}}{L_{i+1}} \right)$$

Where:
- $M_{i-1}, M_{i}, M_{i+1}$: Bending moments at supports $i-1$, $i$, and $i+1$.
- $L_i, L_{i+1}$: Lengths of spans $i$ and $i+1$.
- $A_i, A_{i+1}$: Area of the simply-supported bending moment diagrams for spans $i$ and $i+1$ due to applied external loads.
- $\bar{x}_i$: Distance from the centroid of $A_i$ to support $i-1$.
- $\bar{x}_{i+1}$: Distance from the centroid of $A_{i+1}$ to support $i+2$.

**Right-Hand Side Expressions for Common Loads:**
- **Uniform Load ($w$)**: The RHS term for a single span evaluates to $\frac{w L^3}{4}$.
- **Point Load ($P$)**: At a distance $a$ from the left support (and $b$ from the right support, $a+b=L$), the RHS term for the span to the left of the center support is $\frac{P a}{L}(L^2 - a^2)$, and for the span to the right, it is $\frac{P b}{L}(L^2 - b^2)$.

## 3. Software Architecture and Logic

The software is divided into three major modules:
1. `math_engine.py`: Defines the `BeamSolver` class that constructs a system of linear equations ($Ax = B$) based on the Three-Moment Equation. It solves for the moments using `numpy.linalg.solve` and generates arrays for the internal forces across the beam's length.
2. `visualization.py`: Utilizes `matplotlib.pyplot` to render the computed Shear Force Diagram (SFD) and Bending Moment Diagram (BMD) using standard civil engineering plotting conventions.
3. `gui.py` & `main.py`: Provides the Tkinter interface to gather user inputs (spans and loads) and bridges the backend computations with the visualization.

### Flowchart of Software Logic

```mermaid
flowchart TD
    A[User launches main.py] --> B[Tkinter GUI initializes]
    B --> C{User Action}
    C -->|Add Span| D[Update internal span list]
    C -->|Add Load| E[Update internal load dictionary]
    D --> F[Redraw GUI Canvas Grid]
    E --> F
    C -->|Click 'Solve & Plot'| G[Initialize BeamSolver]
    G --> H[Generate Matrix A & Vector B from loads/spans]
    H --> I[Solve Ax = B via numpy]
    I --> J[Compute support moments]
    J --> K[Calculate V(x) and M(x) across all points]
    K --> L[Pass arrays to matplotlib]
    L --> M[Render SFD and BMD]
    M --> C
```

## 4. Conclusion
The implementation successfully bridges the theoretical concepts of structural analysis and programming. The Three-Moment Equation serves as a robust solver for statically indeterminate continuous beams, and the matrix manipulation handled by `numpy` provides near-instantaneous computation times.
