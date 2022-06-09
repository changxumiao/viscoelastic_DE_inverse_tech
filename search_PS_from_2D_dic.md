#  Introduction
Inspired by the [EUCLID project][EUCLID], a image (2D)-based inverse technique will be developed here. The inverse technique will enable a search for **Prony Series** from loading tests on arbitrary geometry, which contains info of 2D deformation and boundary forces.

The technique addresses issues such as:
- standard test cannot be carried out
- working machine cannot be ceased to test carried material
- temperature-sensitive working environment, while prediction is simultaneously needed.
- to add

# Assumption
To limit our scope, the material is assumed to a double network: permanent network and dynamic network. Both networks are Neo-Hookean material. Their shares of the cross-link density are known.

# Workflow
1. Time-related deformation and loading info on arbitrary geometry
    - Displacement data on nodes
    - Boundary forces on nodes
2. Strain field
3. Predict stress field with overdetermined selection of Prony Series
4. Build cost function
    - Nodal force residual on free degree of freedom
    - Reaction force residual on loading degree of freedom
5. Select Prony Series with accuracy and sparsity.

# Python coding
Ideas:
Each process should be a mini-service in form of Class. Each test is a class instance.
Given:
- Coordinate of nodes
- Displacement data on nodes
- Connectivity of nodes: from nodes to elements
- Boundary forces on nodes
- Material properties: Neo-Hookian
1. A class of extracting geometry
    - input: data frame
    - output: nodes ID, coordinates, displacement
    - Sort coordinates of nodes
2. A class of elements
    - input: data frame of connectivity, nodes ID, coordinate, displacement
    - output: elements, connectivity, deformation gradient, velocity gradient
    - form elements from nodes and connectivity
    - deformation gradient
    - velocity gradient
3. A class of constitutive law
    - input: deformation gradient
    - output: PK2 stress tensor
    - alternative strain field
        - inverse of deformation gradient
        - right Cauchy-Green deformation tensor
        - Green-Lagrange strain tensor
        - invariant of Cauchy-Green deformation tensor
    - Neo-Hookian
    - Linear elasticity

    - stress alternatives (from PK2 to PK1 and Cauchy)
    - alternative stress field
        - from PK2 to PK1
        - from PK2 to Cauchy
4. A class of dynamics law
    - input: velocity gradient, Prony Series
    - Viscousity
    - Dynamic network (Prony Series)

5. A module listing element family
    - element dimension
    - number of nodes per element



[EUCLID]: https://github.com/EUCLID-code