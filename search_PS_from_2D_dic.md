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



[EUCLID]: https://github.com/EUCLID-code