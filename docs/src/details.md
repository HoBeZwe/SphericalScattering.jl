
# Further Details

---
## Units

SI units are employed everywhere.


---
## [Duality Relations](@id dualityRelations)

In many places the duality relations [[1, p. 120]](@ref refs)

| electric current | ``\rightarrow`` | magnetic current | 
|:---------------: | :-------------: | :--------------: |
| ``\bm e``        | ``\rightarrow`` | ``\bm h``        |
| ``\bm h``        | ``\rightarrow`` | ``-\bm e``       |
| ``\varepsilon``  | ``\rightarrow`` | ``\mu``          |
| ``\mu``          | ``\rightarrow`` | ``\varepsilon``  |

are employed in order to compute, e.g., the field of the magnetic current counterparts to the electric currents.


---
## [Rotations](@id rotationDetails)

Due to the rotational symmetry the fields for different orientations of the excitations can be computed via rotations.

#### Plane Wave

A plane wave excitation with arbitrary direction ``\hat{\bm k}`` and polarization ``\hat{\bm p}`` (forming a valid combination) 
can be related to the case ``\hat{\bm k} = \hat{\bm e}_z`` and polarization ``\hat{\bm p} = \hat{\bm e}_x`` by a rotation matrix
```math
\bm R = \begin{bmatrix} \hat{\bm p} & \hat{\bm k} \times \hat{\bm p} & \hat{\bm k}  \end{bmatrix} \,.
```
It constitutes a change from one right-handed orthogonal basis ``(\hat{\bm e}_x^\prime , \hat{\bm e}_y^\prime, \hat{\bm e}_z^\prime )`` to another right-handed orthogonal basis ``(\hat{\bm e}_x = \hat{\bm p}, \hat{\bm e}_y = \hat{\bm k} \times \hat{\bm p}, \hat{\bm e}_z = \hat{\bm k}``.
Hence, ``\bm R`` relates vectors (and points) ``\bm v'`` in the primed coordinate system ``(x',y',z')`` to vectors (and points) ``\bm v`` in the unprimed one ``(x,y,z)`` by 
```math
\bm v = \bm R  \bm v' \qquad \text{and} \qquad \bm v' = \bm R^\mathrm{T}  \bm v
```
leveraging that ``\bm R^{-1} = \bm R^\mathrm{T}``.
```@raw html
<br/>
```

#### Ring-Currents and Dipoles

Ring-currents and dipoles are characterized by only one vector ``\hat{\bm p}`` (defining the orientation).
An arbitrary ``\hat{\bm p}`` can be related to an initial one ``\hat{\bm p}_0`` (e.g., ``\hat{\bm p}_0 = \hat{\bm e}_z^\prime``) by a 
right-handed rotation around a rotation axis ``\hat{\bm a} = (a_x, a_y, a_z)`` by an angle ``\alpha``.
The rotation axis is given by
```math
\hat{\bm a} = \cfrac{\hat{\bm p}_0 \times \hat{\bm p}}{|\hat{\bm p}_0 \times \hat{\bm p}|}
```
and the angle ``\alpha`` is encoded in
```math
\cos(\alpha) = \hat{\bm p}_0 \cdot \hat{\bm p} \qquad \text{and} \qquad \sin(\alpha) = |\hat{\bm p}_0 \times \hat{\bm p}|\,.
```
The rotation matrix is then found by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) as
```math
\bm R = \mathbf{I} + \sin(\alpha)\bm{K} + (1 - \cos(\alpha)) \bm{K}^2
```
with
```math
\bm K = \begin{bmatrix} 0 & -a_z & a_y \\ a_z & 0 & -a_x \\ -a_y & a_x & 0 \end{bmatrix} \,.
```
Vectors (and points) ``\bm v'`` in the initial coordinate system ``(x',y',z')`` (with, e.g., ``\hat{\bm e}_z^\prime = \hat{\bm p}_0``) are then related to vectors (and points) ``\bm v`` in the unprimed one ``(x,y,z)`` by 
```math
\bm v = \bm R  \bm v' \qquad \text{and} \qquad \bm v' = \bm R^\mathrm{T}  \bm v
```
leveraging that ``\bm R^{-1} = \bm R^\mathrm{T}``.
