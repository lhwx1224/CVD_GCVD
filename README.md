# Characteristic Value Decomposition (CVD)

## Graphical Abstract
![](https://github.com/lhwx1224/CVD_GCVD/blob/main/Figures/Docs/Fig_CVD_GraphicalAbstract.jpg)

## Formulation
### General form of the CVD
Given two data matrices $\mathbf{Y}_1\in\mathbb{R}^{m\times n}$ and $\mathbf{Y}_2\in\mathbb{R}^{m\times n}$, the CVD finds the common modes and corresponding coordinates through simultaneously decomposing
$\mathbf{Y}_1 = \mathbf{U}\mathbf{\Sigma_1}\mathbf{V}^H$  and  $\mathbf{Y}_2 = \mathbf{U}\mathbf{\Sigma_2}\mathbf{V}^H$. 

### Practical form of the CVD
Given a generic nonlinear dynamical system 

$$\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u}, \mathbf{w}, t;\mathbf{p})$$

and a nonlinear observation model

$$\mathbf{y} = \mathbf{h}(\mathbf{x}, \mathbf{u}, \mathbf{n})$$

Based on the measurement of the system $\mathbf{Y}$ and the corresponding rate estimate $\dot{\mathbf{Y}}$, the best-approximate system mode estimates can be conducted through the CVD 
$\mathbf{Y} = {\mathbf{U}}\mathbf{\Sigma}{\mathbf{V}}^H$  and  $\dot{\mathbf{Y}} = {\mathbf{U}}\dot{\mathbf{\Sigma}}{\mathbf{V}}^H$. 

The estimated mode shape of the system is given by $\hat{\mathbf{\Phi}} = \mathbf{V}$, the estimated modal coordinates are given by $\hat{\mathbf{\Xi}} = {\mathbf{U}}$, and the poles of the system is estimated from $\hat{\mathbf{\Lambda}}=\mathbf{Sigma}\dot{\mathbf{Sigma}}^{-1}$. 
In conventional modal analysis, linear lumped-parameter dynamical systems either from assumptions or from a discretized technique are usually considered with the following form

$$\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$$

and a linear observer, e.g., sensors that work in their linear range and subjected to additive measurement noise

$${\mathbf{y}} = \mathbf{H}\mathbf{x} + \mathbf{n}$$

CVD can then be used based on meausrements, $\mathbf{Y} = [\mathbf{y}_1\\, \mathbf{y}_2\\, \dots\\, \mathbf{y}_m]$, a stack of the measurements into a matrix and its derivative $\dot{\mathbf{Y}}$, estimated from gradient approximations.
When the states $\mathbf{Y} = \mathbf{X}$ and state rates $\dot{\mathbf{Y}} = \dot{\mathbf{X}}$ are in use, the identified modal parameters will be the true modal parameters of the underlying linear system.

The solution of the CVD is obtained from the eigendecomposition of the quotient matrix $\mathbf{R} = \mathbf{K}_2\mathbf{K}_1^{-1}$,

$$\mathbf{R}\mathbf{V} = \mathbf{V}\mathbf{S}$$

where $\mathbf{K}_2 \triangleq \dot{\mathbf{Y}}^H\mathbf{Y}$ and $\mathbf{K}_1 \triangleq \mathbf{Y}^H\mathbf{Y}$; the characteristic modes are columns of $\mathbf{V}$, the corresponding modal coordinates (unnormalized) are in $\mathbf{P}_1 = \mathbf{Y}\mathbf{V^{-H}}$. The characteristic coordinates are defined as the normalized modal coordinates $\mathbf{U} = \mathbf{P}_1\mathbf{S}_1^{-1}$, with $\mathbf{S}_1 = \mathrm{diag}([\\|\mathbf{p}_1\\|\\,\\|\mathbf{p}_2\\|\\,\dots\\,\\|\mathbf{p}_n\\|])$ and $\\|\cdot\\|$ being the $l_2$ norm. Finally, the characteristic values are obtained via $\mathbf{S}_2\mathbf{S}_1^{-1}$.

### Generalization (GCVD)
Generalized CVD generalizes the quotient matrix to a balanced one, which can be taylored to fit CVD for modal analysis under various forcing conditions.

## Examples
Find the provided examples in the Example folder to use the CVD/GCVD.

## Cite the work
```
Add bibliography information here
```
