# Characteristic Value Decomposition (CVD)

## Graphical Abstract
![](https://github.com/lhwx1224/CVD_GCVD/blob/main/Figures/Docs/Fig_CVD_GraphicalAbstract.jpg)

## Formulation
### General form of the CVD
Given two data matrices $\mathbf{Y}_1\in\mathbb{R}^{m\times n}$ and $\mathbf{Y}_2\in\mathbb{R}^{m\times n}$, the CVD finds the common modes and corresponding coordinates through simultaneously decomposing
$\mathbf{Y}_1 = \mathbf{U}\mathbf{\Sigma_1}\mathbf{V}^H$  and  $\mathbf{Y}_2 = \mathbf{U}\mathbf{\Sigma_2}\mathbf{V}^H$. 

### Practical form of the CVD
Given the measurement of the system $\mathbf{Y}$ and the corresponding rate estimate $\dot{\mathbf{Y}}$, the system mode estimates $\hat{\mathbf{\Phi}}$, the modal coordinates $\hat{\mathbf{U}}$, and their corresponding modal amplitudes $\mathbf{Sigma}$ and $\dot{\mathbf{Sigma}}$ are obtained through 
$\mathbf{Y} = \hat{\mathbf{U}}\mathbf{\Sigma}\hat{\mathbf{\Phi}}^H$  and  $\dot{\mathbf{Y}} = \hat{\mathbf{U}}\dot{\mathbf{\Sigma}}\hat{\mathbf{\Phi}}^H$. 

When the states and state rates are in use, the identified modal parameters will be the true modal parameters of the underlying linear system.

### Generalization (GCVD)


## Examples
Find the provided examples in the Example folder to use the CVD/GCVD.

## Cite the work
```
Add bibliography information here
```
