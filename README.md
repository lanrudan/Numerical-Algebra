# Numerical-Algebra
# 数值代数
记录矩阵计算代码  
$\bullet$ [fwd_sub](fwd_sub.m)：前代法解**下三角方程组**。  
$\bullet$ [back_sub](back_sub.m)：回代法解**上三角方程组**。  
$\bullet$ [naive_gaussian_elimination](naive_gaussian_elimination.m)：朴素高斯消元法，用于求解线性方程组。  
$\bullet$ [gaussian_elimination_complete_pivoting](gaussian_elimination_complete_pivoting.m)：全主元高斯消元法，用于对矩阵进行带有行、列置换的**LU分解**。  
$\bullet$ [gaussian_elimination_pivoting](gaussian_elimination_pivoting.m)：列主元高斯消元法，用于求解**一般线性方程组**。  
$\bullet$ [cholesky_solve](cholesky_solve.m)：Cholesky分解算法，用于高效地求解**对称正定线性方程组**。  
$\bullet$ [ldlt_decomposition](ldlt_decomposition.m)： $LDL^T$ 分解（LDL Transpose Decomposition）算法，对**对称正定矩阵**进行分解，并利用分解结果来高效地求解线性方程组 $Ax=b$。  
$\bullet$ [householder_transform](householder_transform.m)：**Householder变换**，用于根据输入向量构造豪斯霍尔德矩阵。  
$\bullet$ [gauss_seidel_iteration](gauss_seidel_iteration.m)：**高斯-赛德尔迭代法**，用于迭代求解线性方程组。
