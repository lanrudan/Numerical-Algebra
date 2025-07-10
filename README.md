# Numerical-Algebra
# 数值代数
### 算法来源于徐树方、高立、张平文等编著的《数值线性代数（第二版）》，北京师范大学出版社。
$\bullet$ [fwd_sub](fwd_sub.m)：**前代法**解下三角方程组。  
$\bullet$ [back_sub](back_sub.m)：**回代法**解上三角方程组。  
$\bullet$ [naive_gaussian_elimination](naive_gaussian_elimination.m)：**朴素高斯消元法**，用于求解线性方程组。  
$\bullet$ [gaussian_elimination_complete_pivoting](gaussian_elimination_complete_pivoting.m)：**全主元高斯消元法**，用于对矩阵进行带有行、列置换的LU分解。  
$\bullet$ [gaussian_elimination_pivoting](gaussian_elimination_pivoting.m)：**列主元高斯消元法**，用于求解一般线性方程组。  
$\bullet$ [cholesky_solve](cholesky_solve.m)：**Cholesky分解算法**，用于高效地求解对称正定线性方程组。  
$\bullet$ [ldlt_decomposition](ldlt_decomposition.m)：**$LDL^T$ 分解（LDL Transpose Decomposition）算法**，对对称正定矩阵进行分解，并利用分解结果来高效地求解线性方程组 $Ax=b$。  
$\bullet$ [hager_norm_estimator](hager_norm_estimator.m)：用于迭代估计**矩阵逆的1-范数**。  
$\bullet$ [householder_transform](householder_transform.m)：**Householder变换**，用于根据输入向量构造豪斯霍尔德矩阵。  
$\bullet$ [givens_rotation](givens_rotation.m)：**Givens旋转**的核心计算，用于根据两个数值生成旋转矩阵的余弦和正弦因子。  
$\bullet$ [qr_decomposition_householder](qr_decomposition_householder.m)：**基于Householder变换的QR分解**，并利用分解结果求解线性方程组。  
$\bullet$ [jacobi_iteration](jacobi_iteration.m)：**Jacobi迭代法**，用于迭代求解线性方程组。  
$\bullet$ [gauss_seidel_iteration](gauss_seidel_iteration.m)：**高斯-赛德尔迭代法**，用于迭代求解线性方程组。  
$\bullet$ [sor_iteration](sor_iteration.m)：**逐次超松弛 (SOR) 迭代法**，用于迭代求解线性方程组。  
$\bullet$ [steepest_descent_method](steepest_descent_method.m)：**最速下降法**，用于迭代求解对称正定线性方程组。  
$\bullet$ [conjugate_gradient_method](conjugate_gradient_method.m)：**共轭梯度法 (CG)**，用于迭代求解对称正定线性方程组。  
$\bullet$ [conjugate_gradient_solver_practical](conjugate_gradient_solver_practical.m)：**实用共轭梯度法**，用于迭代求解对称正定线性方程组，并采用相对误差作为收敛准则。  
$\bullet$ [preconditioned_conjugate_gradient](preconditioned_conjugate_gradient.m)：**预处理共轭梯度法 (PCG)**，通过使用对角预处理器迭代求解对称正定线性方程组。  
$\bullet$ [power_method](power_method.m)：**幂法**，用于估计矩阵的最大模特征值及其对应的特征向量。  
$\bullet$ [inverse_power_method](inverse_power_method.m)：**反幂法**，通过对矩阵的逆应用幂法，用于估计矩阵的最小模特征值及其对应的特征向量。  
$\bullet$ [hessenberg_decomposition_householder](hessenberg_decomposition_householder.m)：**基于Householder变换的上Hessenberg分解**，将一般方阵相似变换为上Hessenberg矩阵。  
$\bullet$ [double_shift_qr_iteration](double_shift_qr_iteration.m)：**双重步长位移的QR迭代**，用于计算矩阵（通常为上Hessenberg矩阵）的特征值，并通过相似变换将其进一步简化。  
$\bullet$ [implicit_qr](implicit_qr.m)：**隐式QR算法**，通过迭代地对上Hessenberg矩阵应用双重步长QR变换和矩阵分裂技术，以计算其特征值。  
$\bullet$ [symmetric_tridiagonalization_householder](symmetric_tridiagonalization_householder.m)：**基于Householder变换的对称矩阵三对角分解**，将对称矩阵相似变换为三对角矩阵。  
$\bullet$ [implicit_symmetric_qr_algorithm](implicit_symmetric_qr_algorithm.m)：**隐式对称QR算法**，通过迭代地对三对角矩阵应用双重步长QR变换和分裂技术，以计算其特征值。  
$\bullet$ [implicit_symmetric_qr_wilkinson](implicit_symmetric_qr_wilkinson.m)：**带Wilkinson位移的隐式对称QR迭代**，通过一系列Givens旋转将三对角矩阵进一步简化，以加速其特征值的收敛。




