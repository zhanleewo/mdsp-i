<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第六讲-线性预测编码" style="text-align: center; margin-bottom: 2rem; border-bottom: none;">第六讲 线性预测编码</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

## 1. 背景知识

我们研究的信号被认为在一个较短的时间段内具有**宽平稳性**，记为 $\{X(n)\}$。

### 1.1 宽平稳随机过程

宽平稳随机过程满足：均值恒定，自相关函数仅依赖于时间差。即
$$
R_X(n, m) = \mathbb{E}[X(n) X(m)] = R_X(n-m),  \tag{6.1}$$
其中 $R_X(k)$ 是自相关函数，$k = n-m$。这意味着信号的二阶统计量不随时间平移而改变。

### 1.2 语音编码

在语音编码中，利用信号的短时平稳性：在一段时间内，我们可以用一个固定参数的模型（例如线性预测模型）来描述信号。只需要在编码端传输一次模型参数（系数），而不传输原始样本，从而大大减少数据传输量，达到压缩目的。只有当模型即将失效（即不能再假设平稳）时，才重新传输一组新的参数。

### 1.3 预测模型

给定过去 $k$ 个观测值 $X(n-1), X(n-2), \dots, X(n-k)$，我们想预测当前值 $X(n)$。这称为 **$k$ 阶预测**。

### 1.4 线性预测

假设预测值是过去 $k$ 个值的线性组合： $$
\hat{X}(n) = \sum_{i=1}^k \alpha_i X(n-i).  $$
记预测系数向量 $$\alpha = (\alpha_1, \alpha_2, \dots, \alpha_k)^\top \tag{6.2}，$$观测向量
$$
X(n-1:n-k) = \big( X(n-1), X(n-2), \dots, X(n-k) \big)^\top.  \tag{6.3}$$
则预测误差为 $$
e(n) = X(n) - \hat{X}(n) = X(n) - \alpha^\top X(n-1:n-k).  \tag{6.4}$$

定义均方预测误差
$$
h(\alpha) = \mathbb{E}\big[ e(n)^2 \big] = \mathbb{E}\big[ (X(n) - \alpha^\top X(n-1:n-k))^2 \big].  \tag{6.5}$$
我们的目标是选择 $\alpha$ 使 $h(\alpha)$ 最小化。

**梯度推导**：展开 $h(\alpha)$： $$
h(\alpha) = \mathbb{E}[X(n)^2] - 2\alpha^\top \mathbb{E}\big[ X(n-1:n-k) X(n) \big] + \alpha^\top \mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big] \alpha.  \tag{6.6}$$
对 $\alpha$ 求梯度，并令梯度为零，得到正则方程（Wiener-Hopf 方程）：
$$
\mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big] \, \alpha = \mathbb{E}\big[ X(n) X(n-1:n-k) \big].  \tag{6.7}$$
记 $$
R = \mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big], \quad r = \mathbb{E}\big[ X(n) X(n-1:n-k) \big],  \tag{6.8}$$
则上述方程为 $R \alpha = r$。若 $R$ 正定，则最优系数为
$$
\alpha_{\text{opt}} = R^{-1} r = \mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big]^{-1} \mathbb{E}\big[ X(n) X(n-1:n-k) \big].  \tag{6.9}$$
此时最小均方预测误差为 $$
h_{\min} = \mathbb{E}[X(n)^2] - r^\top R^{-1} r.  \tag{6.10}$$

> **注**：这里的符号 $X(n-1:n-k)$ 表示列向量，$X^\top$ 表示行向量，完全遵循我们的约定。梯度表达式中出现的 $\mathbb{E}[X(n-1:n-k) X^\top(n-1:n-k)]$ 正是 $R$，而 $\mathbb{E}[X(n) X(n-1:n-k)]$ 正是 $r$。上述推导与维纳滤波完全一致，只是将目标变量换成了 $X(n)$，观测向量换成了过去 $k$ 个样本。在实际计算中，$R$ 是 Toeplitz 矩阵，可以通过 Levinson‑Durbin 递推高效求解，而不必直接求逆。

如果能够对 $R$ 直接求逆，那么我们就能立即得到 $\alpha_{\text{opt}} = R^{-1} r$，从而一步确定最优预测系数。然而，在实际的硬件实现中（如嵌入式系统、数字信号处理器或低功耗芯片），矩阵求逆的计算复杂度极高。具体来说：
- 普通矩阵求逆的复杂度为 $O(k^3)$；
- 若利用矩阵的对称性（$R$ 对称正定），可降至 $O(k^2)$；
- 进一步利用 $R$ 的 **Toeplitz** 结构（每条对角线元素相同），则求逆复杂度可降至 $O(k)$。

当预测阶数 $k$ 取 10~20 时，$O(k^3)$ 与 $O(k)$ 的差异可达两个数量级。在实时语音编码中（如 8kHz 采样率，每帧 10~20ms），频繁求解 $R\alpha = r$ 若使用普通求逆将严重占用计算资源，增加功耗和延迟，难以满足实时处理的需求。因此，我们必须充分利用 Toeplitz 矩阵的特殊结构，设计更高效的求解算法，避免显式求逆。

下一节我们将详细推导 **Levinson-Durbin 递推**，该算法利用 Toeplitz 矩阵的规律，通过递推逐阶更新预测系数，以 $O(k^2)$ 的复杂度（比 $O(k)$ 稍高，但远低于 $O(k^3)$）高效求得 $\alpha_{\text{opt}}$，为语音编码等实时应用提供了可行的工程方案。

---

## 2. 算法推导

求解过程大致分为四个步骤:
- **齐次化**: 
- **增广**:
- **递推**:
- **Toeplitz**:

下面我们逐步展开这些步骤，并最终给出完整的 Levinson-Durbin 递归算法。

### 2.1 Wiener-Hopf 方程写成矩阵形式

根据前一节的推导，最优预测系数 $\alpha$ 满足正则方程：
$$
\mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big] \, \alpha = \mathbb{E}\big[ X(n) X(n-1:n-k) \big].  \tag{6.11}$$
将矩阵显式写出。设信号的自相关函数为 $R_X(m) = \mathbb{E}[X(n) X(n-m)]$，由宽平稳性，$R_X(-m)=R_X(m)$。

观测向量 $X(n-1:n-k) = (X(n-1), X(n-2), \dots, X(n-k))^\top$，则自相关矩阵为： $$ 
\begin{aligned}
   R_k &= \mathbb{E}\big[ X(n-1:n-k) X^\top(n-1:n-k) \big] \\
   &= \begin{pmatrix}
    R_X(0) & R_X(1) & R_X(2) & \cdots & R_X(k-1) \\
    R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-2) \\
    R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-3) \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    R_X(k-1) & R_X(k-2) & R_X(k-3) & \cdots & R_X(0)
    \end{pmatrix}
\end{aligned}
  \tag{6.12}$$ 

互相关向量为： $$ 
\begin{aligned}
r_k &= \mathbb{E}\big[ X(n) X(n-1:n-k) \big] \\
&=
\begin{pmatrix}
R_X(1) \\ R_X(2) \\ R_X(3) \\ \vdots \\ R_X(k)
\end{pmatrix}
\end{aligned}
  \tag{6.13}$$

预测系数向量 $\alpha = (\alpha_1, \alpha_2, \dots, \alpha_k)^\top$。于是方程 $R_k \alpha = r_k$ 写为：

$$
\begin{pmatrix}
R_X(0) & R_X(1) & \cdots & R_X(k-1) \\
R_X(1) & R_X(0) & \cdots & R_X(k-2) \\
\vdots & \vdots & \ddots & \vdots \\
R_X(k-1) & R_X(k-2) & \cdots & R_X(0)
\end{pmatrix}
\begin{pmatrix}
\alpha_{k,1} \\ \alpha_{k,2} \\ \vdots \\ \alpha_{k,k}
\end{pmatrix}
=
\begin{pmatrix}
R_X(1) \\ R_X(2) \\ \vdots \\ R_X(k)
\end{pmatrix}
  \tag{6.14}$$

这就是 **Yule‑Walker 方程**的标准矩阵形式。利用 $R_k$ 的 Toeplitz 结构，可以设计更高效的递推算法，求解过程大致分为以下四个步骤：

- **齐次化**：引入反向预测问题。定义反向预测系数 $\beta_k = (\beta_{k,1}, \beta_{k,2}, \dots, \beta_{k,k})^\top$ 满足 $R_k \beta_k = \tilde{r}_k$，其中 $\tilde{r}_k = (R_X(k), R_X(k-1), \dots, R_X(1))^\top$。通过同时求解前向（原问题）和反向问题，建立两者之间的对称关系，为递推准备。

- **增广**：将 $k$ 阶方程扩展到 $k+1$ 阶。利用 Toeplitz 矩阵的嵌套结构，将 $R_{k+1}$ 和 $r_{k+1}$ 用低阶矩阵和新增元素表示： $$
  R_{k+1} = \begin{pmatrix} R_k & \tilde{r}_k \\ \tilde{r}_k^\top & R_X(0) \end{pmatrix},
  \qquad
  r_{k+1} = \begin{pmatrix} r_k \\ R_X(k+1) \end{pmatrix}.
    \tag{6.15}$$
  利用低阶解 $\alpha_k$ 和 $\beta_k$ 作为初始猜测，构造高阶解的待定形式。

- **递推**：通过反射系数（PARCOR）$\gamma_k$ 递推更新预测系数。递推公式为： $$
  \alpha_{k+1} = \begin{pmatrix} \alpha_k \\ 0 \end{pmatrix} + \gamma_k \begin{pmatrix} -\beta_k \\ 1 \end{pmatrix},
    \tag{6.16}$$
  其中反射系数 $\gamma_k$ 由预测误差功率 $E_k$ 和互相关新息确定： $$
  \gamma_k = \frac{R_X(k+1) - \tilde{r}_k^\top \alpha_k}{E_k},
    \tag{6.17}$$
  前向预测误差功率递推：$E_{k+1} = E_k (1 - \gamma_k^2)$。同时反向预测系数 $\beta_{k+1}$ 也可类似更新。每步计算量为 $O(k)$，总复杂度 $O(k^2)$。

- **Toeplitz**：整个算法成功的关键在于 $R_k$ 的 **Toeplitz 结构**（每条对角线元素相同）。这一性质保证了增广矩阵具有简单的分块形式，使得低阶解可以线性组合成高阶解，并且反射系数能够通过标量计算得到。如果没有 Toeplitz 性质，递推关系将不成立，算法复杂度无法降低。

通过以上四步，Levinson‑Durbin 递推以 $O(k^2)$ 的复杂度高效求解 Yule‑Walker 方程，在语音编码、谱估计等实时处理中得到广泛应用。

### 2.2 齐次化

首先将原方程 (6.14) 改写为齐次形式。将方程右边的 $r_k$ 移到左边，并引入一个额外的 $-1$ 分量，得到如下的齐次方程： $$
\begin{pmatrix} 
R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-1) \\
R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-2) \\
\vdots & \vdots & \ddots & \vdots \\
R_X(k) & R_X(k-1) & R_X(k-2) & \cdots & R_X(0)
\end{pmatrix}
\begin{pmatrix}
-1 \\ \alpha_{k,1} \\ \alpha_{k,2} \\ \vdots \\ \alpha_{k,k}
\end{pmatrix}
= 0
  \tag{6.18}$$
这个齐次方程将原问题与反向预测联系起来，为后续的增广步骤做好了准备。

### 2.3 增广

观察 (6.18)，左边的矩阵是 $k \times (k+1)$ 的（非方阵）。为了得到方阵以便求解，我们在上面增加一行（对应自相关函数 $R_X(0)$ 到 $R_X(k)$），同时将右边补上一个非零向量，从而得到： $$
\begin{pmatrix}
R_X(0) & R_X(1) & R_X(2) & \cdots & R_X(k) \\
R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-1) \\
R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-2) \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
R_X(k) & R_X(k-1) & R_X(k-2) & \cdots & R_X(0)
\end{pmatrix}
\begin{pmatrix}
-1 \\ \alpha_{k,1} \\ \alpha_{k,2} \\ \vdots \\ \alpha_{k,k}
\end{pmatrix} \\
= \begin{pmatrix}
-E_k \\ 0 \\ 0 \\ \vdots \\ 0 
\end{pmatrix}
  \tag{6.19}$$
现在左边矩阵变成了 $(k+1) \times (k+1)$ 的方阵，右边第一行是 $-E_k$，其余为 $0$。这里 $E_k$ 是 $k$ 阶前向预测误差功率，其具体表达式将在后文给出。

### 2.4 递推

假设我们已经知道 $k-1$ 阶的解 $\alpha_{k-1} = (\alpha_{k-1,1}, \dots, \alpha_{k-1,k-1})^\top$。那么 $k$ 阶的解可以通过递推得到，其形式为： $$
\alpha_{k} = \begin{pmatrix} \alpha_{k-1} \\ 0 \end{pmatrix} + \gamma_{k-1} \begin{pmatrix} -\beta_{k-1} \\ 1 \end{pmatrix},
  \tag{6.20}$$
其中 $\beta_{k-1}$ 是 $k-1$ 阶的反向预测系数，$\gamma_{k-1}$ 是反射系数。为了验证这一递推，我们写出 $k$ 阶增广方程的低阶部分。考虑将 $k$ 阶增广矩阵作用于一个试探向量，该向量由 $k-1$ 阶解补零构成，结果如下： $$ 
\begin{pmatrix} 
R_X(0) & R_X(1) & R_X(2) & \cdots & R_X(k-1) & R_X(k) \\
R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-2) & R_X(k-1) \\
R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-3) & R_X(k-2) \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
R_X(k-1) & R_X(k-2) & R_X(k-3) & \cdots & R_X(2) & R_X(1) \\
R_X(k) & R_X(k-2) & R_X(k-3) & \cdots & R_X(1) & R_X(0) 
\end{pmatrix}
\begin{pmatrix}
-1 \\ \alpha_{k-1,1} \\ \alpha_{k-1,2} \\ \vdots \\ \alpha_{k-1,k-1} \\ 0
\end{pmatrix} \\
= \begin{pmatrix}
-E_{k-1} \\ 0 \\ 0 \\ \vdots \\ 0 \\ A_k
\end{pmatrix}
  \tag{6.21}$$
这里 $A_k$ 是一个非零量（因为补零导致最后一行的方程不满足）。$E_k$ 和 $A_k$ 的表达式为： $$
E_k = R_X(0) - \sum_{i=1}^{k} \alpha_{k,i} R_{i},
  \tag{6.22}$$
 $$
A_k = -R_X(k) + \sum_{i=1}^{k-1} \alpha_{k-1,i} R_{k-i}.
  \tag{6.23}$$
多出来的 $A_k$ 怎么办？我们需要通过调整最后一行的方程来消除它。这就要用到 Toeplitz 矩阵的反转性质。

### 2.5 利用 Toeplitz 性质

首先给出 Toeplitz 矩阵的一个重要性质：若 $R$ 是 Toeplitz 矩阵，且满足 $R x = y$，则 $R \tilde{x} = \tilde{y}$，其中 $\tilde{x}$ 表示向量的反转，即 $\tilde{x} = (x_n, \dots, x_1)$，$\tilde{y}$ 类似。证明如下：

设 $R$ 是 Toeplitz 矩阵，$R_{ij} = R_{|i-j|}$。则 $(Rx)_i = \sum_j R(i-j) x_j = y_i$。将下标反转，令 $i' = n+1-i$，$j' = n+1-j$，则 $$
(R \tilde{x})_{i'} = \sum_j R(i'-j') x_{n+1-j} = \sum_j R((n+1-i)-(n+1-j)) x_{n+1-j} = \sum_j R(j-i) x_{n+1-j}.
  \tag{6.24}$$
由于 $R(j-i)=R(i-j)$（自相关函数的偶对称性），并且由原方程，$\sum_j R(i-j) x_j = y_i$，经变量替换可得 $(R \tilde{x})_{i'} = y_{n+1-i'} = \tilde{y}_{i'}$，因此 $R \tilde{x} = \tilde{y}$。

有了这个性质，我们可以从 (6.21) 推导出另一个方程，只需将左边的待乘向量和右边的结果向量同时反转： $$
\begin{pmatrix}
R_X(0) & R_X(1) & R_X(2) & \cdots & R_X(k-1) & R_X(k) \\
R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-2) & R_X(k-1) \\
R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-3) & R_X(k-2) \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
R_X(k-1) & R_X(k-2) & R_X(k-3) & \cdots & R_X(2) & R_X(1) \\
R_X(k) & R_X(k-2) & R_X(k-3) & \cdots & R_X(1) & R_X(0) 
\end{pmatrix}
\begin{pmatrix}
0 \\ \alpha_{k-1,k-1} \\ \alpha_{k-1,k-2} \\ \vdots \\ \alpha_{k-1,1} \\ -1
\end{pmatrix} \\
= \begin{pmatrix}
A_k \\ 0 \\ 0 \\ \vdots \\ 0 \\ -E_{k-1}
\end{pmatrix}
  \tag{6.25}$$

现在，我们将 (6.21) 和 (6.25) 中的两个解进行线性组合，以消除最后一行的非零项。取组合系数 $\rho_k$，构造新的向量： $$
\left(
\begin{pmatrix}
-1 \\ \alpha_{k-1,1} \\ \alpha_{k-1,2} \\ \vdots \\ \alpha_{k-1,k-1} \\ 0
\end{pmatrix}
+
\rho_k
\begin{pmatrix}
0 \\ \alpha_{k-1,k-1} \\ \alpha_{k-1,k-2} \\ \vdots \\ \alpha_{k-1,1} \\ -1
\end{pmatrix}
\right)
  \tag{6.26}$$
乘以$R_k$ 后，右边为： $$
\begin{pmatrix}
-E_{k-1} \\ 0 \\ 0 \\ \vdots \\ 0 \\ A_k
\end{pmatrix}
+\rho_k
\begin{pmatrix}
A_k \\ 0 \\ 0 \\ \vdots \\ 0 \\ -E_{k-1}
\end{pmatrix}
  \tag{6.27}$$

将左边向量合并，得到： $$ 
\begin{pmatrix} 
R_X(0) & R_X(1) & R_X(2) & \cdots & R_X(k-1) & R_X(k) \\
R_X(1) & R_X(0) & R_X(1) & \cdots & R_X(k-2) & R_X(k-1) \\
R_X(2) & R_X(1) & R_X(0) & \cdots & R_X(k-3) & R_X(k-2) \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
R_X(k-1) & R_X(k-2) & R_X(k-3) & \cdots & R_X(2) & R_X(1) \\
R_X(k) & R_X(k-2) & R_X(k-3) & \cdots & R_X(1) & R_X(0) 
\end{pmatrix}
\times
\begin{pmatrix}
-1 \\ \alpha_{k-1,1} + \rho_k \alpha_{k-1,k-1} \\ \alpha_{k-1,2} + \rho_k \alpha_{k-1,k-2} \\ \vdots \\ \alpha_{k-1,k-1} + \rho_k \alpha_{k-1,1} \\ -\rho_k
\end{pmatrix} \\
= \begin{pmatrix}
-E_{k-1} + \rho_k A_k \\ 0 \\ 0 \\ \vdots \\ 0 \\ A_k - \rho_k E_{k-1}
\end{pmatrix}
  \tag{6.28}$$

为了使得最后一行（第 $k$ 行）变为 $0$，我们选择 $\rho_k$ 满足： $$
A_k - \rho_k E_{k-1} = 0 \quad \Longrightarrow \quad \rho_k = \frac{A_k}{E_{k-1}}.
  \tag{6.29}$$
此时右边第一行变为： $$
-E_{k-1} + \rho_k A_k = -E_{k-1} + \frac{A_k^2}{E_{k-1}} = -\frac{E_{k-1}^2 - A_k^2}{E_{k-1}} = -E_k,
  \tag{6.30}$$
其中 $E_k$ 定义为 $E_k = E_{k-1} - \frac{A_k^2}{E_{k-1}} = (1-\rho_k^2) E_{k-1}$。

于是，我们得到了 $k$ 阶预测系数的递推公式： $$
\alpha_{k, i} = \alpha_{k-1, i} + \rho_k \alpha_{k-1, k-i}, \quad i = 1, \dots, k-1,
  \tag{6.31}$$
 $$
\alpha_{k, k} = -\rho_k.
  \tag{6.32}$$

至此，我们完成了 Levinson-Durbin 递推的推导。下面总结算法步骤。

### 2.6 Levinson-Durbin 算法总结

1. **初始化（$k=1$）**：  
   直接解一阶 Yule-Walker 方程： $$ 
   \begin{pmatrix} 
   R_X(0) & R_X(1) \\
   R_X(1) & R_X(0) 
   \end{pmatrix}
   \begin{pmatrix}
   -1 \\ \alpha_{1,1} 
   \end{pmatrix}
   = \begin{pmatrix}
   -E_1 \\ 0
   \end{pmatrix}
     \tag{6.33}$$
   解得： $$
   \alpha_{1,1} = \frac{R_X(1)}{R_X(0)},\qquad E_1 = R_X(0) - \alpha_{1,1} R_X(1) = R_X(0) - \frac{R_X^2(1)}{R_X(0)}.
     \tag{6.34}$$

2.**递推（$k = 2, 3, \dots$）**：  
   - 已知 $k-1$ 阶的系数 $\alpha_{k-1,1},\dots,\alpha_{k-1,k-1}$ 以及预测误差功率 $E_{k-1}$。
   - 计算 $A_k = -R_X(k) + \sum_{i=1}^{k-1} \alpha_{k-1,i} R_X(k-i)$。
   - 计算反射系数 $\rho_k = A_k / E_{k-1}$。
   - 更新预测误差功率 $E_k = (1 - \rho_k^2) E_{k-1}$。
   - 计算 $k$ 阶预测系数： 
  
  $$
  \alpha_{k,i} = \alpha_{k-1,i} + \rho_k \alpha_{k-1,k-i},\quad i = 1,\dots,k-1,
       \tag{6.35}
  $$
       
  $$
  \alpha_{k,k} = -\rho_k. \tag{6.36}
  $$

这样，我们从一阶开始，逐步递推到所需阶数 $p$，即可得到最优线性预测系数。该算法利用了 Toeplitz 矩阵的结构，计算复杂度为 $O(p^2)$，远低于直接矩阵求逆的 $O(p^3)$，是实时语音编码中的标准方法。

---

## 3. 前向-后向预测：反转的统计意涵

上一节我们从线性代数与矩阵论的角度看到了 Levinson‑Durbin 算法中“反转”操作的关键作用：利用 Toeplitz 矩阵的反转不变性，将低阶解线性组合得到高阶解。这一操作看似是纯粹的代数技巧，实际上蕴含着深刻的统计意义——它对应着前向预测与后向预测之间的对偶关系。

### 3.1 预测模型

给定过去 $k$ 个观测值 $X(n-1), X(n-2), \dots, X(n-k)$，我们想预测当前值 $X(n)$。这称为 **$k$ 阶预测**。

### 3.2 前向预测与后向预测

从这里我们就引入递推思想：  
给定过去 $k-1$ 个观测值 $X(n-1), X(n-2), \dots, X(n-k+1)$，我们想预测当前值 $X(n)$。这称为 **$k-1$ 阶预测**。

预测可不可以看成是投影？即 $$
\operatorname{Proj}_{X(n-1:n-k)} X_n = \operatorname{Proj}_{X(n-1:n-k+1)} X(n) + \operatorname{Proj}_{X(n-k)} X(n),
  \tag{6.37}$$
但在前面的 Wiener 滤波学习中我们知道，这在一般情况下需要正交性，即 $X(n-k)$ 与 $X(n-1:n-k+1)$ 正交。因此，我们需要对观测空间进行正交化分解： $$
X(n-1:n-k) = X(n-1:n-k+1) \oplus \left(X(n-k) - \underbrace{\operatorname{Proj}_{X(n-1:n-k+1)} X(n-k)}_{\text{后向预测}\;\text{的一维部分}}\right).
  \tag{6.38}$$
于是，预测可以写作： $$
\operatorname{Proj}_{X(n-1:n-k)} X(n) = \operatorname{Proj}_{X(n-1:n-k+1)} X(n) + \rho_k \Bigl(X(n-k) - \underbrace{\operatorname{Proj}_{X(n-1:n-k+1)} X(n-k)}_{\text{后向预测误差}}\Bigr).
  \tag{6.39}$$
其中 $\rho_k$ 是投影系数，反映新息 $X(n-k)$ 对预测 $X(n)$ 的贡献。

前向预测与后向预测在这里同时出现：  
- 前向预测：用 $X(n-1),\dots,X(n-k+1)$ 预测 $X(n)$；  
- 后向预测：用 $X(n-1),\dots,X(n-k+1)$ 预测 $X(n-k)$。

正交化后得到后向预测的误差 $\epsilon_k = X(n-k) - \operatorname{Proj}_{X(n-1:n-k+1)} X(n-k)$，然后利用这个误差提升前向预测的结果。

这个预测的本质是假设预测值是过去 $k$ 个值的线性组合： $$
\begin{aligned}
\hat{X}(n) &= \sum_{i=1}^k \alpha_{k,i} X(n-i) \\
&= \sum_{i=1}^{k-1}\alpha_{k-1, i} X(n-i) + \rho_k X(n-k) -\rho_k \sum_{i=1}^{k-1}\beta_{k-1,i} X(n-i) \\
&= \sum_{i=1}^{k-1} \bigl(\alpha_{k-1, i} - \rho_k\beta_{k-1,i}\bigr) X(n-i) + \rho_k X(n-k).
\end{aligned}
  \tag{6.40}$$
其中 $\beta_{k-1,i}$ 是 $k-1$ 阶后向预测系数。因此，如果能清楚认识 $\beta$，我们将会对这个过程有更深刻的统计认识。

反向预测系数 $\beta_k$ 满足如下 Yule-Walker 方程： 

$$
\begin{pmatrix}
R_X(0) & R_X(1) & \cdots & R_X(k-1) \\
R_X(1) & R_X(0) & \cdots & R_X(k-2) \\
\vdots & \vdots & \ddots & \vdots \\
R_X(k-1) & R_X(k-2) & \cdots & R_X(0)
\end{pmatrix}
\begin{pmatrix}
\beta_{k,1} \\ \beta_{k,2} \\ \vdots \\ \beta_{k,k}
\end{pmatrix}
=
\begin{pmatrix}
R_X(k-1) \\ R_X(k-2) \\ \vdots \\ R_X(1)
\end{pmatrix}
  \tag{6.41}$$

由 Toeplitz 结构的对称性，可以证明反向系数与前向系数反转后的关系： 

$$
\begin{pmatrix}
\beta_{k,1} \\ \beta_{k,2} \\ \vdots \\ \beta_{k,k}
\end{pmatrix}
=
\begin{pmatrix}
\alpha_{k,k} \\ \alpha_{k,k-1} \\ \vdots \\ \alpha_{k,1}
\end{pmatrix}
  \tag{6.42}$$
这正是“反转即后向预测”的统计含义。

接下来，我们还需要搞清楚 $\rho_k$ 的统计意义。定义后向预测残差： $$
b_k = X(n-k) - \operatorname{Proj}_{X(n-1:n-k+1)} X(n-k) = X(n-k) - \sum_{i=1}^{k-1} \beta_{k-1,i} X(n-i).
  \tag{6.43}$$
根据正交投影原理，$b_k$ 与 $X(n-1),\dots,X(n-k+1)$ 正交。那么投影系数 $\rho_k$ 即为 $X(n)$ 在 $b_k$ 上的投影： $$
\rho_k = \frac{\mathbb{E}[X(n) b_k]}{\mathbb{E}[b_k^2]}.
  \tag{6.44}$$
计算分子： $$
\mathbb{E}[X(n) b_k] = \mathbb{E}\left[ X(n) \left( X(n-k) - \sum_{i=1}^{k-1} \beta_{k-1,i} X(n-i) \right) \right]
= R_X(k) - \sum_{i=1}^{k-1} \beta_{k-1,i} R_X(i).
  \tag{6.45}$$
由 (3.2.1)，$\beta_{k-1,i} = \alpha_{k-1,k-i}$，因此 $$
\mathbb{E}[X(n) b_k] = R_X(k) - \sum_{i=1}^{k-1} \alpha_{k-1,k-i} R_X(i).
  \tag{6.46}$$
这与前面定义中的 $A_k$（差一个符号）密切相关：实际上 $A_k = -R_X(k) + \sum_{i=1}^{k-1} \alpha_{k-1,i} R_X(k-i)$，通过变量替换 $i' = k-i$ 可知 $A_k = -\mathbb{E}[X(n) b_k]$。所以 $A_k$ 就是数据 $X(n)$ 与残差 $b_k$ 之间的互相关（取反）。

计算分母： $$
\mathbb{E}[b_k^2] = \mathbb{E}\left[ \left( X(n-k) - \sum_{i=1}^{k-1} \beta_{k-1,i} X(n-i) \right)^2 \right].
  \tag{6.47}$$
由于最优估计的残差一定正交于观测数据，因此 $$
\mathbb{E}[b_k^2] = \mathbb{E}\left[ X(n-k) \left( X(n-k) - \sum_{i=1}^{k-1} \beta_{k-1,i} X(n-i) \right) \right] = R_X(0) - \sum_{i=1}^{k-1} \beta_{k-1,i} R_X(k-i).
  \tag{6.48}$$
利用反转关系，$\beta_{k-1,i} = \alpha_{k-1,k-i}$，则 $$
\mathbb{E}[b_k^2] = R_X(0) - \sum_{i=1}^{k-1} \alpha_{k-1,k-i} R_X(k-i).
  \tag{6.49}$$
而前面我们已知前向预测误差功率 $E_{k-1} = R_X(0) - \sum_{i=1}^{k-1} \alpha_{k-1,i} R_X(i)$，通过变量替换 $i' = k-i$ 可知两者相等，即 $\mathbb{E}[b_k^2] = E_{k-1}$。因此，$E_{k-1}$ 实际上是后向预测残差的能量。

综上，反射系数 $$
\rho_k = \frac{\mathbb{E}[X(n) b_k]}{\mathbb{E}[b_k^2]} = \frac{-A_k}{E_{k-1}}.
  \tag{6.50}$$
而我们的递推中取 $\rho_k = A_k / E_{k-1}$，两者仅差一个符号。通常习惯将反射系数定义为 $\rho_k = \frac{A_k}{E_{k-1}}$，这并不影响算法的本质。这里 $A_k$ 的符号取决于定义，实际递推时通过 $A_k = -R_X(k) + \sum \alpha_{k-1,i} R_X(k-i)$ 计算即可。

从统计观点看，$\rho_k$ 就是 $X(n)$ 与后向预测误差 $b_k$ 之间的相关系数，它度量了在剔除掉低阶预测信息后，$X(n)$ 与 $X(n-k)$ 之间的偏相关程度。**反射系数 $\rho_k$ 的绝对值不超过 1，且 $|\rho_k|$ 越小表示第 $k$ 步相关性越弱**。

反转操作在统计上对应着后向预测，即利用未来预测过去。通过正交化分解，前向预测可以分解为低阶前向预测加上一个修正项，这个修正项正比于后向预测误差。这正是 Levinson-Durbin 递推的统计本质，也是该算法能够高效求解 Yule-Walker 方程的原因。

## 4. 课后总结

本节对线性预测编码（LPC）的核心内容和关键结论进行梳理，帮助快速回顾整篇文章的要点。

---

### 4.1 宽平稳性与线性预测模型

- **宽平稳性假设**：信号 $\{X(n)\}$ 的自相关函数 $R_X(m) = \mathbb{E}[X(n)X(n-m)]$ 仅依赖于时间差 $m$，不依赖于绝对时间 $n$。
- **语音编码动机**：利用短时平稳性，用固定模型参数（预测系数）代替传输原始样本，实现数据压缩。
- **预测模型**：用过去 $k$ 个样本 $X(n-1),\dots,X(n-k)$ 线性预测当前值 $X(n)$： $$
  \hat{X}(n) = \sum_{i=1}^k \alpha_i X(n-i)
    \tag{6.51}$$
- **均方预测误差**：$h(\alpha) = \mathbb{E}[(X(n) - \alpha^\top X_{n-1})^2]$，最小化得到 Yule‑Walker 方程 $R_k \alpha = r_k$。
- **最优系数**：$\alpha_{\text{opt}} = R_k^{-1} r_k$，最小预测误差 $h_{\min} = \mathbb{E}[X(n)^2] - r_k^\top R_k^{-1} r_k$。

---

### 4.2 复杂度分析与 Toeplitz 结构

- **矩阵求逆复杂度**：
  - 普通矩阵：$O(k^3)$
  - 对称矩阵：$O(k^2)$
  - Toeplitz 矩阵：$O(k)$（理论可逆，但工程更常用递推）
- **Toeplitz 矩阵**：每条对角线元素相同，如 $R_k$ 中 $[R_k]_{ij} = R_X(i-j)$。
- **利用 Toeplitz 结构**：设计 Levinson-Durbin 递推，将求解复杂度降至 $O(k^2)$，同时避免显式求逆。

---

### 4.3 Levinson‑Durbin 算法核心步骤

1. **齐次化**：将原方程改写为齐次形式，引入反向预测。
2. **增广**：将 $k$ 阶方程扩展为 $(k+1)$ 阶，得到方阵。
3. **递推**：利用反射系数 $\rho_k$ 从 $k-1$ 阶系数递推 $k$ 阶系数： $$
   \alpha_{k,i} = \alpha_{k-1,i} + \rho_k \alpha_{k-1,k-i},\quad i=1,\dots,k-1,\qquad \alpha_{k,k} = -\rho_k
     \tag{6.52}$$
4. **Toeplitz 性质**：利用反转不变性 $R \tilde{x} = \tilde{y}$，证明递推公式的正确性。

- **反射系数** $\rho_k = A_k / E_{k-1}$，其中 $A_k = -R_X(k) + \sum_{i=1}^{k-1} \alpha_{k-1,i} R_X(k-i)$，$E_{k-1}$ 是前向预测误差功率。
- **误差功率递推**：$E_k = (1-\rho_k^2) E_{k-1}$，且 $E_k$ 单调不增。

---

### 4.4 统计视角：前向‑后向预测

- **正交化分解**：将观测空间分解为低阶观测子空间与正交新息分量（后向预测误差）。
- **前向预测** = 低阶前向预测 + 反射系数 × 后向预测误差。
- **反转与后向预测**：反向预测系数 $\beta_k$ 正好是前向预测系数 $\alpha_k$ 的反转，即 $\beta_{k,i} = \alpha_{k,k+1-i}$。
- **反射系数的统计意义**：$\rho_k$ 是 $X(n)$ 与后向预测误差 $b_k$ 之间的偏相关系数，反映第 $k$ 步新息对预测的贡献，且 $|\rho_k| \le 1$。

---

### 4.5 算法流程总结

1. **输入**：自相关函数 $R_X(0), R_X(1), \dots, R_X(p)$。
2. **初始化**：$\alpha_{1,1} = R_X(1)/R_X(0)$，$E_1 = R_X(0) - \alpha_{1,1} R_X(1)$。
3. **循环** $k = 2$ 到 $p$：
   - 计算 $A_k = -R_X(k) + \sum_{i=1}^{k-1} \alpha_{k-1,i} R_X(k-i)$
   - $\rho_k = A_k / E_{k-1}$
   - $E_k = (1-\rho_k^2) E_{k-1}$
   - 对 $i=1$ 到 $k-1$：$\alpha_{k,i} = \alpha_{k-1,i} + \rho_k \alpha_{k-1,k-i}$
   - $\alpha_{k,k} = -\rho_k$
4. **输出**：$p$ 阶预测系数 $\alpha_{p,1},\dots,\alpha_{p,p}$ 和最终预测误差功率 $E_p$。

---

### 4.6 应用与意义

- **语音编码**：低比特率压缩（如 GSM、VoIP）的核心技术。
- **谱估计**：由 LPC 系数构成全极点滤波器，逼近信号功率谱。
- **实时性**：Levinson-Durbin 算法的高效性（$O(p^2)$）使其适合嵌入式实时实现。
- **理论价值**：统一了线性预测、维纳滤波、正交投影、Toeplitz 矩阵和反射系数的数学与统计概念。

---

### 4.7 学习检查清单

- [ ] 能写出线性预测的基本公式 $\hat{X}(n) = -\sum_{i=1}^{p} \alpha_i X(n-i)$，并说明其与 AR 模型的等价性
- [ ] 能写出 Yule-Walker 方程的矩阵形式，并识别其 Toeplitz 结构
- [ ] 能解释正交性原理在线性预测中的作用：预测误差与过去观测正交
- [ ] 能描述 Levinson-Durbin 递推的核心步骤（齐次化、增广、递推），并说明反射系数的统计意义
- [ ] 能写出前向预测误差和后向预测误差的关系，并说明反射系数 $\rho_k$ 是偏相关系数
- [ ] 能解释为什么误差功率的递推公式为 $E_k = (1-\rho_k^2) E_{k-1}$，且单调不增
- [ ] 能说明为什么 LPC 在语音编码中如此重要：全极点模型逼近声道滤波器
- [ ] 能对比 LPC 谱估计与周期图谱估计的区别

### 4.8 思考题

1. **Levinson-Durbin 的"魔法"**：为什么 Toeplitz 结构能让求解复杂度从 $O(p^3)$ 降至 $O(p^2)$？Toeplitz 矩阵的"反转不变性"在递推中起了什么关键作用？如果自相关矩阵不是 Toeplitz 的（比如用有偏和无偏自相关估计的区别），Levinson-Durbin 还能用吗？

2. **反射系数与稳定性**：Burg 算法通过反射系数 $|\rho_k| \le 1$ 保证 AR 滤波器的稳定性。$|\rho_k| = 1$ 意味着什么？在物理上，什么信号会产生 $|\rho_k| \approx 1$ 的反射系数？

3. **前向与后向预测的双重角色**：为什么同时考虑前向和后向预测能改善估计？（提示：这相当于在数据两端都利用了可用信息，Burg 算法正是基于这一点。）这与"数据增强"的思想有何联系？

4. **LPC 的工程取舍**：LPC 将语音信号压缩为少数几个反射系数，实现了高压缩比。但这种全极点模型在表示鼻音（含有零点）时会遇到什么困难？为什么工程中仍然选择全极点模型而非零极点模型？

5. **从 LPC 到谱估计**：用 LPC 系数构造全极点谱 $S(\omega) = \sigma^2 / |1 + \sum \alpha_i e^{-j\omega i}|^2$，得到的谱有什么特点？它为什么能产生尖锐的共振峰？如果阶数选高了会出现什么问题？


<div style="page-break-before: always;"></div>