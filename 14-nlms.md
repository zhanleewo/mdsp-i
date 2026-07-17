<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第十四讲-归一化最小均方" style="text-align: center; margin-bottom: 2rem; border-bottom: none;">第十四讲 归一化最小均方</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

## 1. NLMS 回顾

在上一篇文章中,我们从LMS的尺度敏感性问题出发,一步步走到了NLMS(归一化最小均方，有些场景下也叫归一化最小二乘)。现在,我们把那条线索**原封不动地拽回来**,因为它是NLMS的起点。

### 1.1 标准 LMS 的步长敏感性

标准LMS的迭代步骤我们早就刻在脑子里了：

$$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu X(n) e(n)  \tag{14.1}$$

其中误差$e(n)$定义为： $$
e(n) = d(n) - X^T(n) \theta_n  \tag{14.2}$$

这里$\hat{\theta}_n$是$n$时刻的参数估计(列向量),$X(n)$是$n$时刻的输入向量(列向量),$X^T(n)$是它的转置(行向量),$\mu$是步长。

这个公式看起来干净利落,但它藏着一个**致命的工程隐患**：它对输入信号的幅度极其敏感。

### 1.2 尺度变化：四倍放大的灾难

我们上一篇文章已经算过这笔账了，现在再走一遍，因为这是NLMS诞生的直接导火索。

假设输入信号$X(n)$的幅度整体扩大两倍。为了保持同样的物理模型(即最优参数$\theta$不变),期望输出$d(n)$也必须相应地扩大两倍,因为$d(n) \approx X^T(n)\theta$。

于是：

- 输入向量变为：$2X(n)$
- 期望输出变为：$2d(n)$
- 瞬时误差变为：

$$
e'(n) = 2d(n) - (2X(n))^T \theta_n = 2\left(d(n) - X^T(n)\theta_n\right) = 2e(n)  \tag{14.3}$$

现在,把变化后的输入和误差代进LMS更新公式,看看更新量变成了什么： $$
X'(n) e'(n) = (2X(n)) \cdot (2e(n)) = 4 X(n) e(n)  \tag{14.4}$$

**同一个参数估计问题,$\theta$的最优解纹丝未动,但LMS的更新量却硬生生被放大了四倍。**

这就相当于：你本来调好了一个步长$\mu$,觉得挺稳,结果输入信号稍微大了一点,有效步长就变成了$4\mu$——轻则收敛变慢,重则直接发散。如果你为了防发散把$\mu$调小,那在信号幅度小的时候又慢得像乌龟爬。

这就是LMS尺度敏感性的根源。

### 1.3 工程归一化处理及其局限

面对这种“幅度一变就崩”的局面,工程师们想出了一个直观到几乎粗暴的办法：**既然更新量$X(n)e(n)$被放大了,那我就在更新之前除以输入向量的能量,把它压回去。**

这就是归一化LMS的工程形式：

$$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.5}$$

其中误差$e(n)$照旧： $$
e(n) = d(n) - X^T(n) \theta_n  \tag{14.6}$$

注意这个归一化因子的写法：$\frac{X(n)}{\|X(n)\|^2}$——分子是$X(n)$,分母是$\|X(n)\|^2$(注意这是范数平方,是一个标量)。整个分式仍然是一个$n$维列向量。

验证一下尺度不变性：如果$X(n)$和$d(n)$同时扩大两倍,则$e(n)$扩大两倍,而$X(n)/\|X(n)\|^2$缩小两倍(分子乘2,分母乘4)。两者相乘后更新量**完全不变**。

这就是我们想要的：**对输入尺度免疫**。

但这个归一化因子是人为加上的,它确实解决了工程问题,但它到底有没有理论基础？还是说它只是一个聪明的“补丁”？在Goodwin给出解释之前,答案就是：**纯属瞎搞,调参玄学**。

### 1.4 Goodwin 的洞见：将 NLMS 转化为约束优化

Goodwin等人当年干了一件漂亮的事：他们设计了一个优化问题,然后证明NLMS的更新公式恰恰是这个问题**精确闭式解**。这一下就把NLMS从“工程修补”抬进了“理论殿堂”。

优化问题长这样：

$$
\min_{\hat{\theta}_{n+1}} \|\hat{\theta}_{n+1} - \hat{\theta}_{n}\|^2, \quad \text{s.t.} \quad d(n) = X^T(n) \hat{\theta}_{n+1}  \tag{14.7}$$

这个问题的含义是：**在保证新模型能够完美拟合当前数据点的前提下,参数的变化量尽可能小**。这是典型的信任域思想——我不信任大跨步,每一步都走得小心翼翼。

令$\Delta \hat{\theta}_{n+1} = \hat{\theta}_{n+1} - \hat{\theta}_{n}$,则$\hat{\theta}_{n+1} = \hat{\theta}_n + \Delta \hat{\theta}_{n+1}$。代入约束条件： $$
d(n) = X^T(n) (\hat{\theta}_n + \Delta \hat{\theta}_{n+1})  \tag{14.8}$$

展开：

$$
d(n) = X^T(n) \hat{\theta}_n + X^T(n) \Delta \hat{\theta}_{n+1}  \tag{14.9}$$

移项： $$
X^T(n) \Delta \hat{\theta}_{n+1} = d(n) - X^T(n) \hat{\theta}_n = e(n)  \tag{14.10}$$

于是整个优化问题变成了：

$$
\min_{\Delta \hat{\theta}_{n+1}} \|\Delta \hat{\theta}_{n+1}\|^2, \quad \text{s.t.} \quad e(n) = X^T(n) \Delta \hat{\theta}_{n+1}  \tag{14.11}$$

这就是一个**单等式约束的欠定最小二乘问题**——方程的个数(1个)少于未知数的个数($n$个)。我们上一篇文章已经推导过这种问题的闭式解了： $$
\Delta \hat{\theta}_{n+1} = X(n) \left( X^T(n) X(n) \right)^{-1} e(n) = \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.12}$$

把$\Delta \hat{\theta}_{n+1} = \hat{\theta}_{n+1} - \hat{\theta}_n$代回去：

$$
\hat{\theta}_{n+1} - \hat{\theta}_n = \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.13}$$

移项： $$
\hat{\theta}_{n+1} = \hat{\theta}_n + \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.14}$$

**这就是NLMS的更新公式。**

### 1.5 回顾的意义

通过上面的回顾,我们看清了几件事：

1. **LMS的尺度敏感性**不是小毛病,是结构性的——输入输出的幅度变化会直接放大有效步长。
2. **工程归一化**虽然解决了尺度问题,但一开始缺乏理论支撑。
3. **Goodwin的约束优化**把NLMS从”调参玄学”变成了”带信任域约束的最小范数解”，并且印证了我们在欠定最小二乘中建立的伪逆理论——NLMS的更新方向正是$X^T(n)$的伪逆作用于误差$e(n)$的结果。

但回顾归回顾,NLMS的故事到这里还远没有结束。我们在上一篇文章中只解决了两个问题：“为什么LMS不行”和“为什么NLMS有理论依据”。但接下来的问题才是真正棘手的：

- **收敛速度**：NLMS比LMS快了,但具体快多少？它的收敛速度受什么因素影响？
- **稳态误差**：NLMS收敛后,系数还在抖动,这个抖动有多大？与步长$\mu$是什么关系？
- **步长选择**：$\mu$到底怎么选？有没有理论上的最优范围？
- **与RLS的对比**：NLMS和RLS都在处理输入相关性,但RLS为什么收敛更快？NLMS在什么情况下够用,什么情况下必须上RLS？

这些问题,就是我们这篇文章要啃的硬骨头。

---


## 2. 角度一：在最小扰动下精确拟合当前数据

### 2.1 问题的提出

上面的回顾就是我们的第一个角度，从最直观的场景出发。假设我们有一个旧模型——参数为\(\hat{\theta}_n\)，它过去的表现还行。现在来了一个新样本\((X(n), d(n))\)。

这个我们做个直观的比喻：**让高年级的做低年级的题，要在全部能做对的前提下，年级差异最小的，这是一个角度。**

我们希望在**不把模型改得面目全非**的前提下，让新模型\(\hat{\theta}_{n+1}\)在$X(n)$ 上**完全正确**。

什么叫“完全正确”？就是新模型的输出与期望输出完全一致，即误差为零：

$$
d(n) = X^T(n) \hat{\theta}_{n+1}  \tag{14.15}$$

什么叫“不改得面目全非”？就是新旧参数之间的变化量\(\|\hat{\theta}_{n+1} - \hat{\theta}_n\|\)尽可能小。

这就是一个约束优化问题： $$
\min_{\hat{\theta}_{n+1}} \|\hat{\theta}_{n+1} - \hat{\theta}_{n}\|^2, \quad \text{s.t.} \quad d(n) = X^T(n) \hat{\theta}_{n+1}  \tag{14.16}$$

---

## 3. 角度二：从"最大化差距"视角理解 NLMS

上一节我们用了“高年级做低年级题”的比喻,讲的是**在确保所有题都做对的前提下,年级差异最小**。换句话说就是：**求稳，不折腾**。

现在我们把硬币翻过来。另一个同样合理、甚至更直观的思路是：

### 3.1 重新表述：差距越大越好

还是那个场景。低年级的学生(旧的参数$\hat{\theta}_n$)和高年级的学生(新的参数$\hat{\theta}_{n+1}$)都来做同一套题(当前的数据$X(n)$),看看谁的答案更接近标准答案$d(n)$。

我们定义：

- **低年级的误差**(旧模型的误差)：

$$
e_1(n) = d(n) - X^T(n) \hat{\theta}_n  \tag{14.17}$$

- **高年级的误差**(新模型的误差)： $$
e_2(n) = d(n) - X^T(n) \hat{\theta}_{n+1}  \tag{14.18}$$

现在的问题是：高年级学生的目标是什么？站在高年级的角度——**我不能比低年级差,而且我赢你赢得越多越好**。

也就是说,我们希望让两个误差之间的**差距尽可能大**,同时保证高年级至少不比低年级差：

$$
\max \; |e_2(n) - e_1(n)|, \quad \text{s.t.} \quad |e_2(n)| < |e_1(n)|  \tag{14.19}$$

这个约束$|e_2(n)| < |e_1(n)|$的意思是：**高年级的错误率必须严格低于低年级**,即新模型在当前数据点上的表现必须比旧模型更好。在这个前提下,我们尽可能地“甩开”差距。

### 3.2 展开差距表达式

两个误差的差是什么？ $$
e_2(n) - e_1(n) = \left(d(n) - X^T(n) \hat{\theta}_{n+1}\right) - \left(d(n) - X^T(n) \hat{\theta}_n\right)  \tag{14.20}$$

中间$d(n)$直接抵消,剩下：

$$
e_2(n) - e_1(n) = - X^T(n) \hat{\theta}_{n+1} + X^T(n) \hat{\theta}_n = - X^T(n) \left( \hat{\theta}_{n+1} - \hat{\theta}_n \right)  \tag{14.21}$$

令$\Delta \hat{\theta}_{n+1} = \hat{\theta}_{n+1} - \hat{\theta}_n$,则： $$
e_2(n) - e_1(n) = - X^T(n) \Delta \hat{\theta}_{n+1}  \tag{14.22}$$

于是我们的最大化问题变成了：

$$
\max \; |e_2(n) - e_1(n)| = \max \; |X^T(n) \Delta \hat{\theta}_{n+1}|  \tag{14.23}$$

### 3.3 应用柯西-施瓦茨不等式

现在我们要求$|X^T(n) \Delta \hat{\theta}_{n+1}|$的最大值。这和$\Delta \hat{\theta}_{n+1}$的方向有关——它不是随便取的,因为我们还要受$|e_2(n)| < |e_1(n)|$这个约束的限制。

不过我们先看无约束情况下怎么取到最大值。由柯西不等式： $$
|X^T(n) \Delta \hat{\theta}_{n+1}| \leq \|X(n)\| \cdot \|\Delta \hat{\theta}_{n+1}\|  \tag{14.24}$$

等号成立的条件是：**$X(n)$与$\Delta \hat{\theta}_{n+1}$平行**,即两者之间的夹角为0,相关系数$\cos \theta = 1$。

也就是说,要让两个误差的差距最大,参数的变化量$\Delta \hat{\theta}_{n+1}$必须沿着输入向量$X(n)$的方向。

因此我们可以写出：

$$
\Delta \hat{\theta}_{n+1} = k X(n)  \tag{14.25}$$

其中$k$是一个待定的标量系数,可正可负。

### 3.4 确定系数$k$

把$\Delta \hat{\theta}_{n+1} = k X(n)$代入误差差的表达式： $$
e_2(n) - e_1(n) = - X^T(n) \cdot k X(n) = - k \|X(n)\|^2  \tag{14.26}$$

另一方面,我们要求$|e_2(n)| < |e_1(n)|$。我们试试把$e_2(n)$写成$e_1(n)$的某个倍数：

令$e_2(n) = (1 - \mu) e_1(n)$,其中$\mu \in (0, 1)$。这样：

- 当$\mu = 0$时,$e_2 = e_1$,高年级和低年级一样,没有任何进步；
- 当$\mu \to 1$时,$e_2 \to 0$,高年级接近完美；
- $\mu$越接近1,差距越大。

于是误差差为：

$$
e_2(n) - e_1(n) = (1-\mu)e_1(n) - e_1(n) = -\mu e_1(n)  \tag{14.27}$$

取绝对值： $$
|e_2(n) - e_1(n)| = \mu |e_1(n)|  \tag{14.28}$$

但我们在上一步已经知道$|e_2(n) - e_1(n)| = |k| \|X(n)\|^2$。因此：

$$
|k| \|X(n)\|^2 = \mu |e_1(n)|  \tag{14.29}$$

于是： $$
|k| = \mu \frac{|e_1(n)|}{\|X(n)\|^2}  \tag{14.30}$$

所以$k$可以写成一个带符号的形式：

$$
k = \operatorname{sgn}(k) \cdot \mu \frac{|e_1(n)|}{\|X(n)\|^2}  \tag{14.31}$$

再注意到$\operatorname{sgn}(k)$与$\operatorname{sgn}(e_1(n))$可以合并为一个新的符号因子$\tilde{\mu}$(因为我们后面会把这两者的乘积吸收进步长参数里),于是可以更紧凑地写成： $$
k = \tilde{\mu} \frac{e_1(n)}{\|X(n)\|^2}  \tag{14.32}$$

其中$\tilde{\mu} = \operatorname{sgn}(k) \cdot \operatorname{sgn}(e_1(n)) \cdot \mu$。

### 3.5 回到参数更新

由于$\Delta \hat{\theta}_{n+1} = k X(n)$,把上面求出的$k$代进去：

$$
\Delta \hat{\theta}_{n+1} = \tilde{\mu} \frac{X(n)}{\|X(n)\|^2} e_1(n)  \tag{14.33}$$

而$\Delta \hat{\theta}_{n+1} = \hat{\theta}_{n+1} - \hat{\theta}_n$,因此： $$
\hat{\theta}_{n+1} - \hat{\theta}_n = \tilde{\mu} \frac{X(n)}{\|X(n)\|^2} e_1(n)  \tag{14.34}$$

移项：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \tilde{\mu} \frac{X(n)}{\|X(n)\|^2} e_1(n)  \tag{14.35}$$

**这就是NLMS的更新公式。** 其中$\tilde{\mu}$就是我们通常写的步长参数$\mu$(区别在于这里的$\tilde{\mu}$被限制在$(0,1)$之间,而标准NLMS中$\mu$通常取在$(0,2)$之间,这里的符号吸收稍微有点差异,但本质一致)。

### 3.6 两个角度的统一

现在我们有了**两个角度**来看NLMS：

| 角度 | 问题描述 | 结果 | 直观含义 |
|---|---|---|---|
| **角度一** (Goodwin) | $\min \|\Delta \hat{\theta}\|^2$, s.t. $e_2 = 0$ | $\Delta \hat{\theta} = \frac{X}{\|X\|^2} e_1$ | **求稳**：在能彻底消灭误差的方向上,步子迈得最小 |
| **角度二** (本节) | $\max \|e_2 - e_1\|$, s.t. $\|e_2\| < \|e_1\|$ | $\Delta \hat{\theta} = \mu \frac{X}{\|X\|^2} e_1$ | **进取**：在不退步的前提下,一步跨得最远 |

两个角度得到相同结果。

角度一强调的是**约束下的最小变化**——既然旧模型已经不错了,就别大动干戈,能解决问题就行。角度二强调的是**约束下的最大进步**——既然必须变,那就变得最有意义,在保证不倒退的前提下,一步到位。

两者是同一枚硬币的两面。恰好对应到步长$\mu$这个参数上：

- $\mu$趋近于0：更偏向角度一——变化极小,几乎是原地踏步;
- $\mu$趋近于1：更偏向角度二——变化极大,恨不得一步跨到最优解。

最妙的是,这两个看似对立的“好”标准,最终推出了**同一个公式**,只是步长$\mu$在第一个角度里是硬解出来的(等于1),在第二个角度里是自由参数(在$(0,1)$之间)。步长这个参数,恰恰就是在这两种哲学之间的滑动调节器。

---
## 4. 从正交投影的角度理解NLMS

前面我们用了两个角度来认识NLMS：角度一是“最小扰动+精确拟合”的约束优化，角度二是“最大化差距”的极值问题。现在我们从第三个角度——**正交投影**——来重新审视这个问题。

这个角度比前两个更接近LMS和RLS的“家族血统”，因为我们一直在说：**最小二乘就是做正交投影。** NLMS也不例外。

### 4.1 定义残差

设$\hat{\theta}_{\text{opt}}$是上帝视角的最优参数——即我们永远无法直接得到、但理论上存在的那个完美系数。定义第$n$时刻的系数残差为： $$
\epsilon(n) = \hat{\theta}_n - \hat{\theta}_{\text{opt}}  \tag{14.36}$$

这个$\epsilon(n)$衡量的是：**当前模型离完美模型还有多远。** 如果$\epsilon(n) = 0$，我们就已经到达终点了。

### 4.2 误差与残差的关系

第$n$时刻的瞬时误差定义为：

$$
e(n) = d(n) - X^T(n) \hat{\theta}_n  \tag{14.37}$$

把$d(n)$用最优参数表示出来。根据之前学过的内容，我们很容易得出：最优估计满足$d(n) = X^T(n) \hat{\theta}_{\text{opt}}$（在无噪声的理想情况下）。于是： $$
e(n) = X^T(n) \hat{\theta}_{\text{opt}} - X^T(n) \hat{\theta}_n
= X^T(n) \left( \hat{\theta}_{\text{opt}} - \hat{\theta}_n \right)
= - X^T(n) \epsilon(n)  \tag{14.38}$$

或者写成：

$$
e(n) = X^T(n) \epsilon(n)  \tag{14.39}$$

（符号取决于$\epsilon(n)$的定义方式，这里我们用$\epsilon(n) = \hat{\theta}_n - \hat{\theta}_{\text{opt}}$，所以$e(n) = -X^T(n)\epsilon(n)$。但$\epsilon(n)$可以重新定义为$\hat{\theta}_{\text{opt}} - \hat{\theta}_n$来消掉负号，不影响本质。）

关键结论：**瞬时误差$e(n)$就是残差$\epsilon(n)$在输入方向$X(n)$上的投影（乘以范数）**。

换句话说，$e(n)$告诉我们：**当前的系数误差$\epsilon(n)$还有多少分量沿着$X(n)$的方向没有被消除。**

### 4.3 我们希望残差与数据正交

如果我们能做到$e(n) = 0$，那就意味着$X^T(n)\epsilon(n) = 0$——即**残差与输入数据正交**。

这显然是理想状态：残差与所有观测方向正交，意味着系数误差在任何一个数据方向上都没有分量，那就是真·最优解。

但实际情况是，我们只有一个样本$X(n)$，残差$\epsilon(n)$在这个方向上大概率有分量。所以我们要做的，就是把$\epsilon(n)$中沿着$X(n)$方向的那个分量**减掉**。这就是正交化的思想——将$\epsilon(n)$与当前数据$X(n)$正交化。

### 4.4 逐次正交化：Gram-Schmidt的精神

我们的目标是从$\epsilon(1)$开始，每来一个样本$X(n)$，就从当前的残差中减去它在$X(n)$方向上的投影分量： $$
\epsilon(2) = \epsilon(1) - \text{Proj}_{X(1)} \epsilon(1)  \tag{14.40}$$

$$
\epsilon(3) = \epsilon(2) - \text{Proj}_{X(2)} \epsilon(2)  \tag{14.41}$$ 

$$
\vdots  \tag{14.42}
$$

$$
\epsilon(n+1) = \epsilon(n) - \text{Proj}_{X(n)} \epsilon(n)  \tag{14.43}$$

这个投影算子的表达式是（一维投影，因为$X(n)$是一个列向量，张成一维子空间）： $$
\text{Proj}_{X(n)} \epsilon(n) = X(n) \left( X^T(n) X(n) \right)^{-1} X^T(n) \epsilon(n)  \tag{14.44}$$

因为$X^T(n) X(n) = \|X(n)\|^2$是一个标量，所以：

$$
\text{Proj}_{X(n)} \epsilon(n) = \frac{X(n) X^T(n)}{\|X(n)\|^2} \epsilon(n)  \tag{14.45}$$

代入递推式： $$
\epsilon(n+1) = \epsilon(n) - \frac{X(n) X^T(n)}{\|X(n)\|^2} \epsilon(n)  \tag{14.46}$$

### 4.5 换回$\hat{\theta}$的表示

用$\epsilon(n) = \hat{\theta}_n - \hat{\theta}_{\text{opt}}$代回去。左边：

$$
\epsilon(n+1) = \hat{\theta}_{n+1} - \hat{\theta}_{\text{opt}}  \tag{14.47}$$

右边第一项： $$
\epsilon(n) = \hat{\theta}_n - \hat{\theta}_{\text{opt}}  \tag{14.48}$$

右边第二项的分子$X(n) X^T(n) \epsilon(n)$：

$$
X^T(n) \epsilon(n) = X^T(n) \left( \hat{\theta}_n - \hat{\theta}_{\text{opt}} \right)
= X^T(n) \hat{\theta}_n - X^T(n) \hat{\theta}_{\text{opt}}
= - \left( d(n) - X^T(n) \hat{\theta}_n \right) = - e(n)  \tag{14.49}$$

其中$e(n) = d(n) - X^T(n)\hat{\theta}_n$是瞬时误差。

于是： $$
\frac{X(n) X^T(n)}{\|X(n)\|^2} \epsilon(n) = \frac{X(n)}{\|X(n)\|^2} \cdot \left( X^T(n) \epsilon(n) \right)
= \frac{X(n)}{\|X(n)\|^2} \cdot \left( - e(n) \right)
= - \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.50}$$

代回$\epsilon(n+1)$的递推式：

$$
\hat{\theta}_{n+1} - \hat{\theta}_{\text{opt}} = \left( \hat{\theta}_n - \hat{\theta}_{\text{opt}} \right) - \left( - \frac{X(n)}{\|X(n)\|^2} e(n) \right)  \tag{14.51}$$

即： $$
\hat{\theta}_{n+1} - \hat{\theta}_{\text{opt}} = \hat{\theta}_n - \hat{\theta}_{\text{opt}} + \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.52}$$

两边同时消去$\hat{\theta}_{\text{opt}}$：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.53}$$

**这就是NLMS的更新公式。**

### 4.6 这个角度告诉了我们什么

正交投影这个角度揭示了NLMS的本质：**每一步都在从当前的系数误差中减去它在当前输入方向上的投影。**

这就好比一个逐次正交化过程：

- 旧参数$\hat{\theta}_n$离最优$\hat{\theta}_{\text{opt}}$还差一个$\epsilon(n)$；
- 当前数据$X(n)$“探测”到了这个残差的一部分——就是瞬时误差$e(n)$；
- 我们把这个被探测到的部分从残差中**减掉**（也就是从旧参数中减去），得到新参数；
- 于是在$X(n)$这个方向上，残差被清零了。

如果输入序列的各方向足够丰富（即$X(1), X(2), \ldots$在各个方向上都有分量），那么经过足够多步之后，残差在各方向上的分量都会被依次消除，最终$\epsilon(n) \to 0$，$\hat{\theta}_n \to \hat{\theta}_{\text{opt}}$。

这就是**NLMS收敛性的几何解释**：它不是一个“随机近似”，而是一个**逐次正交化过程**。每一步都是精确的（在单样本意义下），只是我们不知道$\hat{\theta}_{\text{opt}}$，所以用当前残差在$X(n)$方向上的投影来表示。

这也解释了为什么NLMS的收敛比LMS快得多——LMS是在误差曲面上做随机梯度下降，而NLMS直接对残差本身做正交投影（在单样本子空间上），方向更准，步子更稳。

### 4.7 三个角度的统一

至此，我们从三个完全不同的起点出发，得到了同一个公式：

| 角度 | 起点 | 方法 | 结果 |
|---|---|---|---|
| **角度一** | 约束优化 | 最小化$\|\Delta \hat{\theta}\|^2$，约束$e_2(n)=0$ | $\Delta \hat{\theta} = \frac{X}{\|X\|^2} e$ |
| **角度二** | 极值问题 | 最大化$\|e_2 - e_1\|$，约束$\|e_2\|<\|e_1\|$ | $\Delta \hat{\theta} = \mu \frac{X}{\|X\|^2} e$ |
| **角度三** | 正交投影 | 从$\epsilon$中减去在$X$方向上的投影分量 | $\Delta \hat{\theta} = \frac{X}{\|X\|^2} e$ |

三个角度得到相同结果，说明NLMS的更新公式是**结构性的必然**，而不是某个人的灵感或修补。

## 5. 用牛顿法解释NLMS的更新公式对角加载

下面就是NLMS的更新公式： $$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.54}$$

这个更新公式有个问题，$\|X(n)\|^2$在分母上，当$\|X(n)\|$趋近于0，或者等于0时就很麻烦，所以在真正的操作上，我们用的是下面这个公式：
$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \frac{X(n)}{\epsilon + \|X(n)\|^2} e(n)  \tag{14.55}$$
在分母加了一个数，我们前面的章节也介绍过这实际上是对角加载，也可以理解为正则化。

今天我们要用牛顿法和正则化这两个工具来说明理论上这样做是有道理的。

### 5.1 牛顿法视角下的 NLMS 推导

在前面几节中，我们从三个不同的角度（约束优化、极值问题、正交投影）得到了NLMS的更新公式。现在我们换一个完全不同的视角——**牛顿法**（Newton's Method）。

牛顿法是优化理论中比梯度下降更高级的工具：它不仅仅利用目标函数的一阶信息（梯度），还利用了二阶信息（Hessian矩阵）。在误差曲面的曲率信息已知的情况下，牛顿法可以一步跳到二次模型的极小值点。

问题在于：当我们把牛顿法直接套用到NLMS的瞬时目标函数上时，会发生什么？

### 5.1.1 牛顿法的通用形式

考虑无约束优化问题： $$
\min_{\theta} f(\theta)  \tag{14.56}$$

假设当前在点$\hat{\theta}_n$处，我们对目标函数$f$做二阶泰勒展开（只保留到二次项）：

$$
f(\hat{\theta}_n + \Delta) = f(\hat{\theta}_n) + \nabla_{\theta}^T f(\hat{\theta}_n) \Delta + \Delta^T H(f(\hat{\theta}_n)) \Delta  \tag{14.57}$$

其中$\nabla_{\theta} f(\hat{\theta}_n)$是梯度向量，$H(f(\hat{\theta}_n))$是Hessian矩阵（二阶偏导数矩阵）。

为了找到这个二次近似模型的最小值点，我们对$\Delta$求梯度并令其为零： $$
\nabla_{\Delta} \left( \nabla_{\theta}^T f(\hat{\theta}_n) \Delta + \Delta^T H(f(\hat{\theta}_n)) \Delta \right)
= \nabla_{\theta} f(\hat{\theta}_n) + H(f(\hat{\theta}_n)) \Delta = 0  \tag{14.58}$$

如果Hessian矩阵$H(f(\hat{\theta}_n))$可逆，我们就可以解出最优步长方向：

$$
\Delta = - H(f(\hat{\theta}_n))^{-1} \nabla_{\theta} f(\hat{\theta}_n)  \tag{14.59}$$

于是牛顿法的迭代公式为： $$
\hat{\theta}_{n+1} = \hat{\theta}_n - H(f(\hat{\theta}_n))^{-1} \nabla_{\theta} f(\hat{\theta}_n)  \tag{14.60}$$

这就是牛顿法的标准形式。它的几何含义是：**在当前位置用二次曲面逼近目标函数，然后直接跳到这个二次曲面的谷底。**

### 5.1.2 把牛顿法套到瞬时最小二乘上

现在我们把牛顿法应用到NLMS的瞬时目标函数上。在NLMS中，每一时刻我们只关心当前样本的瞬时误差，因此目标函数为：

$$
f(\theta) = \|d(n) - X^T(n)\theta\|^2  \tag{14.61}$$

展开这个平方： $$
f(\theta) = |d(n)|^2 - 2d(n) X^T(n)\theta + \theta^T (X(n) X^T(n)) \theta  \tag{14.62}$$

注意：这里$X(n)$是列向量，$X^T(n)$是行向量，所以$X(n)X^T(n)$是一个$n \times n$的**秩一矩阵**。

**计算梯度**：对$\theta$求导

$$
\nabla_{\theta} f(\theta) = -2d(n) X(n) + 2X(n)X^T(n)\theta  \tag{14.63}$$

把$d(n) - X^T(n)\theta$提出来： $$
\nabla_{\theta} f(\theta) = -2 X(n) \left( d(n) - X^T(n)\theta \right)  \tag{14.64}$$

**计算Hessian矩阵**：对梯度再求一次导

$$
H(f(\theta)) = 2 X(n) X^T(n)  \tag{14.65}$$

### 5.1.3 问题：Hessian奇异

到这里问题就暴露了：$X(n)X^T(n)$是一个**秩一矩阵**。

如果$X(n)$的维度$m > 1$，那么$X(n)X^T(n)$的秩最多为1，必然奇异（不满秩）。这意味着Hessian矩阵$H(f(\theta)) = 2X(n)X^T(n)$**不可逆**。

Hessian不可逆，牛顿法的更新公式$\Delta = - H^{-1} \nabla f$就无法直接使用——因为$H^{-1}$根本不存在。

这就是为什么标准牛顿法不能直接套用到单样本瞬时误差上的根本原因：**单个样本提供的二阶曲率信息是退化的，它只在一个方向（$X(n)$方向）上有曲率，在其他方向上曲率为零。**

### 5.1.4 引入正则项：让Hessian可逆

既然Hessian奇异是因为$X(n)X^T(n)$秩一，那就给它加一个对角项$\epsilon I$来“撑”满秩。这就是**正则化**（Tikhonov正则化）的思想——在目标函数中加入一个惩罚项： $$
f(\theta) = \|d(n) - X^T(n)\theta\|^2 + \epsilon \|\theta\|^2  \tag{14.66}$$

这个额外的$\epsilon \|\theta\|^2$项鼓励参数向零收缩，相当于在Hessian矩阵上加了$\epsilon I$，使其变成$\epsilon I + X(n)X^T(n)$——这是一个满秩矩阵（只要$\epsilon > 0$）。

重新计算梯度和Hessian：

$$
\nabla_{\theta} f(\theta) = -2 X(n) \left( d(n) - X^T(n)\theta \right) + 2\epsilon \theta  \tag{14.67}$$ 

$$
H(f(\theta)) = 2 X(n)X^T(n) + 2\epsilon I  \tag{14.68}$$

现在Hessian可逆了。于是带正则项的牛顿法更新公式为：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu (\epsilon I + X(n)X^T(n))^{-1} \left( X(n) \underbrace{(d(n) - X^T(n)\hat{\theta}_n)}_{e(n)} - \epsilon \hat{\theta}_n \right)  \tag{14.69}$$

这里我们引入了一个步长$\mu$来控制更新的幅度（标准牛顿法$\mu=1$，但在自适应滤波中我们通常保留它作为调节参数）。

### 5.1.5 用矩阵求逆引理化简

现在我们用Woodbury矩阵求逆引理来化简$(\epsilon I + X(n)X^T(n))^{-1}$。

Woodbury恒等式： $$
(A + BCD)^{-1} = A^{-1} - A^{-1} B (DA^{-1}B + C^{-1})^{-1} D A^{-1}  \tag{14.70}$$

令$A = \epsilon I$，$B = X(n)$，$C = 1$，$D = X^T(n)$。则：

$$
(\epsilon I + X(n) X^T(n))^{-1} = \epsilon^{-1} I - \epsilon^{-1} X(n) \left( X^T(n) \epsilon^{-1} X(n) + 1 \right)^{-1} X^T(n) \epsilon^{-1}  \tag{14.71}$$

因为$X^T(n) \epsilon^{-1} X(n) = \frac{\|X(n)\|^2}{\epsilon}$，所以括号内的标量为： $$
1 + \frac{\|X(n)\|^2}{\epsilon} = \frac{\epsilon + \|X(n)\|^2}{\epsilon}  \tag{14.72}$$

它的倒数为$\frac{\epsilon}{\epsilon + \|X(n)\|^2}$。

于是：

$$
(\epsilon I + X(n) X^T(n))^{-1} = \frac{1}{\epsilon} I - \frac{1}{\epsilon^2} X(n) X^T(n) \cdot \frac{\epsilon}{\epsilon + \|X(n)\|^2}  \tag{14.73}$$

化简第二项：$\frac{1}{\epsilon^2} \cdot \epsilon = \frac{1}{\epsilon}$，所以： $$
(\epsilon I + X(n) X^T(n))^{-1} = \frac{1}{\epsilon} I - \frac{1}{\epsilon} \cdot \frac{X(n) X^T(n)}{\epsilon + \|X(n)\|^2}  \tag{14.74}$$

提公因式$\frac{1}{\epsilon}$：

$$
(\epsilon I + X(n) X^T(n))^{-1} = \frac{1}{\epsilon} \left( I - \frac{X(n) X^T(n)}{\epsilon + \|X(n)\|^2} \right)  \tag{14.75}$$

### 5.1.6 代回更新公式

现在我们计算$(\epsilon I + X(n)X^T(n))^{-1} X(n)$，因为更新公式中需要它乘以$e(n)$。

用上面得到的表达式： $$
(\epsilon I + X(n) X^T(n))^{-1} X(n) = \frac{1}{\epsilon} \left( I - \frac{X(n) X^T(n)}{\epsilon + \|X(n)\|^2} \right) X(n)  \tag{14.76}$$

展开：

$$
= \frac{1}{\epsilon} X(n) - \frac{1}{\epsilon} \cdot \frac{X(n) X^T(n) X(n)}{\epsilon + \|X(n)\|^2}  \tag{14.77}$$

因为$X^T(n) X(n) = \|X(n)\|^2$，所以$X(n) X^T(n) X(n) = \|X(n)\|^2 X(n)$： $$
= \frac{1}{\epsilon} X(n) - \frac{1}{\epsilon} \cdot \frac{\|X(n)\|^2 X(n)}{\epsilon + \|X(n)\|^2}  \tag{14.78}$$

提公因式$\frac{1}{\epsilon} X(n)$：

$$
= \frac{X(n)}{\epsilon} \left( 1 - \frac{\|X(n)\|^2}{\epsilon + \|X(n)\|^2} \right)
= \frac{X(n)}{\epsilon} \cdot \frac{\epsilon}{\epsilon + \|X(n)\|^2}
= \frac{X(n)}{\epsilon + \|X(n)\|^2}  \tag{14.79}$$

于是： $$
(\epsilon I + X(n)X^T(n))^{-1} X(n) = \frac{X(n)}{\epsilon + \|X(n)\|^2}  \tag{14.80}$$

### 5.1.7 最终结果

把上面化简的结果代回牛顿法更新公式：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \left[ \frac{X(n)}{\epsilon + \|X(n)\|^2} e(n) - (\epsilon I + X(n)X^T(n))^{-1} \epsilon \hat{\theta}_n \right]  \tag{14.81}$$

这就是从带正则项的牛顿法推导出来的更新公式。

### 5.1.8 结果的评价

说实话，**不漂亮**。

它比NLMS多了一项：$-\mu (\epsilon I + X(n)X^T(n))^{-1} \epsilon \hat{\theta}_n$。这一项是正则项$\epsilon \|\theta\|^2$带来的“副作用”——它把旧参数$\hat{\theta}_n$往回拉，方向是指向原点的收缩。

如果我们想要得到纯正的NLMS公式，就必须把这个多出来的项**扔掉**。问题是：凭什么扔？

做这个推导的人找了很多理由来解释为什么要扔掉这一项，比如：

- “$\epsilon$很小，所以这一项可以忽略”——但$\epsilon$小的话，$(\epsilon I + XX^T)^{-1}$中的$\epsilon I$又不足以解决Hessian的奇异性问题，自相矛盾。
- “这一项是正则化带来的偏差，我们在无偏估计中可以忽略”——但正则化本身就是引入偏差来换取数值稳定性，你说忽略就忽略？
- “在稳态附近$\hat{\theta}_n$变化缓慢，所以这一项可以视为常数”——这种近似在严格推导中站不住脚。

说白了，**这些理由都挺牵强的**。

这个推导路径的历史意义在于：它展示了牛顿法套到单样本最小二乘上会遇到什么问题，以及正则化如何“勉强”解决这个问题。但最终得到的公式并不比我们前几节通过约束优化和正交投影得到的结果更简洁或更深刻。

真正有价值的结论其实很简单：

1. **单样本最小二乘的Hessian矩阵是秩一的，不可逆。**
2. **加正则项可以强行让它可逆，但会多出一个收缩项。**
3. **扔掉这个收缩项，就回到了NLMS。**
4. **这说明NLMS本质上是一个“近似牛顿法”——它只保留了牛顿法中与数据方向$X(n)$相关的部分，丢弃了正则化带来的收缩部分。**

换句话说，NLMS可以看作是一个**截断的、正则化牛顿法的近似**。但这个近似能成立，需要$\epsilon$足够大来保证Hessian可逆，同时又足够小来让收缩项可忽略——这个“足够大又足够小”的区间是否存在，完全取决于具体数据，没有一个普适的保证。

所以，这条路径最终的价值不是得到一个漂亮的公式，而是**告诉我们为什么直接套牛顿法行不通，以及NLMS在牛顿法的框架下处于什么位置。** 它没有前三个角度那么干净利落，但也算是一个有意义的补充视角。

## 6. NLMS的推广：仿射投影自适应滤波算法（APA）

### 6.1 NLMS 直观理解回顾

**让高年级的做低年级的题，高年级要尽可能考的分数更高。**

这是单样本约束——我们只要求新模型在当前这一个样本上表现更好（或完全正确）。NLMS每一步只“盯”着当前这一个数据点，步子迈得再准，视野也有限。

### 6.2 推广：从单约束到多约束

现在这个推广是：**让高年级的做所有比自己低的年级的题，都要做得更好。**

也就是说，我们不再满足于新模型只拟合最近的一个样本，而是要求新模型在**过去若干个样本上同时表现更好**。我们要求新模型$\hat{\theta}_{n+1}$同时完美拟合过去$n-1$个样本。

约束条件变成了一组等式： $$
d(k) = X^T(k) \hat{\theta}_{n+1}, \quad k = 1, 2, \ldots, n-1  \tag{14.82}$$

把这$n-1$个约束写成矩阵形式。定义数据矩阵和期望输出向量：

$$
\tilde{X} = (X_1, X_2, \ldots, X_{n-1})  \tag{14.83}$$

注意：这里的$\tilde{X}$是一个$m \times (n-1)$矩阵，每一列$X_k$是一个$m$维列向量。$\tilde{d} = (d_1, d_2, \ldots, d_{n-1})^T$是一个$(n-1)$维列向量。

于是约束可以简洁地写成： $$
\tilde{d} = \tilde{X}^T \hat{\theta}_{n+1}  \tag{14.84}$$

我们希望新模型在满足所有这些约束的同时，与旧模型的差异尽可能小：

$$
\min_{\hat{\theta}_{n+1}} \| \hat{\theta}_{n+1} - \hat{\theta}_{n} \|^2, \quad \text{s.t.} \ \tilde{d} = \tilde{X}^T \hat{\theta}_{n+1}  \tag{14.85}$$

这就是**多约束版本的信任域优化问题**。下面按步骤推导出闭式解。

---

### 6.3 完整推导步骤

#### 6.3.1 步骤一：拉格朗日乘子法

令$\Delta \hat{\theta}_{n+1} = \hat{\theta}_{n+1} - \hat{\theta}_n$，则$\hat{\theta}_{n+1} = \hat{\theta}_n + \Delta \hat{\theta}_{n+1}$。代入约束： $$
\tilde{d} = \tilde{X}^T (\hat{\theta}_n + \Delta \hat{\theta}_{n+1})  \tag{14.86}$$

移项：

$$
\tilde{d} - \tilde{X}^T \hat{\theta}_n = \tilde{X}^T \Delta \hat{\theta}_{n+1}  \tag{14.87}$$

定义广义误差向量： $$
\tilde{e}(n) = \tilde{d} - \tilde{X}^T \hat{\theta}_n  \tag{14.88}$$

这是旧模型在$n-1$个历史样本上的误差向量（每一个分量对应一个历史样本的误差）。于是约束简化为：

$$
\tilde{e}(n) = \tilde{X}^T \Delta \hat{\theta}_{n+1}  \tag{14.89}$$

优化问题变为： $$
\min_{\Delta \hat{\theta}_{n+1}} \| \Delta \hat{\theta}_{n+1} \|^2, \quad \text{s.t.} \ \tilde{e}(n) = \tilde{X}^T \Delta \hat{\theta}_{n+1}  \tag{14.90}$$

现在引入拉格朗日乘子向量$\lambda \in \mathbb{R}^{n-1}$（注意：因为有$n-1$个等式约束，所以$\lambda$是一个$n-1$维列向量），构造拉格朗日函数：

$$
L(\Delta \hat{\theta}_{n+1}, \lambda) = \| \Delta \hat{\theta}_{n+1} \|^2 + \lambda^T \left( \tilde{e}(n) - \tilde{X}^T \Delta \hat{\theta}_{n+1} \right)  \tag{14.91}$$

展开： $$
L(\Delta \hat{\theta}_{n+1}, \lambda) = \Delta \hat{\theta}_{n+1}^T \Delta \hat{\theta}_{n+1} + \lambda^T \tilde{e}(n) - \lambda^T \tilde{X}^T \Delta \hat{\theta}_{n+1}  \tag{14.92}$$

注意第三项：$\lambda^T \tilde{X}^T \Delta \hat{\theta}_{n+1} = (\tilde{X} \lambda)^T \Delta \hat{\theta}_{n+1}$，这是一个标量。

#### 6.3.2 步骤二：求梯度，算出 $\Delta \hat{\theta}_{n+1}$

对$\Delta \hat{\theta}_{n+1}$求梯度（向量求导）：

$$
\nabla_{\Delta \hat{\theta}} L = 2 \Delta \hat{\theta}_{n+1} - \tilde{X} \lambda = 0  \tag{14.93}$$

移项： $$
2 \Delta \hat{\theta}_{n+1} = \tilde{X} \lambda  \tag{14.94}$$

于是：

$$
\Delta \hat{\theta}_{n+1} = \frac{1}{2} \tilde{X} \lambda  \tag{14.95}$$

这就是最优参数变化量的表达式——它必须是$\tilde{X}$各列的线性组合，即落在$\tilde{X}$的列空间中。这和NLMS单样本情形下$\Delta$必须沿着$X(n)$方向完全一致，只是现在空间变大了。

#### 6.3.3 步骤三：代入约束，求出 $\lambda$

将$\Delta \hat{\theta}_{n+1} = \frac{1}{2} \tilde{X} \lambda$代入约束$\tilde{e}(n) = \tilde{X}^T \Delta \hat{\theta}_{n+1}$： $$
\tilde{e}(n) = \tilde{X}^T \left( \frac{1}{2} \tilde{X} \lambda \right) = \frac{1}{2} \tilde{X}^T \tilde{X} \lambda  \tag{14.96}$$

两边乘以2：

$$
2 \tilde{e}(n) = \tilde{X}^T \tilde{X} \lambda  \tag{14.97}$$

如果$\tilde{X}^T \tilde{X}$可逆（即数据矩阵列满秩），则： $$
\lambda = 2 (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.98}$$

这正是你写的：

$$
\lambda = -(\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.99}$$

（符号差异取决于拉格朗日函数中约束项的符号约定。若用 $+\lambda^T(\tilde{e} - \tilde{X}^T \Delta)$，则 $\lambda$ 前面是正的 $2$。若用 $-\lambda^T(\tilde{e} - \tilde{X}^T \Delta)$，则出现负号。两者等价，只是 $\lambda$ 的定义差一个符号。）

#### 6.3.4 步骤四：求解参数更新表达式

将$\lambda$代回$\Delta \hat{\theta}_{n+1} = \frac{1}{2} \tilde{X} \lambda$： $$
\Delta \hat{\theta}_{n+1} = \frac{1}{2} \tilde{X} \cdot 2 (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n) = \tilde{X} (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.100}$$

因此：

$$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \tilde{X} (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.101}$$

这里你把$\tilde{X}$写成了$X(n)$，把$\tilde{e}(n)$写成了$e(n)$，形式完全一致： $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + X(n) (X(^T(n)X(n))^{-1} e(n)  \tag{14.102}$$

注意：在这个推广中，$X(n)$不再是单个列向量，而是$m \times (n-1)$的数据矩阵（每一列是一个历史输入向量），$e(n)$是$(n-1)$维误差向量。这个写法直接继承了NLMS的“形”，但内涵已经大不相同。

#### 6.3.5 步骤五：写出参数增量

从步骤4直接得到：

$$
\Delta \hat{\theta}_{n+1} = X(n) (X^T(n)X(n))^{-1} e(n)  \tag{14.103}$$

这就是参数的变化量。它的几何含义是：将误差向量$e(n)$通过伪逆$X(n)(X^T(n)X(n))^{-1}$映射到参数空间，得到最优修正方向。

#### 6.3.6 步骤六：误差与信息的关系

验证这个等式的链条：

第一项是误差的定义： $$
e(n) = \tilde{d}(n) - X^T(n) \hat{\theta}_{n}  \tag{14.104}$$

把约束$\tilde{d}(n) = X^T(n) \hat{\theta}_{n+1}$代入：

$$
e(n) = X^T(n) \hat{\theta}_{n+1} - X^T(n) \hat{\theta}_{n}  \tag{14.105}$$

提取公因式$X^T(n)$： $$
e(n) = X^T(n) (\hat{\theta}_{n+1} - \hat{\theta}_{n}) = X^T(n) \Delta \hat{\theta}_{n+1}  \tag{14.106}$$

这是一个漂亮的恒等式：**多约束情形下，误差向量$e(n)$就是参数变化量$\Delta \hat{\theta}_{n+1}$在数据矩阵$X(n)$的行空间上的投影。**

### 6.4 优化问题的最终形式

综合以上所有步骤，我们的优化问题呼之欲出：

$$
\min \|\Delta \hat{\theta}_{n+1}\|^2 \quad \text{s.t.} \ e(n) = \tilde{X}^T \Delta \hat{\theta}_{n+1}  \tag{14.107}$$

其中$\tilde{X} = (X_1, X_2, \ldots, X_{n-1})$是数据矩阵，$e(n) = \tilde{d} - \tilde{X}^T \hat{\theta}_n$是广义误差向量。

### 6.5 推广方法的命名

这个推广就是**仿射投影自适应滤波算法（Affine Projection Algorithm, APA）**。

- 当约束数量$L = 1$时（即只约束当前一个样本），APA退化为NLMS。
- 当约束数量$L = n$（即约束全部历史样本）时，APA退化为RLS（递归最小二乘）。

所以APA站在了NLMS和RLS的中间：

| 算法 | 约束数量 | 每步复杂度 | 收敛速度 | 数值稳定性 |
|---|---|---|---|---|
| **NLMS** | $L = 1$ | $O(m)$ | 慢 | 一般 |
| **APA** | $1 < L < n$ | $O(L^2 + Lm)$ | 中等 | 较好 |
| **RLS** | $L = n$ | $O(m^2)$ | 快 | 依赖于$P(n)$的条件数 |

仿射投影算法的名字来源于它的几何解释：每一步，我们把旧参数$\hat{\theta}_n$投影到由$L$个约束超平面张成的仿射子空间上（这些约束对应$L$个历史样本的完美拟合条件）。投影点就是新参数$\hat{\theta}_{n+1}$。

这就是NLMS的自然推广——从“盯住一个点”到“盯住一片区域”，视野开阔了，步子更稳了，计算量也上去了。后面我们会详细分析APA的收敛性能和工程实现。
### 6.6 欠定情形下的闭式解与正则化

现在我们来求解这个优化问题。需要指出的是，前面我们一直默认$\tilde{X}^T \tilde{X}$是可逆的，这在数据矩阵列满秩时成立。但在很多实际场景中（例如约束数量$L$大于参数维度$m$时），$\tilde{X}^T \tilde{X}$可能是奇异的。

所以我们要假设这个问题是**欠定**的，从伪逆的角度来重新审视。

### 6.6.1 欠定情形下的伪逆表达

由步骤4我们知道参数变化量的闭式解为： $$
\Delta \hat{\theta}_{n+1} = \tilde{X} (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.108}$$

这个表达式是在$\tilde{X}^T \tilde{X}$可逆的前提下得到的。当$\tilde{X}^T \tilde{X}$不可逆（即$\tilde{X}$列不满秩）时，这个公式无法直接使用。

但如果从伪逆的角度看，这个公式其实就是欠定情形下伪逆的写法。回顾我们在第二章欠定最小二乘中得到的结论：

$$
\Delta \hat{\theta} = \tilde{X}^T (\tilde{X} \tilde{X}^T)^{-1} \tilde{e}(n)  \tag{14.109}$$

这是超定情形（$L > m$，即约束个数多于参数个数）下的表达。而现在我们假设的是欠定情形（$L < m$，即约束个数少于参数个数），$\tilde{X} \tilde{X}^T$是$L \times L$的满秩矩阵，而$\tilde{X}^T \tilde{X}$是$m \times m$的奇异矩阵。

所以正确的伪逆表达应该是： $$
\Delta \hat{\theta}_{n+1} = \tilde{X}^T (\tilde{X} \tilde{X}^T)^{-1} \tilde{e}(n)  \tag{14.110}$$

但是，如果我们强行使用$\tilde{X} (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)$这个形式，当$\tilde{X}^T \tilde{X}$奇异时，$(\tilde{X}^T \tilde{X})^{-1}$不存在。因此我们需要引入**对角加载**（Diagonal Loading）来使其可逆。

### 6.6.2 对角加载：用正则化使矩阵可逆

引入对角加载的做法是在目标函数中加入一个正则项$\epsilon \|\theta\|^2$：

$$
L(\theta) = \|\tilde{d} (n) - \tilde{X}^T \theta \|^2 + \epsilon \|\theta\|^2  \tag{14.111}$$

这里$\epsilon > 0$是一个小的正则化参数。加了这一项之后，Hessian矩阵从$\tilde{X}^T \tilde{X}$变成了$\epsilon I + \tilde{X}^T \tilde{X}$。这个矩阵一定是满秩的（因为$\epsilon I$正定，$\tilde{X}^T \tilde{X}$半正定，两者之和正定）。

于是用牛顿法（实际上是正则化最小二乘的闭式解），更新公式变为： $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + (\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T \tilde{e}(n)  \tag{14.112}$$

注意：这里的"牛顿法"是一个广义的说法。严格来说，这是对带正则项的二次目标函数直接求梯度得到的闭式解，其形式与牛顿法一致——Hessian矩阵$\epsilon I + \tilde{X}^T \tilde{X}$取代了梯度下降中的步长，起到了对梯度方向进行“预条件”的作用。

### 6.6.3 利用Hessian的对称性化简

由于Hessian矩阵$H = \epsilon I + \tilde{X}^T \tilde{X}$是对称矩阵（因为$\tilde{X}^T \tilde{X}$对称，$\epsilon I$对称），所以$H^{-1}$也是对称的。

利用这个对称性，我们可以把$H^{-1} \tilde{X}^T$改写为$\tilde{X}^T H^{-1}$。验证如下：

$$
(\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T = \tilde{X}^T (\epsilon I + \tilde{X}^T \tilde{X})^{-1}  \tag{14.113}$$

这个等式的成立依赖于Hessian矩阵的对称性。只要Hessian对称，$H^{-1}$对称，那么$H^{-1} \tilde{X}^T = (\tilde{X} H^{-1})^T$，而$\tilde{X} H^{-1}$与$H^{-1} \tilde{X}^T$之间的关系需要仔细核对。

实际上我们来看两个表达式：

- $(\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T$是$m \times L$矩阵
- $\tilde{X}^T (\epsilon I + \tilde{X}^T \tilde{X})^{-1}$是$m \times L$矩阵

它们维度一致。但是我们用恒等式$(\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T = \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1}$更准确。这个恒等式可以通过Woodbury引理或直接验证得到。

不过，你写的等式 $$
(\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T = \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1}  \tag{14.114}$$

实际上才是正确的恒等式（它依赖于Hessian的对称性）。我们验证一下：两边左乘$(\epsilon I + \tilde{X}^T \tilde{X})$：

左边：$\tilde{X}^T$

右边：$(\epsilon I + \tilde{X}^T \tilde{X}) \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1} = \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)(\epsilon I + \tilde{X} \tilde{X}^T)^{-1} = \tilde{X}^T$

验证通过。所以：

$$
(\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T = \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1}  \tag{14.115}$$

因此更新公式可以写成： $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1} \tilde{e}(n)  \tag{14.116}$$

### 6.6.4 两种形式的对比

现在我们有了两个等价的正则化解：

| 形式 | 公式 | 适用场景 |
|---|---|---|
| **形式一** | $\hat{\theta}_{n+1} = \hat{\theta}_{n} + (\epsilon I + \tilde{X}^T \tilde{X})^{-1} \tilde{X}^T \tilde{e}(n)$ | 参数维度$m$不太大时 |
| **形式二** | $\hat{\theta}_{n+1} = \hat{\theta}_{n} + \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1} \tilde{e}(n)$ | 约束数量$L$不太大时 |

形式二的好处是：如果约束数量$L$远小于参数维度$m$（即$L \ll m$），那么$\epsilon I + \tilde{X} \tilde{X}^T$是一个$L \times L$的矩阵，求逆代价远低于对$m \times m$的矩阵求逆。

这就是在APA算法中常用的**数值优化技巧**：用低维矩阵的逆替代高维矩阵的逆，极大降低了计算复杂度。$\epsilon I + \tilde{X} \tilde{X}^T$中的对角加载$\epsilon I$保证了矩阵的正定性，同时$\epsilon$的选择也控制了正则化强度——$\epsilon$越大，步长越小，收敛越慢但越稳定；$\epsilon$越小，收敛越快但对奇异越敏感。

## 7. 变种：符号误差算法

### 7.1 用符号函数代替误差

在标准的NLMS中，更新量的方向由误差向量$\tilde{e}(n)$决定，大小由$\mu$控制：

$$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \tilde{X}^T \tilde{e}(n)  \tag{14.117}$$

一个很自然的变种是：**只保留误差的符号信息，丢掉幅度信息**。 $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \tilde{X}^T \operatorname{sgn}(\tilde{e}(n))  \tag{14.118}$$

其中$\operatorname{sgn}(\cdot)$是逐分量作用的符号函数：

$$
\operatorname{sgn}(x) = \begin{cases} 1, & x > 0 \\ 0, & x = 0 \\ -1, & x < 0 \end{cases}  \tag{14.119}$$

**直观解释**：只要知道误差的方向就够了——是正还是负，决定了参数应该往哪个方向调。至于误差具体有多大，全部被步长$\mu$吸收掉了。这相当于把误差信息从“连续值”压缩成了“三态值”（$+1, 0, -1$）。

**工程好处**：符号函数在硬件上极其便宜——不需要乘法，只需要比较器和选择器。在低功耗、低成本DSP或FPGA实现中，省掉一个乘法器就是省掉一大片硅片。

**理论代价**：收敛速度变慢，稳态误差变大。这是用性能换硬件的典型交易。

---

### 7.2 理论基础：最小一乘（Least Absolute Deviations）

为什么用符号函数是有道理的？因为它对应着一个不同的优化准则——**最小一乘**。

#### 7.2.1 从最小二乘到最小一乘

标准LMS/NLMS最小化的是**平方误差**： $$
J_2(\theta) = \|\tilde{d} - \tilde{X}^T \theta\|^2 = \sum_{i} (e_i)^2  \tag{14.120}$$

而符号误差算法对应的是最小化**绝对误差**：

$$
J_1(\theta) = \|\tilde{d} - \tilde{X}^T \theta\|_1 = \sum_{i} |e_i|  \tag{14.121}$$

这个差别看起来只是把平方换成了绝对值，但后果极其深远。

#### 7.2.2 为什么最小一乘是稳健的（Robust）

最小二乘对 **异常值（outlier）** 极其敏感——一个离谱的误差点，平方之后会主导整个损失函数，把估计结果拉偏。

最小一乘只惩罚误差的绝对值，不是平方。一个异常值的影响是线性的，而不是平方级的。所以最小一乘对异常值不那么敏感，是一种**稳健估计**方法。

**一句话**：最小二乘容易被极端值带偏，最小一乘扛得住。

---

### 7.3 概率角度：中位数 vs 均值

从概率统计的角度看，最小一乘和最小二乘对应着完全不同的“中心”概念。

#### 7.3.1 均值是最小二乘的解

对于一组数据$\{X_i\}$，最小二乘估计$\hat{a}$使平方误差最小： $$
\min_a \mathbb{E}[(X - a)^2] \quad \Longrightarrow \quad a = \mathbb{E}[X] = \mu  \tag{14.122}$$

所以**均值（mean）是最小二乘意义上的最优中心**。

样本均值：

$$
\bar{X} = \frac{1}{n} \sum_{i=1}^n X_i  \tag{14.123}$$

#### 7.3.2 中位数是最小一乘的解

对于一组数据$\{X_i\}$，最小一乘估计$\hat{a}$使绝对误差最小： $$
\min_a \mathbb{E}[|X - a|] \quad \Longrightarrow \quad a = \text{median}(X)  \tag{14.124}$$

**中位数（median）是最小一乘意义上的最优中心**。

样本中位数：对$X$排序后取中间值

$$
\text{median}(X) = X_{(n/2)}  \tag{14.125}$$

#### 7.3.3 对比

| 准则 | 最优中心 | 对异常值敏感度 | 稳健性 |
|---|---|---|---|
| 最小二乘（$L_2$） | 均值 | 高 | 差 |
| 最小一乘（$L_1$） | 中位数 | 低 | 好 |

这就是为什么符号误差算法更稳健——它本质上是在做$L_1$范数优化，而$L_1$范数的最优点是中位数，中位数对极端值不敏感。

---

### 7.4 次梯度：处理不可导的$L_1$范数

$J_1(\theta) = \sum_i |e_i|$在$e_i = 0$处**不可导**。用标准梯度下降法没法直接处理。

这时候需要**次梯度（Subgradient）**。

#### 7.4.1 次梯度的定义

对于凸函数$f: \mathbb{R}^n \to \mathbb{R}$，在点$x$处的次梯度$g$满足： $$
f(y) \geq f(x) + g^T (y - x), \quad \forall y  \tag{14.126}$$

所有满足这个条件的$g$组成的集合称为**次微分**$\partial f(x)$。

#### 7.4.2 绝对值函数的次梯度

对于$f(x) = |x|$：

- 当$x > 0$时，$\partial f(x) = \{1\}$（就是普通导数）
- 当$x < 0$时，$\partial f(x) = \{-1\}$（就是普通导数）
- 当$x = 0$时，$\partial f(0) = [-1, 1]$（整个区间都是次梯度）

也就是说，在尖点处，**从$-1$到$1$之间的任何一个数都可以作为次梯度**。

#### 7.4.3 最小一乘的次梯度

对于$J_1(\theta) = \|\tilde{d} - \tilde{X}^T \theta\|_1 = \sum_i |e_i|$，其中$e_i = d_i - X_i^T \theta$：

对$\theta$求次梯度（链式法则）：

$$
\partial J_1(\theta) = - \sum_i X_i \cdot \partial |e_i|  \tag{14.127}$$

其中： $$
\partial |e_i| = \begin{cases}
\{1\}, & e_i > 0 \\
\{-1\}, & e_i < 0 \\
[-1, 1], & e_i = 0
\end{cases}  \tag{14.128}$$

如果我们**在$e_i = 0$时随便选一个次梯度**（比如选0），并且忽略掉那些$e_i = 0$的项，那么：

$$
\partial J_1(\theta) \approx - \sum_i X_i \cdot \operatorname{sgn}(e_i)  \tag{14.129}$$

**次梯度就是符号函数**。

#### 7.4.4 用次梯度做更新

次梯度下降的更新公式是： $$
\theta_{n+1} = \theta_n - \mu \cdot g, \quad g \in \partial J_1(\theta_n)  \tag{14.130}$$

把$g \approx - \sum_i X_i \operatorname{sgn}(e_i)$代进去：

$$
\theta_{n+1} = \theta_n + \mu \sum_i X_i \operatorname{sgn}(e_i)  \tag{14.131}$$

写成矩阵形式就是： $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \tilde{X}^T \operatorname{sgn}(\tilde{e}(n))  \tag{14.132}$$

**这就是符号误差算法的次梯度解释。**

---

### 7.5 符号函数的三个变种

文献中通常把基于符号的自适应算法分为三类：

| 变种 | 更新公式 | 对什么取符号 |
|---|---|---|
| **符号误差（Signed Error）** | $\theta_{n+1} = \theta_n + \mu X(n) \operatorname{sgn}(e(n))$ | 只对误差取符号 |
| **符号回归量（Signed Regressor）** | $\theta_{n+1} = \theta_n + \mu \operatorname{sgn}(X(n)) e(n)$ | 只对输入取符号 |
| **符号-符号（Sign-Sign）** | $\theta_{n+1} = \theta_n + \mu \operatorname{sgn}(X(n)) \operatorname{sgn}(e(n))$ | 对误差和输入都取符号 |

其中**符号-符号算法**计算量最小——连乘法都没有了，只剩下符号位的异或操作。但收敛也最慢。

你写的$\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \tilde{X}^T \operatorname{sgn}(e(n))$是**符号误差**变种在APA框架下的推广。

---

### 7.6 CCITT标准中的符号算法

符号误差类算法在通信标准中有实际应用。CCITT（现ITU-T）在32 kbps ADPCM（自适应差分脉冲编码调制）标准中采用了基于符号的自适应预测/滤波算法。

#### 7.6.1 为什么CCITT用符号算法？

电话语音信号的特点是：
- 语音信号的动态范围极大（从耳语到喊叫能差40-60dB）
- 但短期内语音信号的统计特性相对平稳

符号算法在这种场景下的优势：
- **对信号幅度不敏感**：符号函数把幅度信息全部丢掉，只保留方向。在大动态范围的语音信号中，这反而成了一种“自动增益控制”——不会因为信号忽大忽小导致步长失效。
- **硬件便宜**：电话交换机里有成百上千条线路同时做回声消除，每个通道都要跑自适应滤波。符号算法省下来的乘法器，在硬件上就是实打实的成本。
- **符合G.165/G.168标准**：ITU-T G.165和G.168是回声消除器的国际标准，要求算法在特定条件下达到一定的收敛速度和稳态性能。符号类算法虽然收敛慢，但在电话回声消除的场景下（回声路径相对稳定，不需要极速收敛）完全够用。

#### 7.6.2 具体是哪一种符号变种？

CCITT 32 kbps ADPCM标准中使用的自适应滤波器，实际上是一种**符号-符号算法的变种**——对误差取符号，对输入信号也做某种形式的“剪裁”（clipping）。

具体来说，ADPCM中的自适应量化器步长更新用的是：

$$
\Delta(n+1) = \Delta(n) \cdot M(n)  \tag{14.133}$$

其中$M(n)$由量化误差的符号和 magnitude 共同决定。这种更新策略本质上就是在用符号信息来控制步长的调整方向，和符号误差算法的精神一脉相承。

#### 7.6.3 为什么叫“伪对数”（Pseudo-Logarithmic）？

ADPCM标准中还有一个值得注意的细节：信号采用的是**A-law或μ-law伪对数编码**。

这种编码的本质是：**对小信号精细量化，对大信号粗糙量化**——在数字域用8位就实现了相当于12-13位线性量化的动态范围。

符号算法和伪对数编码**天然匹配**：
- 伪对数编码已经把信号的幅度信息压缩了
- 符号算法彻底丢掉幅度信息，只留方向
- 两者叠加，相当于在信号进入滤波器之前做了一次“粗粒度”处理

用符号误差替代真实误差，实际上就是在**牺牲精度换简化**——这在电话质量的语音通信中是合理的，因为人耳对相位（方向）比幅度更敏感。

---

### 7.7 本节总结

符号误差算法不是瞎搞，它有三层理论支撑：

| 层次 | 内容 |
|---|---|
| **优化准则** | 最小一乘（$L_1$范数）替代最小二乘（$L_2$范数） |
| **统计含义** | 中位数替代均值，对异常值更稳健 |
| **数学工具** | 次梯度替代梯度，符号函数是$\|x\|$的次梯度 |

工程上，符号算法用**计算量的降低**换取了**收敛速度的牺牲**，在语音通信这类对成本敏感、对精度容忍度高的场景中找到了自己的位置。CCITT ADPCM标准就是它的“身份证”——不是实验室里的玩具，是跑在现实电话网络里的东西。

## 8. 课后总结

本章从LMS的尺度敏感性问题出发，逐步建立了NLMS的完整理论框架。现在我们以公式为核心，系统回顾每一节的关键结论。

---

### 8.1 LMS的尺度敏感性

标准LMS更新公式： $$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu X(n) e(n), \qquad e(n) = d(n) - X^T(n) \hat{\theta}_n  \tag{14.134}$$

当$X(n)$和$d(n)$同时扩大2倍时，更新量$X(n)e(n)$扩大4倍：

$$
X'(n) e'(n) = (2X(n)) \cdot (2e(n)) = 4 X(n) e(n)  \tag{14.135}$$

**核心结论**：同一问题，最优解$\theta$不变，但有效步长被放大了4倍。

---

### 8.2 Goodwin约束优化（角度一）

**优化问题**： $$
\min_{\hat{\theta}_{n+1}} \|\hat{\theta}_{n+1} - \hat{\theta}_{n}\|^2, \quad \text{s.t.} \quad d(n) = X^T(n) \hat{\theta}_{n+1}  \tag{14.136}$$

**直观比喻**：让高年级的做低年级的题，要全部能做对的前提下，年级差异最小。

**拉格朗日乘子法求解**：

令$\Delta = \hat{\theta}_{n+1} - \hat{\theta}_n$，约束变为$e(n) = X^T(n)\Delta$。

拉格朗日函数：

$$
L(\Delta, \lambda) = \|\Delta\|^2 + \lambda(e(n) - X^T(n)\Delta)  \tag{14.137}$$

求梯度并令为零： $$
\nabla_{\Delta}L = 2\Delta - \lambda X(n) = 0 \quad \Longrightarrow \quad \Delta = \frac{\lambda}{2}X(n)  \tag{14.138}$$

代入约束：

$$
e(n) = X^T(n) \cdot \frac{\lambda}{2}X(n) = \frac{\lambda}{2}\|X(n)\|^2
\quad \Longrightarrow \quad \frac{\lambda}{2} = \frac{e(n)}{\|X(n)\|^2}  \tag{14.139}$$

因此： $$
\Delta = \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.140}$$

代回$\Delta = \hat{\theta}_{n+1} - \hat{\theta}_n$：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.141}$$

**结论**：NLMS是带信任域约束的瞬时误差最小化问题的精确闭式解。

---

### 8.3 最大化差距（角度二）

**核心思想**：高年级做低年级的题，得分差距要尽可能拉大。

定义低年级误差$e_1(n)$和高年级误差$e_2(n)$： $$
e_1(n) = d(n) - X^T(n) \hat{\theta}_n, \qquad e_2(n) = d(n) - X^T(n) \hat{\theta}_{n+1}  \tag{14.142}$$

最大化差距：

$$
\max |e_2(n) - e_1(n)|, \quad \text{s.t.} \quad |e_2(n)| < |e_1(n)|  \tag{14.143}$$

展开差距： $$
e_2(n) - e_1(n) = -X^T(n) \Delta \hat{\theta}_{n+1}  \tag{14.144}$$

由柯西不等式：

$$
|X^T(n) \Delta \hat{\theta}_{n+1}| \leq \|X(n)\| \cdot \|\Delta \hat{\theta}_{n+1}\|  \tag{14.145}$$

等号成立条件：$\Delta \hat{\theta}_{n+1}$与$X(n)$平行。令$\Delta \hat{\theta}_{n+1} = k X(n)$。

结合约束$|e_2(n)| < |e_1(n)|$，令$e_2(n) = (1-\mu)e_1(n)$，$\mu \in (0,1)$，得到： $$
k = \mu \frac{e_1(n)}{\|X(n)\|^2}  \tag{14.146}$$

因此：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.147}$$

**结论**：两个看似对立的视角——"最小变化"和"最大差距"——得到相同结果，指向同一个公式。

---

### 8.4 正交投影（角度三）

定义残差： $$
\epsilon(n) = \hat{\theta}_n - \hat{\theta}_{\text{opt}}  \tag{14.148}$$

误差与残差的关系：

$$
e(n) = -X^T(n) \epsilon(n)  \tag{14.149}$$

逐次正交化（Gram-Schmidt精神）： $$
\epsilon(n+1) = \epsilon(n) - \text{Proj}_{X(n)} \epsilon(n)  \tag{14.150}$$

一维投影算子：

$$
\text{Proj}_{X(n)} \epsilon(n) = \frac{X(n) X^T(n)}{\|X(n)\|^2} \epsilon(n)  \tag{14.151}$$

代入递推式，换回$\hat{\theta}$表示： $$
\hat{\theta}_{n+1} = \hat{\theta}_n + \frac{X(n)}{\|X(n)\|^2} e(n)  \tag{14.152}$$

**结论**：NLMS每一步都在从系数残差中减去它在当前输入方向上的投影，是一个逐次正交化过程。

---

### 8.5 牛顿法视角（角度四）

单样本瞬时最小二乘目标函数：

$$
f(\theta) = \|d(n) - X^T(n)\theta\|^2  \tag{14.153}$$

梯度： $$
\nabla_{\theta} f(\theta) = -2X(n)(d(n) - X^T(n)\theta)  \tag{14.154}$$

Hessian矩阵：

$$
H(f(\theta)) = 2X(n)X^T(n)  \tag{14.155}$$

当$m > 1$时，$X(n)X^T(n)$秩最多为1，Hessian**缺秩**，牛顿法无法直接使用。

引入正则项$\epsilon\|\theta\|^2$后的Hessian： $$
H(f(\theta)) = 2(\epsilon I + X(n)X^T(n))  \tag{14.156}$$

带正则项的牛顿法更新公式：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \left[ \frac{X(n)}{\epsilon + \|X(n)\|^2} e(n) - (\epsilon I + X(n)X^T(n))^{-1} \epsilon \hat{\theta}_n \right]  \tag{14.157}$$

**结论**：直接套牛顿法行不通（Hessian缺秩）；加正则项强行可逆后多出一个收缩项；扔掉收缩项才回到NLMS。这条路径的价值是告诉我们为什么直接套牛顿法不可行。

---

### 8.6 仿射投影算法（APA）：NLMS的多约束推广

**单约束**（NLMS）： $$
d(n) = X^T(n) \hat{\theta}_{n+1}  \tag{14.158}$$

**多约束**（APA）：

$$
d(k) = X^T(k) \hat{\theta}_{n+1}, \quad k = 1, 2, \ldots, L  \tag{14.159}$$

写成矩阵形式： $$
\tilde{d} = \tilde{X}^T \hat{\theta}_{n+1}  \tag{14.160}$$

优化问题：

$$
\min_{\hat{\theta}_{n+1}} \|\hat{\theta}_{n+1} - \hat{\theta}_n\|^2, \quad \text{s.t.} \quad \tilde{d} = \tilde{X}^T \hat{\theta}_{n+1}  \tag{14.161}$$

定义广义误差： $$
\tilde{e}(n) = \tilde{d} - \tilde{X}^T \hat{\theta}_n  \tag{14.162}$$

拉格朗日乘子法解得：

$$
\hat{\theta}_{n+1} = \hat{\theta}_n + \tilde{X} (\tilde{X}^T \tilde{X})^{-1} \tilde{e}(n)  \tag{14.163}$$

或正则化形式： $$
\hat{\theta}_{n+1} = \hat{\theta}_n + \tilde{X}^T (\epsilon I + \tilde{X} \tilde{X}^T)^{-1} \tilde{e}(n)  \tag{14.164}$$

**三种算法的统一视角**：

| 算法 | 约束数量$L$ | 每步复杂度 |
|---|---|---|
| NLMS | $L = 1$ | $O(m)$ |
| APA | $1 < L < m$ | $O(L^2 + Lm)$ |
| RLS | $L = n$（全部历史） | $O(m^2)$ |

---

### 8.7 符号误差变种

用符号函数代替误差：

$$
\hat{\theta}_{n+1} = \hat{\theta}_{n} + \mu \tilde{X}^T \operatorname{sgn}(\tilde{e}(n))  \tag{14.165}$$

**理论基础**：最小一乘（$L_1$范数）替代最小二乘（$L_2$范数） $$
J_1(\theta) = \|\tilde{d} - \tilde{X}^T \theta\|_1 = \sum_i |e_i|  \tag{14.166}$$

**概率含义**：中位数替代均值，对异常值更稳健。

**次梯度**：$|x|$在$x=0$处的次梯度是$[-1, 1]$，符号函数是其中的一个代表。

**三个变种**：

| 变种 | 更新公式 |
|---|---|
| 符号误差 | $\theta_{n+1} = \theta_n + \mu X(n) \operatorname{sgn}(e(n))$ |
| 符号回归量 | $\theta_{n+1} = \theta_n + \mu \operatorname{sgn}(X(n)) e(n)$ |
| 符号-符号 | $\theta_{n+1} = \theta_n + \mu \operatorname{sgn}(X(n)) \operatorname{sgn}(e(n))$ |

**工程意义**：CCITT ADPCM标准中采用了符号类算法，用计算量的降低换取硬件的简化，在电话语音通信中被证明是可行的。

---

### 8.9 本章公式汇总

| 公式 | 编号 | 说明 |
|---|---|---|
| $\hat{\theta}_{n+1} = \hat{\theta}_n + \mu X(n) e(n)$ | (1) | 标准LMS |
| $\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \frac{X(n)}{\|X(n)\|^2} e(n)$ | (2) | NLMS（工程形式） |
| $\min \|\Delta\|^2, \text{s.t. } e(n) = X^T(n)\Delta$ | (3) | Goodwin约束优化 |
| $\Delta = \frac{X(n)}{\|X(n)\|^2} e(n)$ | (4) | 角度一的闭式解 |
| $\epsilon(n+1) = \epsilon(n) - \text{Proj}_{X(n)}\epsilon(n)$ | (5) | 正交投影递推 |
| $\hat{\theta}_{n+1} = \hat{\theta}_n + \tilde{X}(\tilde{X}^T\tilde{X})^{-1}\tilde{e}(n)$ | (6) | APA闭式解 |
| $\hat{\theta}_{n+1} = \hat{\theta}_n + \tilde{X}^T(\epsilon I + \tilde{X}\tilde{X}^T)^{-1}\tilde{e}(n)$ | (7) | APA正则化形式 |
| $\hat{\theta}_{n+1} = \hat{\theta}_n + \mu \tilde{X}^T \operatorname{sgn}(\tilde{e}(n))$ | (8) | 符号误差APA |

---

### 👑 五大算法 vs. QR-RLS 对比表

| 算法 | 核心代价函数 | 求解/迭代公式 | 核心参数 | 计算复杂度 | 一句话灵魂（本质） |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **最小二乘 (LS)** (批量) | \(\sum_{i=1}^{N} e_i^2\) | \(\hat{\theta} = (X^T X)^{-1} X^T y\) | 无 | \(O(N^3)\) | **“翻旧账”**：数据攒齐算总账。 |
| **总体最小二乘 (TLS)** (批量) | \(\min \|[\Delta X, \Delta y]\|_F^2\) | 对增广矩阵做 **SVD分解** | 无 | \(O(N^3)\) | **“公平怪”**：输入输出都可能有错，找最靠谱的垂直距离。 |
| **LMS** | \(E[e^2]\) (理论) <br> **实际用**：\(e^2(n)\) | \(\omega(n+1) = \omega(n) + \mu e(n) x(n)\) | **步长 \(\mu\)** | \(O(N)\) | **“近视眼”**：只看脚下这一步，凭感觉瞎摸。 |
| **NLMS** | 同LMS | \(\omega(n+1) = \omega(n) + \frac{\mu}{\|x(n)\|^2+\epsilon} e(n) x(n)\) | **步长 \(\mu\)** | \(O(N)\) | **“戴眼镜的近视眼”**：怕输入太强走过头，主动缩步子。 |
| **标准RLS** | \(\sum_{i=0}^{n} \lambda^{n-i} e^2(i)\) | **\(k(n) = \frac{P(n-1)x(n)}{\lambda + x^T(n)P(n-1)x(n)}\)** <br> **\(\omega(n) = \omega(n-1) + k(n)e(n)\)** | **遗忘因子 \(\lambda\)** | \(O(N^2)\) | **“神算子”**：记住所有历史，算出最优方向。 |
| **QR-RLS** | 同标准RLS | **维护矩阵 \(R\)**（\(X^TX\)的平方根），用**Givens旋转**等正交变换进行**递归更新**<br>最终解：\(R(n)\omega(n) = u(n)\) | **遗忘因子 \(\lambda\)** | \(O(N^2)\) | **“带了计算器的神算子”**：每一步都用更稳健的数学工具（QR分解）精确计算。 |


### 8.10 本章整体脉络

$$
\boxed{\text{LMS尺度敏感} \longrightarrow \text{NLMS工程修补} \longrightarrow \text{Goodwin约束优化} \longrightarrow \text{四角度统一证明}}  \tag{14.167}$$ 

$$
\downarrow  \tag{14.168}$$

$$
\boxed{\text{APA多约束推广} \longrightarrow \text{正则化与伪逆} \longrightarrow \text{符号误差变种与$L_1$范数}}  \tag{14.169}$$

**核心结论**：

1. **LMS的尺度敏感性**是结构性的，根源在于更新量$X(n)e(n)$对输入输出幅度的二次依赖。

2. **NLMS**不是一个拍脑袋的归一化——它同时是约束优化、极值问题和正交投影的必然结果。

3. **四个角度**（约束优化、最大化差距、正交投影、正则化牛顿法）从不同路径汇聚到同一个公式，说明NLMS的更新规则具有数学上的多重必然性。

4. **APA**是NLMS的自然推广：从"盯住一个样本"扩展到"盯住一片历史窗口"，约束数量的增加换来了更快的收敛速度。

5. **符号误差变种**用$L_1$范数替代$L_2$范数，用中位数替代均值，用次梯度替代梯度——以牺牲精度为代价，换来了硬件简化和对异常值的稳健性。

6. 从NLMS到APA再到符号算法，整个发展路径展示了**同一套正交投影思想在不同约束条件和不同硬件约束下的适应性演化**。

---

### 8.10 学习检查清单

- [ ] 能写出 NLMS 的更新公式：$\mathbf{w}(n+1) = \mathbf{w}(n) + \mu \frac{\mathbf{x}(n)}{\|\mathbf{x}(n)\|^2} e(n)$，并解释分母的作用
- [ ] 能说明 LMS 尺度敏感性的根源：更新量 $\mathbf{x}(n)e(n)$ 对输入输出幅度有二次依赖
- [ ] 能从四个等价角度理解 NLMS：（1）约束最小范数更新（2）最大化差距（3）正交投影（4）正则化牛顿法
- [ ] 能写出 Goodwin 约束优化问题：$\min\|\Delta\mathbf{w}\|^2$ s.t. $e(n) = \mathbf{x}^\top(n)\Delta\mathbf{w}$，并推导闭式解
- [ ] 能解释正交投影视角：NLMS 的权值更新 = 权值误差在当前输入方向上的正交投影
- [ ] 能写出 APA（仿射投影算法）的更新公式并说明它与 NLMS 的关系：NLMS 是 APA 在投影阶数 $P=1$ 时的特例
- [ ] 能说明为什么 APA 比 NLMS 收敛更快：同时利用多个历史时刻的信息做多维投影
- [ ] 能区分三种符号算法（符号误差、符号回归量、符号-符号）及其对鲁棒性和计算复杂度的权衡

### 8.11 思考题

1. **NLMS 的四视角统一的哲学意义**：同一个公式可以从约束优化、最大化差距、正交投影和正则化牛顿法四个角度推导出来。这是巧合还是必然？"多重推导的统一"在数学和物理中往往意味着触及了深层本质——NLMS 的深层本质是什么？

2. **NLMS 与 APA 的投影几何**：NLMS 每次投影到当前输入 $\mathbf{x}(n)$ 的正交补上，APA 投影到最近 $P$ 个输入张成的子空间的正交补上。为什么多维投影能加速收敛？投影维度 $P$ 的增大有什么代价——除了计算量增加，是否还有统计上的风险？

3. **正则化 APA 与对角加载**：APA 的正则化形式 $(\epsilon I + \tilde{X}\tilde{X}^\top)^{-1}$ 中的 $\epsilon$ 与岭回归中的 $\lambda$ 作用相同。在自适应滤波中，$\epsilon$ 不仅改善数值稳定性，还能在初始阶段防止系数发散——这个双重角色的物理本质是什么？

4. **符号算法的"L1 哲学"**：符号误差算法本质上是将 $L_2$ 损失换为 $L_1$ 损失（中位数替代均值）。在存在脉冲噪声（outliers）的环境中，为什么中位数比均值更鲁棒？这种鲁棒性是以什么代价换来的（提示：统计效率）？

5. **从 NLMS 到深度学习优化器**：NLMS 的输入归一化思想和 Adam 等自适应学习率方法有相似之处——都是根据数据的二阶统计量来调整更新幅度。两者在设计哲学上的根本区别在哪里？为什么深度学习选择了 Adam 而非 NLMS 式的显式归一化？


把 **QR-RLS** 加入战局后，这张“线性估计算法对比表”就真的完整了。它和标准RLS是“亲兄弟”，但为了解决标准RLS的一个致命弱点——**数值稳定性**——而诞生。

QR-RLS的核心思想是，不再直接处理协方差矩阵 `P`，而是对输入数据矩阵进行**QR分解**（正交三角分解）。它也被称为**平方根RLS（Square Root RLS）**。

<div style="page-break-before: always;"></div>