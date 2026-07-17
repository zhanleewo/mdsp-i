<div style="display: flex; flex-direction: column; justify-content: center; align-items: center; height: 95vh; text-align: center; padding: 0 5%; background: 
    radial-gradient(ellipse at 20% 30%, rgba(56, 189, 248, 0.20) 0%, transparent 55%),
    radial-gradient(ellipse at 80% 70%, rgba(250, 204, 21, 0.10) 0%, transparent 45%),
    linear-gradient(145deg, #0a0f2a 0%, #1a1a4a 45%, #0d1b2a 85%, #050814 100%);">

  <h1 style="font-size: 3.2rem; margin-bottom: 0.3rem; margin-top: 120px; color: white; border-bottom: none;">现代数字信号处理 I 笔记</h1>
  <hr style="width: 120px; height: 2px; border: none; margin: 1.8rem auto; background: linear-gradient(to right, transparent, #38bdf8, #facc15, #38bdf8, transparent); width: 62%">
  <p style="font-size: 2.2rem; margin-top: 1rem; color: white;">我给 AI 画个谱</p>
  
  <div style="flex:1; width:90%; max-width:900px; min-height:200px; margin:8px 0; z-index:2;">
    <img src="assets/00/page-0.svg" alt="现代数字信号处理 I 讲义" style="width:100%; height:100%; display:block; background:transparent;">
  </div>
  <p style="font-size: 1.3rem; margin-top: auto; color: white;">https://github.com/zhanleewo/mdsp-i</p>

</div>

<div style="page-break-before: always;"></div>
<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="引言" style="text-align: center; margin-bottom: 2rem; border-bottom: none;">引言</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

## 0.1 先修知识

阅读本讲义需要具备以下数学基础：
- **微积分**：求导、积分、级数展开，是理解概率密度函数和优化问题的基本工具。
- **线性代数**：矩阵运算、特征值分解、二次型，是处理多维随机变量和信号模型的必备语言。
- **概率论**：至少达到“小脑反应”的熟练度——即看到随机变量能本能想到分布、期望和方差，看到两个变量能本能想到联合、边缘和条件。
- **复变函数**：了解留数定理和傅里叶变换在复数域上的推广，对谱分析有帮助（少量涉及）。
- **信号与系统 / 数字信号处理**：非强制，但若有背景，会对滤波、采样、频域分析等章节有更直观的把握。

## 0.2 参考书目

本讲义并非凭空杜撰，而是以下经典著作核心思想的提炼与串联：
1. **《统计信号处理基础》**（Kay）：估计与检测理论的圣经，本书的核心算法框架来源于此。
2. **《自适应滤波器》**（Haykin）：将统计信号处理推向工程实现的经典，涉及维纳滤波、LMS、RLS。
3. **《现代信号谱分析》**（Stoica & Moses）：从周期图到现代参数化谱估计的完整图谱。
4. **《概率机器学习：一种可引入的方法》**（Murphy）：作为贝叶斯视角的补充读物，可选。

## 0.3 三条格言：贯穿全书的哲学底色

> **No data, no truth.**
>
> **No analytics, no understanding.**
>
> **No programming, no cognition.**

这三句话，不仅是口号，更是本讲义循序渐进的三个层次。

**第一句：No data, no truth（无数据，无真相）**

在信号处理的世界里，真相（信号的幅度、频率、到达时间、是否存在目标）隐藏在物理世界的某个角落。我们永远无法像上帝一样直接窥见真相，我们唯一能抓住的，是传感器采样回来的一串有限长的数字：$x[0], x[1], \dots, x[N-1]$。数据是我们与客观世界之间唯一的、不可逾越的桥梁。脱离数据的“推理”只是空中楼阁。这对应着**统计三要素**中的**“数据”**。

**第二句：No analytics, no understanding（无分析，无理解）**

原始数据只是一堆冰冷的数字。如果没有数学工具去剖析它——如果没有用概率分布去描述噪声的统计规律，如果没有用矩阵分解去提取信号子空间，如果没有用假设检验去量化决策的置信度——那么数据就永远只是数据，无法转化为洞见。分析，是我们赋予数据意义的唯一途径。这对应着**统计三要素**中的**“模型”**与**“决策”**中的数学推导。

**第三句：No programming, no cognition（无编程，无认知）**

理论分析给出了优雅的公式（如维纳-霍夫方程、卡尔曼递推），但真实的信号不会自动代入公式求解。我们需要通过编程（如 MATLAB / Python）将算法实例化，在仿真中验证性能，在实测数据中检验鲁棒性。**编程不是把公式敲进电脑，而是强迫我们将模糊的数学直觉转化为精准的数值逻辑——这个过程本身就是一种认知的升华。** 当算法在屏幕上跑出结果的那一刻，我们才真正“理解”了它。本讲义会给出关键的算法伪代码，但真正的实现留给你亲手完成。

## 0.4 全书的核心主线：数据、模型、决策

将上面的哲学落到地面上，本书所有内容可以用三个词概括：

1. **数据**（$\mathbf{x}$）：观测向量，是我们拥有的一切。
2. **模型**（$p(\mathbf{x};\theta)$ 或 $p(\mathbf{x}|\theta)$）：我们对“数据是如何产生的”这一物理过程的数学假设。
3. **决策**（$\hat{\theta}=g(\mathbf{x})$ 或假设选择）：我们最终要采取的行动，即参数估计或信号检测。

> **数据是上帝给的（不可更改），模型是人造的（可商榷）。**

这句话极其重要。**数据是客观事实**，我们无法挑选、无法篡改（只能尽量保证测量准确）。**模型是主观假设**，它取决于工程师对物理世界的理解——你认为噪声是高斯还是拉普拉斯？你认为信号是确定性正弦波还是随机过程？不同的模型选择，会导出完全不同的估计算法和检测性能。优秀的信号处理工程师，不在于他推导公式有多快，而在于他能否为一个实际问题**构建出既逼近物理真相、又在数学上可处理的模型**。


**约定：不特别说明的情况下，不加转置表示列向量，加转置表示行向量：$a$ 表示列向量，$a^\top$ 表示行向量。**

**课程录像不包含第二讲的内容，这是我自己添加的。**

**课程视频:** https://www.bilibili.com/video/BV1ga4y157L5

| 维度 | 现代数字信号处理 I | 现代数字信号处理 II |
|------|-------------------|-------------------|
| 纪元 | 1950 – 1980 | 1980 – now |
| 系统模型 | 线性 | 非线性 |
| 过程假设 | 宽平稳 | 非平稳 |
| 噪声/分布假设 | 高斯 | 非高斯 |
| 分解理论 | 正交 | 非正交 |