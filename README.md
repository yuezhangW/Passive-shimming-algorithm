# Passive-shimming-algorithm
The module for passive shimming algorithm

本项目基于Python语言进行编写，独立于无源匀场辅助软件，针对已有传感器采集数据进行无源匀场实验，对比不同算法下的匀场效果差异。项目提出了一种2-stage方案，对整数规划进行处理。

This project is written in Python, independent of passive shimming software, to conduct passive leveling experiments for existing sensor acquisition data and compare the difference of shimming effect under different algorithms. The project proposes a 2-stage scheme to process the integer planning. (Translated from DeepL)



## 理论推导

### 数据矩阵

##### 裸磁场数据矩阵

程序中用类Dataloader保存外部裸磁场数据，通过函数getMagneticFieldMatrix返回磁场强度数据矩阵，文件"6009map1.txt"得到的矩阵size为(33,32)，文件"9006map1.txt"得到的矩阵size为(33,24)（phi方向，theta方向）。



##### 核磁共振仪器匀场片分布

NMR的匀场片分布尺寸可由用户指定，以(40,24)为例，指phi方向上有40片，z方向上有24片，共有960片。



##### 裸磁场数据一维向量

在计算灵敏度系数矩阵和匀场算法前，需要将二维的裸磁场数据reshape到一维。行表示phi方向，列表示theta方向。
$$
\left[
	\begin{matrix}
		A_{11} & A_{12} & A_{13} & A_{14} \\
		A_{21} & A_{22} & A_{23} & A_{24} \\
		A_{31} & A_{32} & A_{33} & A_{34} \\
		
		A_{41} & A_{42} & A_{43} & A_{44} 
	\end{matrix} 
\right]
\rightarrow
\left[
	\begin{matrix}
		A_{11}\\
		A_{12}\\
		A_{13}\\
		A_{14}\\
		...\\
		A_{42}\\
		A_{43}\\
		A_{44}\\
	\end{matrix} 
\right]
$$
(33,32)->(1056,)  (33,24)->(796,)



##### 灵敏度系数矩阵

灵敏度系数矩阵指：单位面积的匀场片对DSV区域的每个点产生的z轴磁场强度。列数为NMR尺寸reshape到一维后的长度，如40x24=960，行数为DSV区域采样点reshape到一维后的长度，如33x32=1056。灵敏度系数矩阵的尺寸为(1056,960)
$$
B_z(r_i,z_i)=\frac{\mu{_0}m_z}{4\pi}\left(\frac{3z_i^2}{(r_i^2+z_i^2)^{2.5}}-\frac{1}{(r_i^2+z_i^2)^{1.5}}\right)
$$

这里的$r_i$和$z_i$分别表示P点和Q点在x0y平面的距离和z轴方向上的距离。$r_i^2+z_i^2$表示P点和Q点欧式距离的平方。



灵敏度系数矩阵$A$第$i$行$j$列的值表示，第$j$片匀场片对DSV区域第$i$个采样点z方向上的磁场强度。

匀场片分布$x$是长度为960的向量(960,1)，表示NMR每个匀场片的体积。

因此，$Ax$表示
$$
\left[
	\begin{matrix}
		a_{1,1} & a_{1,2} & a_{1,3} & ... & a_{1,960} \\
		a_{2,1} & a_{2,2} & a_{2,3} & ... & a_{2,960} \\
		a_{3,1} & a_{3,2} & a_{3,3} & ... & a_{3,960} \\
		a_{4,1} & a_{4,2} & a_{4,3} & ... & a_{4,960} \\
		... & ... & ... & ... & ...\\
		a_{1056,1} & a_{1056,2} & a_{1056,3} & ... & a_{1056,960} \\
	\end{matrix} 
\right].
\left[
	\begin{matrix}
		x_{1}\\
		x_{2}\\
		x_{3}\\
		x_{4}\\
		...\\
		x_{960}\\
	\end{matrix} 
\right]=
\left[
	\begin{matrix}
		a_{1,1}x_1+a_{1,2}x_2+a_{1,3}x_3+...+a_{1,960}x_{960}\\
		a_{2,1}x_1+a_{2,2}x_2+a_{2,3}x_3+...+a_{2,960}x_{960}\\
		a_{3,1}x_1+a_{3,2}x_2+a_{3,3}x_3+...+a_{3,960}x_{960}\\
		a_{4,1}x_1+a_{4,2}x_2+a_{4,3}x_3+...+a_{4,960}x_{960}\\
		...\\
		a_{1056,1}x_1+a_{1056,2}x_2+a_{1056,3}x_3+...+a_{1056,960}x_{960}\\
	\end{matrix}
\right]
$$
$Ax$向量长度为(1056,1)，表示960个匀场片区域对第$i$个DSV区域采样点z方向的磁场强度。



### 核心公式

整个匀场优化不等式的核心条件：不均匀度低于$\varepsilon$
$$
H=\frac{|B_m+Ax-B_t|}{B_t}<\varepsilon
$$
$B_m$表示DSV区域采样z方向的裸磁场强度，本项目中size为(1056,1)，$Ax$表示匀场片对DSV区域产生的z方向磁场强度，$B_t$表示DSV区域采样点z方向的目标磁场强度，$|B_m+Ax-B_t|$表示DSV区域采样点实际磁场强度偏差值，与$B_t$相除得到DSV区域的不均匀度向量$H$。

去绝对值：
$$
\begin{cases}
	Ax-(1+\varepsilon)B_t\leq -B_m \\
	-Ax+(1-\varepsilon)B_t\leq B_m

\end{cases}
$$
这里的$B_t$有两种处理方法：视作变量还是定值？如果将$B_t$作为定值来处理(取裸磁场强度数据$B_m$的均值$B_{avg}$)：

|                                  |  待匀场数据1(33x32)  |  待匀场数据2(33x24)   |
| :------------------------------: | :------------------: | :-------------------: |
|          匀场前不均匀度          | 374.1361988011714ppm | 673.6933666712266ppm  |
| Bt为定值时最低不均匀度(非整数解) | 27.96210756804757ppm | 21.643841045385646ppm |

如果将$B_t$作为变量来处理，由于$\varepsilon$是变量，所以目标会变成一个非线性规划问题(存在$\varepsilon{B_t}$)，这时就有两种方案：

1. 固定$\varepsilon$，将$B_t$作为变量，通过步长遍历的方法手动搜索最低的$\varepsilon$值（因为单纯形法可以解出线性规划问题的全局最优解，如果手动设置的$\varepsilon$无法找到全局最优解，则证明无法将不均匀度降低到该值）
2. 将问题转化成非线性规划问题，求得最小的$\epsilon$



A如果采用第一种方案：

|                                  |  待匀场数据1(33x32)   |  待匀场数据2(33x24)  |
| :------------------------------: | :-------------------: | :------------------: |
|          匀场前不均匀度          | 374.1361988011714ppm  | 673.6933666712266ppm |
| Bt为变量时最低不均匀度(非整数解) | 13.760798195898063ppm | 20.08963664231919ppm |

这里需要注意Python的Scipy库的linprog函数，在使用”highs“方法时，某些输入下会出现程序卡死的bug，目前尝试了timeout等方法但无法作用，可能是直接进程卡死了，暂时没有解决，拟通过回避bug参数的方法来避免。（改变Config文件中stepNum的值）

