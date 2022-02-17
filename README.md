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

<div align=center>
<img src="https://latex.codecogs.com/svg.image?\left[\begin{matrix}A_{11}&A_{12}&A_{13}&A_{14}\\A_{21}&A_{22}&A_{23}&A_{24}\\A_{31}&A_{32}&A_{33}&A_{34}\\A_{41}&A_{42}&A_{43}&A_{44}\\\end{matrix}&space;\right]\left[\begin{matrix}A_{11}\\A_{12}\\A_{13}\\A_{14}\\...\\A_{42}\\A_{43}\\A_{44}\\\end{matrix}&space;\right]&space;" title="\left[\begin{matrix}A_{11}&A_{12}&A_{13}&A_{14}\\A_{21}&A_{22}&A_{23}&A_{24}\\A_{31}&A_{32}&A_{33}&A_{34}\\A_{41}&A_{42}&A_{43}&A_{44}\\\end{matrix} \right]\left[\begin{matrix}A_{11}\\A_{12}\\A_{13}\\A_{14}\\...\\A_{42}\\A_{43}\\A_{44}\\\end{matrix} \right] " />
</div>


(33,32)->(1056,)  (33,24)->(796,)



##### 灵敏度系数矩阵

灵敏度系数矩阵指：单位面积的匀场片对DSV区域的每个点产生的z轴磁场强度。列数为NMR尺寸reshape到一维后的长度，如40x24=960，行数为DSV区域采样点reshape到一维后的长度，如33x32=1056。灵敏度系数矩阵的尺寸为(1056,960)

<div align=center>
  <img src="https://latex.codecogs.com/svg.image?B_z(r_i,z_i)=\frac{\mu{_0}m_z}{4\pi}\left(\frac{3z_i^2}{(r_i^2&plus;z_i^2)^{2.5}}-\frac{1}{(r_i^2&plus;z_i^2)^{1.5}}\right)" title="B_z(r_i,z_i)=\frac{\mu{_0}m_z}{4\pi}\left(\frac{3z_i^2}{(r_i^2+z_i^2)^{2.5}}-\frac{1}{(r_i^2+z_i^2)^{1.5}}\right)" />
</div>

这里的![](https://latex.codecogs.com/svg.image?r_i)和![](https://latex.codecogs.com/svg.image?z_i)分别表示P点和Q点在x0y平面的距离和z轴方向上的距离。![](https://latex.codecogs.com/svg.image?r_i^2+z_i^2)表示P点和Q点欧式距离的平方。

灵敏度系数矩阵![](https://latex.codecogs.com/svg.image?A)第![](https://latex.codecogs.com/svg.image?i)行![](https://latex.codecogs.com/svg.image?j)列的值表示，第![](https://latex.codecogs.com/svg.image?j)片匀场片对DSV区域第![](https://latex.codecogs.com/svg.image?i)个采样点z方向上的磁场强度。

匀场片分布![](https://latex.codecogs.com/svg.image?x)是长度为960的向量(960,1)，表示NMR每个匀场片的体积。

因此，![](https://latex.codecogs.com/svg.image?Ax)表示

<div align=center>
<img src="https://latex.codecogs.com/svg.image?\left[\begin{matrix}&space;&space;a_{1,1}&space;&&space;a_{1,2}&space;&&space;a_{1,3}&space;&&space;...&space;&&space;a_{1,960}&space;\\&space;&space;a_{2,1}&space;&&space;a_{2,2}&space;&&space;a_{2,3}&space;&&space;...&space;&&space;a_{2,960}&space;\\&space;&space;a_{3,1}&space;&&space;a_{3,2}&space;&&space;a_{3,3}&space;&&space;...&space;&&space;a_{3,960}&space;\\&space;&space;a_{4,1}&space;&&space;a_{4,2}&space;&&space;a_{4,3}&space;&&space;...&space;&&space;a_{4,960}&space;\\&space;&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...\\&space;&space;a_{1056,1}&space;&&space;a_{1056,2}&space;&&space;a_{1056,3}&space;&&space;...&space;&&space;a_{1056,960}&space;\\&space;\end{matrix}&space;\right].\left[&space;\begin{matrix}&space;&space;x_{1}\\&space;&space;x_{2}\\&space;&space;x_{3}\\&space;&space;x_{4}\\&space;&space;...\\&space;&space;x_{960}\\&space;\end{matrix}&space;\right]=\left[&space;\begin{matrix}&space;&space;a_{1,1}x_1&plus;a_{1,2}x_2&plus;a_{1,3}x_3&plus;...&plus;a_{1,960}x_{960}\\&space;&space;a_{2,1}x_1&plus;a_{2,2}x_2&plus;a_{2,3}x_3&plus;...&plus;a_{2,960}x_{960}\\&space;&space;a_{3,1}x_1&plus;a_{3,2}x_2&plus;a_{3,3}x_3&plus;...&plus;a_{3,960}x_{960}\\&space;&space;a_{4,1}x_1&plus;a_{4,2}x_2&plus;a_{4,3}x_3&plus;...&plus;a_{4,960}x_{960}\\&space;&space;...\\&space;&space;a_{1056,1}x_1&plus;a_{1056,2}x_2&plus;a_{1056,3}x_3&plus;...&plus;a_{1056,960}x_{960}\\&space;\end{matrix}\right]" title="\left[\begin{matrix} a_{1,1} & a_{1,2} & a_{1,3} & ... & a_{1,960} \\ a_{2,1} & a_{2,2} & a_{2,3} & ... & a_{2,960} \\ a_{3,1} & a_{3,2} & a_{3,3} & ... & a_{3,960} \\ a_{4,1} & a_{4,2} & a_{4,3} & ... & a_{4,960} \\ ... & ... & ... & ... & ...\\ a_{1056,1} & a_{1056,2} & a_{1056,3} & ... & a_{1056,960} \\ \end{matrix} \right].\left[ \begin{matrix} x_{1}\\ x_{2}\\ x_{3}\\ x_{4}\\ ...\\ x_{960}\\ \end{matrix} \right]=\left[ \begin{matrix} a_{1,1}x_1+a_{1,2}x_2+a_{1,3}x_3+...+a_{1,960}x_{960}\\ a_{2,1}x_1+a_{2,2}x_2+a_{2,3}x_3+...+a_{2,960}x_{960}\\ a_{3,1}x_1+a_{3,2}x_2+a_{3,3}x_3+...+a_{3,960}x_{960}\\ a_{4,1}x_1+a_{4,2}x_2+a_{4,3}x_3+...+a_{4,960}x_{960}\\ ...\\ a_{1056,1}x_1+a_{1056,2}x_2+a_{1056,3}x_3+...+a_{1056,960}x_{960}\\ \end{matrix}\right]" />
</div>

![](https://latex.codecogs.com/svg.image?Ax)向量长度为(1056,1)，表示960个匀场片区域对第![](https://latex.codecogs.com/svg.image?i)个DSV区域采样点z方向的磁场强度。



### 核心公式

整个匀场优化不等式的核心条件：不均匀度低于![](https://latex.codecogs.com/svg.image?\varepsilon)

<div align=center>
  <img src="https://latex.codecogs.com/svg.image?H=\frac{|B_m&plus;Ax-B_t|}{B_t}<\varepsilon" title="H=\frac{|B_m+Ax-B_t|}{B_t}<\varepsilon" />
</div>

![](https://latex.codecogs.com/svg.image?B_m)表示DSV区域采样z方向的裸磁场强度，本项目中size为(1056,1)，![](https://latex.codecogs.com/svg.image?Ax)表示匀场片对DSV区域产生的z方向磁场强度，![](https://latex.codecogs.com/svg.image?B_t)表示DSV区域采样点z方向的目标磁场强度，![](https://latex.codecogs.com/svg.image?|B_m+Ax-B_t|)表示DSV区域采样点实际磁场强度偏差值，与![](https://latex.codecogs.com/svg.image?B_t)相除得到DSV区域的不均匀度向量![](https://latex.codecogs.com/svg.image?H)。

去绝对值：

<div align=center>
  <img src="https://latex.codecogs.com/svg.image?\begin{cases}&space;Ax-(1&plus;\varepsilon)B_t\leq&space;-B_m&space;\\&space;-Ax&plus;(1-\varepsilon)B_t\leq&space;B_m\end{cases}" title="\begin{cases} Ax-(1+\varepsilon)B_t\leq -B_m \\ -Ax+(1-\varepsilon)B_t\leq B_m\end{cases}" />
</div>
#### 线性规划

这里的![](https://latex.codecogs.com/svg.image?B_t)有两种处理方法：视作变量还是定值？如果将![](https://latex.codecogs.com/svg.image?B_t)作为定值来处理(取裸磁场强度数据![](https://latex.codecogs.com/svg.image?B_m)的均值![](https://latex.codecogs.com/svg.image?B_{avg})，后续简称为FTMF)：

|                                  |  待匀场数据1(33x32)  |  待匀场数据2(33x24)   |
| :------------------------------: | :------------------: | :-------------------: |
|          匀场前不均匀度          | 374.1361988011714ppm | 673.6933666712266ppm  |
| Bt为定值时最低不均匀度(非整数解) | 27.96210756804757ppm | 21.643841045385646ppm |

如果将![](https://latex.codecogs.com/svg.image?B_t)作为变量来处理，由于![](https://latex.codecogs.com/svg.image?\varepsilon)是变量，所以目标会变成一个非线性规划问题(存在![](https://latex.codecogs.com/svg.image?\varepsilon{B_t}))。此时可以通过手动调整![](https://latex.codecogs.com/svg.image?\varepsilon)的方法来保留问题的线性。

固定![](https://latex.codecogs.com/svg.image?\varepsilon)，将![](https://latex.codecogs.com/svg.image?B_t)作为变量，通过步长遍历的方法手动搜索最低的![](https://latex.codecogs.com/svg.image?\varepsilon)值（因为单纯形法可以解出线性规划问题的全局最优解，如果手动设置的![](https://latex.codecogs.com/svg.image?\varepsilon)无法找到全局最优解，则证明无法将不均匀度降低到该值，后续简称为OTMF)

|                                  |  待匀场数据1(33x32)   |  待匀场数据2(33x24)   |
| :------------------------------: | :-------------------: | :-------------------: |
|          匀场前不均匀度          | 374.1361988011714ppm  | 673.6933666712266ppm  |
| Bt为变量时最低不均匀度(非整数解) | 14.108702021411231ppm | 20.249867969269378ppm |

这里需要注意Python的Scipy库的linprog函数，在使用”highs“方法时，某些输入下会出现程序卡死的bug，目前尝试了timeout等方法但无法作用，可能是直接进程卡死了，暂时没有解决，拟通过回避bug参数的方法来避免。（改变Config文件中stepNum的值）

<div align=center>
  <img src="https://github.com/YueZhangX/Passive-shimming-algorithm/blob/main/resources/ImageFiles/6009map1_hom%26x_linear.png" title="数据1的FTMF和OTMF对比"/>
  <br/>
  <strong>
    图1 "6009map1"数据下FTMF和OTMF的对比
  </strong>
  <br/>
  <img src="https://github.com/YueZhangX/Passive-shimming-algorithm/blob/main/resources/ImageFiles/9006map1_hom%26x_linear.png" title="数据2的FTMF和OTMF对比"/>
  <br/>
  <strong>
    图2 "9006map1"数据下FTMF和OTMF的对比
  </strong>
</div>



从上图可以看出，OTMF的方法通过轻微变动Bt值，可以大幅度减小同样不均匀度下所需的最小匀场片厚度，同时也能找到更低的不均匀度匀场方案。

而在耗时方面：由于FTMF是单次线性高维线性规划，而OTMF是多次高维线性规划的搜索算法，因此OTMF相较于FTMF耗时更多，但是也在1分钟之内，是可接受范围。

硬件环境： CPU: Apple M1, 核总数:8（4性能和4能效）

|                                      | 待匀场数据1(33x32)  | 待匀场数据2(33x24)  |
| :----------------------------------: | :-----------------: | :-----------------: |
|             FTMF方法耗时             | 1.6761548519134521s | 1.7174441814422607s |
| OTMF方法耗时(SHIMMING_STEP = 200000) | 31.515713930130005s | 16.34144902229309s  |

#### 非线性规划

在线性规划方案中，如果![](https://latex.codecogs.com/svg.image?B_t)和![](https://latex.codecogs.com/svg.image?\varepsilon)同时为变量，由于不等式中存在![](https://latex.codecogs.com/svg.image?\varepsilon{B_t})，使得问题转化成非线性规划问题。如果直接将本问题当成非线性规划问题来处理，会发现程序无法在限定时间内(timeout)求得有效解。



如果将![](https://latex.codecogs.com/svg.image?B_t)仍然考虑为定值，将目标函数从线性函数变为基于L1范数的最小二乘非线性函数：

<div align=center>
  <img src="https://latex.codecogs.com/svg.image?f_{min}=\|Ax-b\|_2&plus;\lambda{_1}\vert&space;Ax-b\vert&space;_1&plus;\lambda{_2}\vert&space;x-x_0\vert{_1}" title="f_{min}=\|Ax-b\|_2+\lambda{_1}\vert Ax-b\vert _1+\lambda{_2}\vert x-x_0\vert{_1}" />
</div>

测试结果来看，该方法较为不稳定，相较于线性规划算法运算时间较长，效果不理想，且线性规划算法可以设置权值来平衡匀场片厚度与不均匀度。

|                        | 待匀场数据1(33x32)  |  待匀场数据2(33x24)  |
| :--------------------: | :-----------------: | :------------------: |
|       L1范式耗时       | 384.8193471431732s  |  732.1910841464996s  |
| L1范式可匀最低不均匀度 | 37.6942220181026ppm | 25.86320973965214ppm |

无论是匀场效果还是耗时，都非常不稳定。因此本项目舍弃了基于L1范式的非线性规划方法。

#### 整数规划

本问题下如果直接使用整数规划并不可行，原因在于x的维度过高，导致整数规划的复杂度过高，现有的算力无法快速求得全局最优解。如果改用启发式算法来直接进行整数规划，会发现第一组数据结果会收敛到52ppm左右（理论最优为14ppm），效果较差。使用遗传算法效果如图所示，其中种群大小为50，代数为1000。



|                                        |  待匀场数据1(33x32)   |  待匀场数据2(33x24)   |
| :------------------------------------: | :-------------------: | :-------------------: |
|             匀场前不均匀度             | 374.1361988011714ppm  | 673.6933666712266ppm  |
| 直接使用遗传算法做整数规划最小不均匀度 | 52.135839175566744ppm | 56.925301679805685ppm |
| 直接使用遗传算法做整数规划匀场片总厚度 |        4.6224m        |  4.756600000000001m   |

<div align=center>
  <img src="https://github.com/YueZhangX/Passive-shimming-algorithm/blob/main/resources/ImageFiles/6009map1_GA_1-Stage.png" title="数据1的遗传算法图示"/>
  <br/>
  <strong>
    图1 "6009map1"数据下的遗传算法图示
  </strong>
  <br/>
  <img src="https://github.com/YueZhangX/Passive-shimming-algorithm/blob/main/resources/ImageFiles/9006map1_GA_1-Stage.png" title="数据2的遗传算法图示"/>
  <br/>
  <strong>
    图2 "9006map1"数据下的遗传算法图示
  </strong>
</div>

### 实际问题求解

在实际操作中，匀场片的个数必然是整数片，然而由于整数规划方法在本问题下的困难性，多数传统方案都不选择整数规划算法作为匀场工具。在实际操作时，仅对非整数规划结果进行简单的取整处理，这样的处理方法会对不均匀度产生一定影响。

|                              |  待匀场数据1(33x32)   |  待匀场数据2(33x24)   |
| :--------------------------: | :-------------------: | :-------------------: |
|        匀场前不均匀度        | 374.1361988011714ppm  | 673.6933666712266ppm  |
|     线性规划最低不均匀度     | 14.108702021411231ppm | 20.249867969269378ppm |
|     线性规划匀场片总厚度     |  2.8755590444607533m  |  2.163487999115938m   |
|   取整后不均匀度(向上取整)   | 20.249867969269378ppm | 25.228274074270562ppm |
| 取整后匀场片总厚度(向上取整) |        2.8814m        |        2.1677m        |
|   取整后不均匀度(向下取整)   | 18.28996957466761ppm  | 24.86446623664483ppm  |
| 取整后匀场片总厚度(向下取整) |  2.8699000000000003m  |        2.1599m        |

如第一组数据，对非整数部分向上取整时，不均匀度增加了43.53%，向下取整时，增加了29.64%；

第二组数据向上和向下取整，不均匀度分别增加了24.58%和22.79%。

#### 基于启发式算法的2-Stage整数匀场方案

为了得到用于实际操作的整数匀场解，本项目提出了一种基于启发式算法的2-Stage整数匀场方案，该方案解决了高维数据强约束下直接使用传统整数规划算法计算机无法求解的问题，相较于直接使用启发式算法，能更快收敛，得到更优的匀场解。

1. 使用混合线性规划算法求得初始最优解。此处使用“highs”库，与插点法得到的非整数线性规划算法效果相似，但结果是稀疏矩阵，即在实际操作中只需要对少部分内腔插入匀场片，而不用对几乎所有的内腔进行操作。同时，基于混合线性规划得到的非整数x较少，基于本方法得到的整数解不均匀度更低。
2. 基于step1得到的初始解，使用启发式算法搜索最优整数解。在每个![](https://latex.codecogs.com/svg.image?[x_0-x_t,x_0+x_t])范围内寻找最优整数解。

|                                               |  待匀场数据1(33x32)   |  待匀场数据2(33x24)   |
| :-------------------------------------------: | :-------------------: | :-------------------: |
|                匀场前不均匀度                 | 374.1361988011714ppm  | 673.6933666712266ppm  |
|     step1混合线性规划的不均匀度(非整数解)     | 14.108702021411231ppm | 20.249867969269378ppm |
|   step1混合线性规划的匀场片总厚度(非整数解)   |  2.8755590444607533m  |  2.163487999115938m   |
|   step2启发式算法整数规划后的不均匀度(xt=1)   | 15.136251212857795ppm |  21.1689179306274ppm  |
| step2启发式算法整数规划后的匀场片总厚度(xt=1) |        2.8851m        |  2.1773000000000002m  |

此处默认使用遗传算法，且xt默认设置为1。

可以发现，使用该方法后，大幅减小了最优解取整后的不均匀度损失，极大地提高了实际匀场效果。



如果对xt进行调整，可以发现可匀场的不均匀度也会发生变化：

|                |  待匀场数据1(33x32)   |  待匀场数据2(33x24)   |
| :------------: | :-------------------: | :-------------------: |
| 匀场前不均匀度 | 374.1361988011714ppm  | 673.6933666712266ppm  |
|      xt=1      | 15.136251212857795ppm |  21.1689179306274ppm  |
|      xt=2      | 15.057665279621647ppm | 20.297134729296808ppm |
|      xt=3      | 15.744356823468413ppm | 21.54439023328372ppm  |
|      xt=4      | 17.238976301748362ppm | 21.578010558231068ppm |

可以发现，在xt取到2的时候，通常能得到最低的不均匀度。

选择不同的启发式算法效果也有差别(xt=2)：

|      算法名      |      遗传算法      |    差分进化算法     |     粒子群算法      |    模拟退火算法    |
| :--------------: | :----------------: | :-----------------: | :-----------------: | :----------------: |
|   最小不均匀度   |     15.5242ppm     |     19.2280ppm      |     19.8203ppm      |     17.2931ppm     |
| 最小匀场片总厚度 |      2.9129m       | 2.8911883281557045m | 2.8931170924993785m | 2.892387300063346m |
|     算法耗时     | 42.08870029449463s | 21.57760500907898s  | 1.928091049194336s  | 6.786157846450806s |

|      算法名      |      遗传算法       |    差分进化算法    |     粒子群算法     |    模拟退火算法     |
| :--------------: | :-----------------: | :----------------: | :----------------: | :-----------------: |
|   最小不均匀度   |     20.9651ppm      |     25.9972ppm     |     25.7552ppm     |     22.9704ppm      |
| 最小匀场片总厚度 |       2.2184m       | 2.196548588243077m | 2.19728483477403m  | 2.1960629903420914m |
|     算法耗时     | 41.302629232406616s | 20.05046796798706s | 1.513017177581787s | 5.760963201522827s  |

不同的启发式算法效果各不相同，综合来看遗传算法效果更优，但算法耗时与最小匀场片厚度较高；模拟退火算法的最低匀场不均匀度较高，但算法耗时更少，且方案的匀场片总厚度更大。

## 代码结构

### 主函数

主函数为main.py文件，里面包含了多项test函数，用于不同的匀场实验，test文件不需要任何形参，即开即用。

### 参数文件

参数文件存放于Config.py，里面包含了各项常量，基本上需要调整项目的配置参数，只需要调整Config即可，无需在项目内再做修改。

### 裸磁场加载器

裸磁场加载器类存放于Dataloader.py，读取dam文件中的裸磁场数据

### 磁共振仪器参数

磁共振类存放于NMR.py，用于保存磁共振仪器的各项参数

### 灵敏度系数矩阵

灵敏度系数矩阵类存放于SensitivityCoefficientMatrix.py，用于计算灵敏度系数矩阵

### 求解器

求解器类存放于Solver.py，用于测试各项匀场算法
