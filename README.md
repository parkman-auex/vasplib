# vasplib

![Vasplib_intro](/fig/vasplib_intro.png "Magic Gardens")

parkman@buaa.edu.cn

#### 介绍
{**以下是 vasplib 说明，您可以帮忙更改**
vasplib是 BUAA-TCMP 推出的基于 matlab 的第一性原理计算结果分析处理，紧束缚模型与低能有效kp模型的构建，调控与计算。为拓扑材料，磁性材料以及量子输运相关科研工作者提供稳定快速、简易高效的程序开发框架
无论是个人、团队、或是高校，都能够用 vasplib实...。

#### 软件架构
外部程序：spglib
matlab工具箱：
vasplib基本函数，vasplib类，HR类，HK类，Htrig类

#### 安装教程

~/Documents/MATLAB/vasplib 
addpath 放入matlab目录

#### 使用说明

直接使用；编写脚本或命令行使用

#### 短期目标
- [ ] plot类的构建 完美继承matlab的基础绘图的handle类的情况下 完成所有绘图函数的构建（参考<https://github.com/plotly/plotly_matlab><https://github.com/masumhabib/PlotPub> 外部箭头类内部化）
- [ ] 拓扑电路绘图切片能带
- [ ] HR2Hckt
- [ ] Clifford 类的建立 一般性地对哈密顿量进行Clifford拆分
- [ ] 强化 HR Htrig HK 对称性分析功能
- [ ] 向量化 作对称性 Htrig HK


#### 发展目标

- [ ] TB 模型中各种场的加入，拉伸以及缺陷，散射！(PYBINDING)
- [ ] HR addSOC功能的修复 
- [ ] Oper Class 的特征标构建 （waiting）
- [ ] 晶体场 Pycrysrtalfield 的复制（原理已懂）-> HCEF的构建 ->结合 vasplib.experiment 对接实验 （2022）
- [ ] 构建 各类的 机器学习拟合DFT能带的程序 （2021.07.08-07.12）完成80%
- [ ] （待选）学习牛派的半经典方法，后续接入VASP的WAVECAR读取，先做个平面波的方法，再做个VASP的借口
- [ ] （待选）神经网络的matlab实现，对接DFT或拓扑以及拓扑量子化学 （）
- [ ] （待选）自旋模型的对称性构建，自动生成 简单ED模型求解基本强关联特性
- [ ] （待选）DMFT+VASP （Tpython什么）接口或复制
- [ ] （待选）解析DFT磁性结果自动生成格点模型自动给出蒙卡，DMRG，XTRG等（取决于李老师的Qmatgen发展情况）
- [ ] （待选）HatreeFock与DFT原理实现（toy-vasp）


#### 已完成
- [x] Symmetry TB in HR-(2021.07.12-07.20)(Htrig working;HK done)；完成
- [x] 多带Berry phase in HK；Then -> HR Htrig （2021.07）(fyang) 完成
- [x]  主类 vasplib 绘图与输入输出的提炼总结 （2021.07~08）(已基本构建) 完成
- [x]  HK Htrig 用Gamma矩阵（tau sigma）矩阵展开
- [x]  Oper 自动生成对称矩阵; 
- [x]  fyang:学习牛派的半经典方法(BCP,BCD等)，融入vasplib体系 （2021.08）程序尽量General一点;
- [x]  Hckt的建立
- [x] 分析 hspice的数据 建立拓扑电路绘图的基本函数
#### 参与贡献

* parkman 曾旭涛 框架设计 主体程序
* cchen 陈聪 程序启发 Slater-Koster HTRCI  待续
* fyang 杨帆 OpenMX QE等第一性原理软件接口 Mag 半经典（BC，BCD，BCP）待续
* leo 刘彬彬 转角石墨烯 程序测试与开发 待续
* xlsheng 胜献雷 指导、开创，评价与反馈

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request
