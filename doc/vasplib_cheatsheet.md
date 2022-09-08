+++
title = "vasplib cheatsheet"
date = 2022-02-02T18:20:31+08:00
draft = false
author = "parkman"
description = "MATLAB vasplib"
toc = true
math = true
featured = true
categories = [
	'vasplib',
]
tags = [
	'MATLAB'
]
series = [
]
aliases = []
images = [
]
+++
vasplib 快捷命令一览
<!--more-->

## vasplib

vasplib 和第一性原理计算对接，Tools function

### bandplot

{{< expand "Dependence" >}} POSCAR KPOINTS EIGENCAR DOSCAR {{< /expand >}}
直接画能带图

``` MATLAB
bandplot();
```

间接画能带图

``` MATLAB
EIGENCAR = EIGENVAL_read();
bandplot(EIGENCAR);
```

调整能量范围

``` MATLAB
bandplot(EIGENCAR,[-3,3]);
```

指定颜色(键值对)

``` MATLAB
bandplot(EIGENCAR,[-3,3],'Color','r');
```

指定依赖的POSCAR和KPOINTS(键值对)

``` MATLAB
bandplot(EIGENCAR,[-3,3],'Color','r','KPOINTS','KPOINTS','POSCAR','POSCAR');
```

画多个能带

``` MATLAB
bandplot({EIGENCAR1,EIGENCAR2},[-3,3]);
%% or
[fig,ax] = bandplot(EIGENCAR1,[-3,3]);
[fig,ax] = bandplot(EIGENCAR2,[-3,3],'fig',fig,'ax',ax);
```

### bandcompare

比较两个能带

``` MATLAB
vasplib.bandplot(EIGENCAR1,EIGENCAR2,[-3,3]);
```

### export .svg

(输出向量图)
600 指像素 temp为文件名

``` MATLAB
set(gcf,'Renderer','Painters'); print(gcf,'-dsvg','-r600','temp.svg')
```

### export .png

(输出位图)
600 指像素 temp为文件名

``` MATLAB
export_fig temp.png -r600
```

### pbandplot

直接画投影能带图(vaspkit 213)
{{< expand "Dependence" >}} POSCAR KPOINTS EIGENCAR DOSCAR PBAND_X.dat {{< /expand >}}

``` MATLAB
pbandplot();
%% or
vasplib.pbandplot();
```

在有WEIGHTCAR(和EIGENCAR结构一致)的情况下 画投影能带图

``` MATLAB
vasplib.pbandplot(WEIGHTCAR,EIGENCAR,'Ecut',[-3,3]);
```

调整能量范围

``` MATLAB
bandplot(EIGENCAR,[-3,3]);
```

指定颜色(键值对)

``` MATLAB
bandplot(EIGENCAR,[-3,3],'Color','r');
```

指定依赖的POSCAR和KPOINTS(键值对)

``` MATLAB
bandplot(EIGENCAR,[-3,3],'Color','r','KPOINTS','KPOINTS','POSCAR','POSCAR');
```

### dosplot(TO BE MODIFIED)

直接画分波态密度(vaspkit 111 113)
{{< expand "Dependence" >}} tdos.dat PDOS_X.dat {{< /expand >}}

``` MATLAB
dosplot();
```

不看d轨道

``` MATLAB
dosplot('d',false);
```

### gen .win .in (TO BE MODIFIED)

生成 wannier90 wt symmetry 模版文件
{{< expand "Dependence" >}} POSCAR KPOINTS EIGENCAR DOSCAR {{< /expand >}}

spinless

``` MATLAB
init_v2w(0);
```

spinful

``` MATLAB
init_v2w(0);
```

判断能带选择区间

``` MATLAB
init_v2w('check',true);
```

快速生成，以选取石墨一炔12个$p_z$轨道为例

``` MATLAB
init_v2w('Projector_list',ones(1,18));
```

### POSCAR_read

{{< expand "Dependence" >}} POSCAR {{< /expand >}}
读取POSCAR信息

``` MATLAB
POSCAR_read;
%% or
[Rm,sites,Atom_name,Atom_num]=POSCAR_readin();
```

### plot BZ

``` MATLAB
vasplib.BZplot(Rm);
```

More
{{< expand "Dependence" >}} KPOINTS {{< /expand >}}

``` MATLAB
vasplib.BZplot(Rm,'KPOINTS','KPOINTS');
```

### plot fermi_arc data from wt

{{< expand "Dependence" >}}  arc.dat_l/arc.dat_r {{< /expand >}}

```MATLAB
[DOSCAR_l,X,Y]  = DOSCAR_read('arc.dat_l',[301 301],'arc'); % 301 is your NK1 NK2
heatplot(DOSCAR_l,Y,X);
```

### plot surf_dos data from wt

{{< expand "Dependence" >}}  dos.dat_l/dos.dat_r {{< /expand >}}

```MATLAB
[DOSCAR_l,X,Y]  = DOSCAR_read('dos.dat_l',[1000 402]); 
%1000 is your w_num 603 is your total K points along surf kpath
heatplot(DOSCAR_l,Y,X);
```

## HR

**Construction**

### From wannier90_hr.dat

{{< expand "Dependence" >}} wannier90_hr.dat {{< /expand >}}

``` MATLAB
HRobj = HR.from_wannier90('wannier90_hr.dat');
```

### SKTB

POSCAR 需要特别的设置，详见POSCAR 设置：[vasplib_POSCAR](https://parkman.top/en/posts/matlab/vasplib/vasplib_poscar/)
{{< expand "Dependence" >}} POSCAR {{< /expand >}}

#### sparse

##### layer C pz or 3d C pz 体系

``` MATLAB
% parameters = [d0,a0,VppS_0,VppP_0,delta0,GlobalShift];
clear parameters ; % test
SKTB =@(parameters) HR.from_POSCAR_SE("POSCAR",...
    'Type',Type,...
    'r_max',d_max,...
    'level_cut',-1,...
    'search_range',search_range,...
    'Accuracy',1e-3,...
    'chiral',false,...
    'onsite',false,...
    'deltarule',1,...
    'Rd',[parameters(1),parameters(2)],...
    'para',struct('strvar',["VppS","VppP"],'numvar',[parameters(3),parameters(4)],'delta',  parameters(5)),...
    'E0',parameters(6));
```

调参并最终获得 H_TBSK 这个 *HRobj*

``` MATLAB
d0 = 3.1251;
a0 = 1.0128;
VppS_0  = 0.8217;%-0.48%0.48;
VppP_0  = -5.5798;%*10;%2.7;
delta0 = 0.6501;%0.45255;%0.39337;%0.45255;
Fermi = -1.2715;
parameters = [d0,a0,VppS_0,VppP_0,delta0,Fermi];
H_TBSK  = SKTB(parameters);
```

### I/O

#### 导入晶格基失

{{< expand "Dependence" >}} POSCAR {{< /expand >}}

``` MATLAB
HRobj = HRobj.input_Rm('POSCAR');
```

#### 导入orb坐标

{{< expand "Dependence" >}} POSCAR {{< /expand >}}

``` MATLAB
HRobj = HRobj.input_orb_struct('POSCAR');
%% or
HRobj = HRobj < 'POSCAR';
```

#### 载入kpath

{{< expand "Dependence" >}} KPOINTS {{< /expand >}}

``` MATLAB
HRobj = HRobj.kpathgen3D('POSCAR');
%% or
HRobj = HRobj < 'KPOINTS';
```

#### 查看 符号化变量

``` MATLAB
HRobj.symvar_list;
```

#### list hopping

``` MATLAB
HRobj.list;
```

list certain Rvector hopping

``` MATLAB
HRobj.list([0,0,0]);
```

list certain i->j hopping

``` MATLAB
HRobj.list([1,1]);
```

check onsite mat

``` MATLAB
HRobj.HnumL(:,:,HRobj.Line_000);
%% or
HRobj.HcoeL(:,:,HRobj.Line_000);
```

**Modification**

**Solve**

### EIGENCAR & WAVECAR

#### default

``` MATLAB
[EIGENCAR,WAVECAR] = HRobj.EIGENCAR_gen();
```

#### solve on certain kpoints

``` MATLAB
[EIGENCAR,WAVECAR] = HRobj.EIGENCAR_gen('klist',[0,0,0]);
```

#### solve on given klist

``` MATLAB
[EIGENCAR,WAVECAR] = HRobj.EIGENCAR_gen('klist',klist);
```
### band 3D

``` MATLAB
[EIGENCAR_3D,klist1,klist2] = HRobj.EIGENCAR_gen_3D([101,101],[-0.5,-0.5,0;1,0,0;0 1 0]);
```

``` MATLAB
[fig,ax] = Grpahene_TB_n.BZplot(HRobj.Rm,'mode','2D');
vasplib.bandplot_3d(EIGENCAR_3D,klist1,klist2,'fig',fig,'ax',ax);
```


## HK

## Htrig

Initializing Htrig class

``` MATLAB
Graphene = Htrig(2);
```

Defining Htrig with given Hamiltonian

``` MATLAB
Graphene = Graphene +...
    Trig(-t*cos(delta(1,:)*k),sigma_x)+...
    Trig( t*sin(delta(1,:)*k),sigma_y)+...
    Trig(-t*cos(delta(2,:)*k),sigma_x)+...
    Trig( t*sin(delta(2,:)*k),sigma_y)+...
    Trig(-t*cos(delta(3,:)*k),sigma_x)+...
    Trig( t*sin(delta(3,:)*k),sigma_y);
```

Band plot with Htrig

``` MATLAB
a = 1;t =1;
Graphene_n  = Graphene.Subsall();
Graphene_n = Graphene_n<'POSCAR';
Graphene_n = Graphene_n <'KPOINTS';
EIGENCAR = Graphene_n.EIGENCAR_gen();
[klist_l,kpoints_l,kpoints_name] = Graphene_n.kpath_information();
bandplot(EIGENCAR,[-3,3],klist_l,kpoints_l,kpoints_name,'title','Graphene','Color','r');
```

**Other functions such as slab band can be viewed in "vasplib/example/example08-1-HTrig-class/Htrig_C3/Graphene_Htrig.mlx".**

## Oper
