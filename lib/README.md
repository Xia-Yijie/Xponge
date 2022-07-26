# Xponge

## 简介
Xponge是由Python编写的轻量化能高度自定义的分子动力学模拟的前后处理工具。Python版本需大于3.6，安装直接使用pip即可
```bash
pip install Xponge
```
安装成功以后，在Python脚本中直接`import`即可使用
```python
import Xponge
```
Xponge最基本的功能除了依赖于Python的标准库外，还依赖了numpy、pubchempy，这两个包将会在pip install Xponge时将会自动安装，而一些更复杂的功能可能包括其他模块，这些模块需要自行安装。

Xponge目前还处于早期开发版本，此处介绍一些基本功能，不代表最终版本效果。

## 力场生成
### A. 现有力场
#### A1. 现有残基，无初始坐标文件

- 引入现有力场
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB
```
会得到包含参考文献的打印信息
```plain
Reference for ff14SB.py:
  James A. Maier, Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling
    ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB
    Journal of Chemical Theory and Computation 2015 11 (8), 3696-3713
    DOI: 10.1021/acs.jctc.5b00255
```

- 查看现有残基
```python
print(Xponge.ResidueType.types)
```
得到结果：
```plain
{'ACE': Type of Residue: ACE, 'ASH': Type of Residue: ASH, 'CYM': Type of Residue: CYM, 'GLH': Type of Residue: GLH, 'LYN': Type of Residue: LYN, 'NME': Type of Residue: NME, 'NALA': Type of Residue: NALA, 'ALA': Type of Residue: ALA, 'CALA': Type of Residue: CALA, 'NARG': Type of Residue: NARG, 'ARG': Type of Residue: ARG, 'CARG': Type of Residue: CARG, 'NASN': Type of Residue: NASN, 'ASN': Type of Residue: ASN, 'CASN': Type of Residue: CASN, 'NASP': Type of Residue: NASP, 'ASP': Type of Residue: ASP, 'CASP': Type of Residue: CASP, 'NCYS': Type of Residue: NCYS, 'CYS': Type of Residue: CYS, 'CCYS': Type of Residue: CCYS, 'NCYX': Type of Residue: NCYX, 'CYX': Type of Residue: CYX, 'CCYX': Type of Residue: CCYX, 'NGLN': Type of Residue: NGLN, 'GLN': Type of Residue: GLN, 'CGLN': Type of Residue: CGLN, 'NGLU': Type of Residue: NGLU, 'GLU': Type of Residue: GLU, 'CGLU': Type of Residue: CGLU, 'NGLY': Type of Residue: NGLY, 'GLY': Type of Residue: GLY, 'CGLY': Type of Residue: CGLY, 'NHID': Type of Residue: NHID, 'HID': Type of Residue: HID, 'CHID': Type of Residue: CHID, 'NHIE': Type of Residue: NHIE, 'HIE': Type of Residue: HIE, 'CHIE': Type of Residue: CHIE, 'NHIP': Type of Residue: NHIP, 'HIP': Type of Residue: HIP, 'CHIP': Type of Residue: CHIP, 'NHE': Type of Residue: NHE, 'HYP': Type of Residue: HYP, 'CHYP': Type of Residue: CHYP, 'NILE': Type of Residue: NILE, 'ILE': Type of Residue: ILE, 'CILE': Type of Residue: CILE, 'NLEU': Type of Residue: NLEU, 'LEU': Type of Residue: LEU, 'CLEU': Type of Residue: CLEU, 'NLYS': Type of Residue: NLYS, 'LYS': Type of Residue: LYS, 'CLYS': Type of Residue: CLYS, 'NMET': Type of Residue: NMET, 'MET': Type of Residue: MET, 'CMET': Type of Residue: CMET, 'NPHE': Type of Residue: NPHE, 'PHE': Type of Residue: PHE, 'CPHE': Type of Residue: CPHE, 'NPRO': Type of Residue: NPRO, 'PRO': Type of Residue: PRO, 'CPRO': Type of Residue: CPRO, 'NSER': Type of Residue: NSER, 'SER': Type of Residue: SER, 'CSER': Type of Residue: CSER, 'NTHR': Type of Residue: NTHR, 'THR': Type of Residue: THR, 'CTHR': Type of Residue: CTHR, 'NTRP': Type of Residue: NTRP, 'TRP': Type of Residue: TRP, 'CTRP': Type of Residue: CTRP, 'NTYR': Type of Residue: NTYR, 'TYR': Type of Residue: TYR, 'CTYR': Type of Residue: CTYR, 'NVAL': Type of Residue: NVAL, 'VAL': Type of Residue: VAL, 'CVAL': Type of Residue: CVAL, 'HIS': Type of Residue: HIE, 'NHIS': Type of Residue: NHIE, 'CHIS': Type of Residue: CHIE}
```
- 使用现有残基构建初始结构

残基类（Xponge.ResidueType）、残基（Xponge.Residue）、分子（Xponge.Molecule）之间，可以使用`+`进行连接获得分子（Xponge.Molecule），也可以`*`整数（int）进行多次自我连接。连接规则由刚才import的力场文件决定。
```python
protein = ACE + ALA * 3 + NME  
```

- 保存输出

残基类（Xponge.ResidueType）、残基（Xponge.Residue）、分子（Xponge.Molecule）可以使用`Save_PDB`、`Save_SPONGE_Input`、`Save_Mol2`、`Save_NPZ`分别保存为对应的文件格式。
```python
Save_PDB(protein, "protein.pdb")
Save_SPONGE_Input(protein, "protein")
```
使用VMD观看保存的PDB
```bash
vmd protein.pdb
```
或利用vmd的SPONGE插件
```bash
vmd -sponge_mass ./protein_mass.txt -sponge_crd ./protein_coordinate.txt
```
可获得如下结果

![输入图片说明](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/1.png)

- 动力学模拟

此处仅作示例，因此保存的分子较小，人为修改生成的protein_coordinate.txt文件的最后一行，将盒子的长度改为50.0，以避免过小的周期性边界造成的影响。

修改好后，新建"mdin.txt"文件，并填写
```plain
basic test of Xponge
    mode = NVT
    thermostat = Andersen_thermostat
    dt = 2e-3
    constrain_mode = SHAKE
    cutoff = 8.0
    step_limit = 5000
    coordinate_in_file = protein_coordinate.txt
    angle_in_file = protein_angle.txt
    bond_in_file = protein_bond.txt
    charge_in_file = protein_charge.txt
    dihedral_in_file = protein_dihedral.txt
    exclude_in_file = protein_exclude.txt
    LJ_in_file = protein_LJ.txt
    mass_in_file = protein_mass.txt
    nb14_in_file = protein_nb14.txt
    residue_in_file = protein_residue.txt
```
使用SPONGE运行，获得如下结果
```plain
---------------------------------------------------------------------------------------
        step =         1000,         time =        2.000,  temperature =       305.35,
   potential =        47.21,           LJ =        -3.73,          PME =      -232.74,
     nb14_LJ =        26.62,      nb14_EE =       196.30,         bond =        12.74,
       angle =        16.63,     dihedral =        32.51,
---------------------------------------------------------------------------------------
        step =         2000,         time =        4.000,  temperature =       194.94,
   potential =        51.01,           LJ =        -2.02,          PME =      -229.18,
     nb14_LJ =        30.27,      nb14_EE =       188.48,         bond =        12.04,
       angle =        22.31,     dihedral =        30.75,
---------------------------------------------------------------------------------------
        step =         3000,         time =        6.000,  temperature =       246.57,
   potential =        49.60,           LJ =        -4.42,          PME =      -222.75,
     nb14_LJ =        31.28,      nb14_EE =       189.09,         bond =        11.10,
       angle =        16.47,     dihedral =        30.80,
---------------------------------------------------------------------------------------
        step =         4000,         time =        8.000,  temperature =       273.81,
   potential =        55.53,           LJ =        -4.54,          PME =      -231.54,
     nb14_LJ =        30.79,      nb14_EE =       196.70,         bond =         9.60,
       angle =        18.88,     dihedral =        29.40,
---------------------------------------------------------------------------------------
        step =         5000,         time =       10.000,  temperature =       344.12,
   potential =        48.44,           LJ =        -5.35,          PME =      -230.99,
     nb14_LJ =        29.09,      nb14_EE =       198.53,         bond =         6.32,
       angle =        19.48,     dihedral =        31.79,
---------------------------------------------------------------------------------------
```
> 部分内容在后续例子中相同，后续例子只在有必要的时候展示重复操作的结果部分。
 - 该部分完整Python代码
 

```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB
protein = ACE + ALA * 3 + NME
Save_PDB(protein, "protein.pdb")
Save_SPONGE_Input(protein, "protein")
```

#### A2. 现有残基，有初始坐标

目前推荐是从PDB文件中读入

- 引入现有力场
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB
```

- 读入现有文件
> [示例pdb文件](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/0.15_80_10_pH6.5_6lzg.result.pdb)：由H++补氢并添加了SSBOND信息的新冠病毒Spike蛋白与ACE2受体结合的PDB
```python
protein = loadpdb("0.15_80_10_pH6.5_6lzg.result.pdb")
```
- 该部分完整Python代码
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB
protein = loadpdb("0.15_80_10_pH6.5_6lzg.result.pdb")
Save_PDB(protein, "protein.pdb")
Save_SPONGE_Input(protein, "protein")
```

> 注意事项
> 1. PDB文件格式是标准化的，修改前请确认理解PDB的格式
> 2. 从Protein Data Bank下载的PDB文件不包含氢，目前SPONGE不支持自动补氢，请手动或使用其他工具（如[H++](http://biophysics.cs.vt.edu/)）补氢。补氢时需要注意氢的命名需与标准相同
> 3. H++产生的PDB不包含SSBOND，需自行将原始PDB中SSBOND部分添加至H++产生的PDB中，Xponge能自动将二硫键的半胱氨酸残基更名。另外，Xponge能够自动根据质子判断组氨酸的残基名

#### A3. 非现有单独残基，已知信息

只需要构建好残基类（Xponge.ResidueType），将非现有残基转化为现有残基

以gaff力场下的自定义水为例

##### 方案1 从mol2文件中构建

编写mol2文件
```mol2
@<TRIPOS>MOLECULE
TP3
    3     2     1     0     1 
SMALL
USER_CHARGES
@<TRIPOS>ATOM
  1 O       0.000000    0.000000    0.000000 oh    1 WAT    -0.8340 ****
  2 H1      0.957200    0.000000    0.000000 ho    1 WAT     0.4170 ****
  3 H2     -0.239988    0.926627    0.000000 ho    1 WAT     0.4170 ****
@<TRIPOS>BOND
    1     1     2 1
    2     1     3 1
@<TRIPOS>SUBSTRUCTURE
      1  WAT              1 ****               0 ****  **** 
```
然后在Python中
```python
import Xponge
import Xponge.forcefield.AMBER.gaff
TP3 = loadmol2("WAT.mol2")
Save_SPONGE_Input(TP3, "tp3")
```
> Xponge会根据mol2文件中的残基名（此例中的"WAT"）在内部添加残基类（Xponge.ResidueType），而loadmol2最后会返回的（此例中的"TP3"）是包含该mol2中所有残基的分子（Xponge.Molecule）。参与PDB中残基识别只与残基类（Xponge.ResidueType）有关。
##### 方案2 Python代码构建
```python
#构建新的空残基类
import Xponge
import Xponge.forcefield.AMBER.gaff
WAT = Xponge.ResidueType(name = "WAT")
#添加原子信息
WAT.Add_Atom(name = "O", atom_type = Xponge.AtomType.types["oh"], x = 0, y = 0, z = 0)
WAT.O.Update(**{"charge[e]": -0.8340})
WAT.Add_Atom(name = "H1", atom_type = Xponge.AtomType.types["ho"], x = 0.9572, y = 0, z = 0)
WAT.H1.Update(**{"charge[e]": 0.4170})
WAT.Add_Atom(name = "H2", atom_type = Xponge.AtomType.types["ho"], x = -0.239988, y = 0.926627, z = 0)
WAT.H2.Update(**{"charge[e]": 0.4170})
#添加连接信息
WAT.Add_Connectivity(WAT.O, WAT.H1)
WAT.Add_Connectivity(WAT.O, WAT.H2)

Save_SPONGE_Input(WAT, "TP3")
```

#### A4. 非现有单独残基，未知信息
以gaff力场下的2,3-二甲基苯甲酸乙酯为例
##### 方案1 通过IUPAC名或SMILES结构式从PubChem中获取基本的结构
```python
import Xponge

#从PubChem中获取结构
assign = Get_Assignment_From_PubChem("ethyl 2,6-dimethylbenzoate", "name")
#也可以使用SMILES
#assign = Get_Assignment_From_PubChem("CCOC(=O)C1=C(C=CC=C1C)C", "smiles")

#自动推断原子类别
import Xponge.forcefield.AMBER.gaff
assign.Determine_Atom_Type("GAFF")

#保存assignment至mol2文件中
Save_Mol2(assign, "EDF_ASN.mol2")

#通过vmd判断出等价性
equivalence = [[9,10], range(16,22), [3,4], [5,6], [13,14]]
q = assign.Calculate_Charge("RESP", basis = "6-311g*", grid_density = 1, 
    extra_equivalence = equivalence, opt = True)

EDF = assign.To_ResidueType("EDF", q)
Save_Mol2(EDF, "EDF.mol2")
```
最后获得的"EDF.mol2"即可进行下一步计算
```mol2
@<TRIPOS>MOLECULE
EDF
 27 27 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    O    1.826   -0.000    0.361   os     1      EDF  -0.563454
     2   O1    1.320   -0.002   -1.787    o     1      EDF  -0.656809
     3    C   -0.453   -0.000   -0.196   ca     1      EDF  -0.456248
     4   C1   -1.103   -1.216   -0.002   ca     1      EDF   0.283884
     5   C2   -1.101    1.217   -0.003   ca     1      EDF   0.283884
     6   C3   -2.432   -1.197    0.395   ca     1      EDF  -0.284050
     7   C4   -2.431    1.200    0.394   ca     1      EDF  -0.284050
     8   C5   -3.091    0.002    0.592   ca     1      EDF  -0.137956
     9   C6    0.979   -0.001   -0.654    c     1      EDF   1.049933
    10   C7   -0.382   -2.527   -0.219   c3     1      EDF  -0.396791
    11   C8   -0.379    2.527   -0.221   c3     1      EDF  -0.396791
    12   C9    3.218   -0.000    0.060   c3     1      EDF   0.506223
    13  C10    3.968    0.001    1.373   c3     1      EDF  -0.267413
    14    H   -2.953   -2.125    0.550   ha     1      EDF   0.170890
    15   H1   -2.951    2.129    0.548   ha     1      EDF   0.170890
    16   H2   -4.122    0.002    0.899   ha     1      EDF   0.163931
    17   H3   -0.010   -2.607   -1.234   hc     1      EDF   0.118938
    18   H4    0.467   -2.623    0.451   hc     1      EDF   0.118938
    19   H5   -1.041   -3.366   -0.037   hc     1      EDF   0.118938
    20   H6   -0.007    2.605   -1.236   hc     1      EDF   0.118938
    21   H7   -1.037    3.366   -0.041   hc     1      EDF   0.118938
    22   H8    0.470    2.622    0.449   hc     1      EDF   0.118938
    23   H9    3.449   -0.875   -0.533   h1     1      EDF  -0.046089
    24  H10    3.449    0.872   -0.535   h1     1      EDF  -0.046089
    25  H11    5.037    0.001    1.189   hc     1      EDF   0.064160
    26  H12    3.723    0.880    1.958   hc     1      EDF   0.064160
    27  H13    3.723   -0.876    1.960   hc     1      EDF   0.064160
@<TRIPOS>BOND
     1      1      9 1
     2      1     12 1
     3      2      9 1
     4      3      4 1
     5      3      5 1
     6      3      9 1
     7      4      6 1
     8      4     10 1
     9      5      7 1
    10      5     11 1
    11      6      8 1
    12      6     14 1
    13      7      8 1
    14      7     15 1
    15      8     16 1
    16     10     17 1
    17     10     18 1
    18     10     19 1
    19     11     20 1
    20     11     21 1
    21     11     22 1
    22     12     13 1
    23     12     23 1
    24     12     24 1
    25     13     25 1
    26     13     26 1
    27     13     27 1
@<TRIPOS>SUBSTRUCTURE
    1      EDF      1 ****               0 ****  **** 
```

> 目前Xponge不支持自动判断等价原子，因此需要靠vmd来手动判断

![输入图片说明](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/2.png)


> 目前Xponge的assignment只支持只由C、H、O、N、S、P、F、Cl、Br、I元素构成的物质

##### 方案2 其他方式构建assignment的mol2文件

只需将对应的Get_Assignment_From_PubChem更改为Get_Assignment_From_Mol2即可
```python
assign = Get_Assignment_From_Mol2("EDF_ASN.mol2")
```

> 注意，任务（Xponge.assign.Assign）和分子（Xponge.Molecule）均可读mol2，但是它们是不一样的，而且文件中的信息有效性也不同

> 任务（Xponge.assign.Assign）使用Get_Assignment_From_Mol2读取，其中的原子类别为元素符号，且bond部分中的键需要明确键级

> 分子（Xponge.Molecule）使用loadmol2读取，其中的原子类别为力场中的符号，且bond部分中的键不会读取键级

##### 方案3 Python代码构建
关键在于构建任务类(Xponge.assign.Assign)，其他步骤相同。此处以构建水分子为例。
```python
import Xponge
import Xponge.forcefield.AMBER.gaff
assign = Xponge.assign.Assign()
assign.Add_Atom(element = "O", x = 0, y = 0, z = 0)
assign.Add_Atom(element = "H", x = 0.9572, y = 0, z = 0)
assign.Add_Atom(element = "H", x = -0.239988, y = 0.926627, z = 0)
assign.Add_Bond(0,1,1) #0号原子和1号原子之间形成一个单键
assign.Add_Bond(0,2,1) #0号原子和2号原子之间形成一个单键

#注意，在推断原子类型前需要推断环和键的信息，虽然对水没有用
assign.Determine_Ring_And_Bond_Type()

assign.Determine_Atom_Type("GAFF")
q = assign.Calculate_Charge("RESP", extra_equivalence = [[1,2]])

WAT = assign.To_ResidueType("WAT", q)
Save_Mol2(WAT, "WAT.mol2")
```
最终即获得"WAT.mol2":
```mol2
@<TRIPOS>MOLECULE
WAT
 3 2 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    O    0.000    0.000    0.000   oh     1      WAT  -0.799534
     2    H    0.957    0.000    0.000   ho     1      WAT   0.399767
     3   H1   -0.240    0.927    0.000   ho     1      WAT   0.399767
@<TRIPOS>BOND
     1      1      2 1
     2      1      3 1
@<TRIPOS>SUBSTRUCTURE
    1      WAT      1 ****               0 ****  **** 
```
#### A5. 非现有非单独残基
##### 方案1 mol2文件构建
以ff14SB力场下丙氨酸二肽的构建为例
> 实际上，ff14SB力场中包含丙氨酸二肽，此处只为示例

先构建dipeptide.mol2文件
```mol2
@<TRIPOS>MOLECULE
ACE
 22 21 3 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1   H1    0.466   -8.051    1.242   HC     1      ACE   0.112300
     2  CH3    0.289   -7.339    0.436   CT     1      ACE  -0.366200
     3   H2   -0.703   -6.901    0.548   HC     1      ACE   0.112300
     4   H3    0.352   -7.853   -0.524   HC     1      ACE   0.112300
     5    C    1.325   -6.213    0.457    C     1      ACE   0.597200
     6    O    2.209   -6.195    1.311    O     1      ACE  -0.567900
     7    N    1.275   -5.267   -0.433    N     2      ALA  -0.415700
     8    H    0.563   -5.249   -1.150    H     2      ALA   0.271900
     9   CA    2.256   -4.201   -0.413   CX     2      ALA   0.033700
    10   HA    2.199   -3.671    0.538   H1     2      ALA   0.082300
    11   CB    3.668   -4.752   -0.580   CT     2      ALA  -0.182500
    12  HB1    3.744   -5.277   -1.532   HC     2      ALA   0.060300
    13  HB2    4.384   -3.930   -0.561   HC     2      ALA   0.060300
    14  HB3    3.888   -5.443    0.234   HC     2      ALA   0.060300
    15    C    2.009   -3.207   -1.539    C     2      ALA   0.597300
    16    O    1.072   -3.368   -2.318    O     2      ALA  -0.567900
    17    N    2.791   -2.179   -1.682    N     3      NME  -0.415700
    18    H    3.571   -2.011   -1.063    H     3      NME   0.271900
    19  CH3    2.555   -1.233   -2.754   CT     3      NME  -0.149000
    20 HH31    1.686   -1.546   -3.331   H1     3      NME   0.097600
    21 HH32    2.374   -0.244   -2.333   H1     3      NME   0.097600
    22 HH33    3.429   -1.195   -3.405   H1     3      NME   0.097600
@<TRIPOS>BOND
     1      1      2 1
     2      2      3 1
     3      2      4 1
     4      2      5 1
     5      5      6 1
     6      5      7 1
     7      7      8 1
     8      7      9 1
     9      9     10 1
    10      9     11 1
    11      9     15 1
    12     11     12 1
    13     11     13 1
    14     11     14 1
    15     15     16 1
    16     15     17 1
    17     17     18 1
    18     17     19 1
    19     19     20 1
    20     19     21 1
    21     19     22 1
@<TRIPOS>SUBSTRUCTURE
    1      ACE      1 ****               0 ****  **** 
    2      ALA      7 ****               0 ****  **** 
    3      NME     17 ****               0 ****  **** 
```
然后在python中读入力场参数后加载mol2文件即可正常使用该残基，并会自动按照mol2的前后连接关系在后续的PDB读入中连接。
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

protein = loadmol2("dipeptide.mol2")
```
如果不需要使用`+`和`*`构建新分子，则到处就已经满足使用了。如果希望正确使用`+`和`*`，则需要再设置连接信息
```python
res = Xponge.ResidueType.types["ALA"]
res.head = "N"  #残基的头部（与前一个残基连接的原子）atom name是N
res.head_length = 1.3  #与前一个残基连接的键长是1.3埃
res.head_next = "CA" #残基的主链头部第二个原子atom name是CA
res.tail = "C" #残基的尾部（与后一个残基连接的原子）atom name是N
res.tail_length = 1.3 #与后一个残基连接的键长是1.3埃
res.tail_next = "CA" #残基的主链尾部第二个原子atom name是CA

#清空ff14SB中已有的信息，普通构建未有残基不需要
while res.head_link_conditions:
    res.head_link_conditions.pop()
while res.head_link_conditions:
    res.tail_link_conditions.pop()
import numpy as np
#与前一个残基连接时，CA-N-前一个残基的tail原子形成的角为120/180 * np.pi
res.head_link_conditions.append({"atoms":["CA", "N"], "parameter": 120/180 * np.pi})
#与前一个残基连接时，H-CA-N-前一个残基的tail原子形成的二面角为-np.pi
res.head_link_conditions.append({"atoms":["H", "CA", "N"], "parameter": -np.pi})
#与后一个残基连接时，CA-C-后一个残基的head原子形成的角为120/180 * np.pi
res.tail_link_conditions.append({"atoms":["CA", "C"], "parameter": 120/180 * np.pi})
#与后一个残基连接时，O-CA-C-后一个残基的head原子形成的二面角为-np.pi
res.tail_link_conditions.append({"atoms":["O", "CA", "C"], "parameter": -np.pi})     
```
> 刚体连接有6个自由度，分别是1个bond，2个angle，3个dihedral。其中1个bond由键长决定，1个主链dihedral（tail_next - tail - head - head_next）设定为180度防止重叠，剩余4个自由度由前一个残基和后一个残基分别定义一个angle和dihedral来实现

如果希望在PDB中，端基不同，则可产生不同的残基类(Xponge.ResidueType)，如NALA、CALA，然后在全局设置中设置
```python
Xponge.GlobalSetting.Add_PDB_Residue_Name_Mapping("head", "ALA", "NALA")
Xponge.GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", "ALA", "CALA")
```

特别地，对于组氨酸判别不同质子态和半胱氨酸判别二硫键，需要额外设置
```python
Xponge.GlobalSetting.HISMap["DeltaH"] = "HD1"
Xponge.GlobalSetting.HISMap["EpsilonH"] = "HE2"
Xponge.GlobalSetting.HISMap["HIS"].update({"HIS": {"HID":"HID", "HIE":"HIE", "HIP":"HIP"}, 
                                    "CHIS":{"HID":"CHID", "HIE":"CHIE", "HIP":"CHIP"},
                                    "NHIS":{"HID":"NHID", "HIE":"NHIE", "HIP":"NHIP"}})
Xponge.ResidueType.types["CYX"].connect_atoms["ssbond"] = "SG"
```
##### 方案2 Python代码构建
类似A3方案2构建残基类(Xponge.ResidueType)，然后再按A5方案1设置连接信息即可。
### B. 非现有力场
#### B1. 新力场参数
##### 方案1 字符串/文本文件/npz文件读取
以在AMBER力场形式下构建OPLSAA力场的硝基甲烷为例。数据来源于[GMXTPP](https://jerkwin.github.io/prog/gmxtop.html)。
```python
import Xponge
import Xponge.forcefield.AMBER

#构建新的原子类
Xponge.AtomType.New_From_String("""name  mass  charge[e] LJtype
opls_760  14.01   0.54   opls_760
opls_761  16.00  -0.37   opls_761
opls_762  12.01   0.02   opls_762
opls_763  1.008   0.06   opls_763""")

#构建新的LJ参数信息
Xponge.forcefield.BASE.LJ.LJType.New_From_String("""name  sigma[nm]  epsilon[kJ/mol]
opls_760-opls_760 0.3250000    0.5020800
opls_761-opls_761 0.2960000    0.7112800
opls_762-opls_762 0.3500000    0.2761440
opls_763-opls_763 0.2500000    0.0627600  
""")

#构建新的成键信息
#New_From_NPZ待更新

temp_dict =  {"opls_762-opls_763": {"b[nm]":0.109, "k[kJ/mol·nm^-2]": 284512},
              "opls_762-opls_760": {"b[nm]":0.14900, "k[kJ/mol·nm^-2]": 313800},
              "opls_760-opls_761": {"b[nm]":0.1225, "k[kJ/mol·nm^-2]": 460240}}
Xponge.forcefield.BASE.BOND.BondType.New_From_Dict(temp_dict)

Xponge.forcefield.BASE.ANGLE.AngleType.New_From_File("test_angle.txt")

Xponge.forcefield.BASE.DIHEDRAL.ProperType.New_From_String("""name phi0[degree] k[kJ/mol] periodicity  reset
opls_763-opls_762-opls_760-opls_761  0 0 0 1
""")

Xponge.forcefield.BASE.DIHEDRAL.ImproperType.New_From_String("""name phi0[degree] k[kJ/mol] periodicity
opls_761-opls_761-opls_760-X 180.0     43.93200   2
""")

#OPLS力场的nb14的修正系数与amber不同
Xponge.forcefield.BASE.NB14.NB14Type.New_From_String("""name kLJ kee
opls_763-opls_761 0.5 0.5
""")

NIM = Xponge.ResidueType(name = "NIM")
#添加原子信息
NIM.Add_Atom(name = "CT", atom_type = Xponge.AtomType.types["opls_762"], x = 1.107, y = 1.13, z = 1.117)
NIM.Add_Atom(name = "HC1", atom_type = Xponge.AtomType.types["opls_763"], x = 1.143, y = 1.029, z = 1.117)
NIM.Add_Atom(name = "HC2", atom_type = Xponge.AtomType.types["opls_763"], x = 1.143, y = 1.18, z = 1.03)
NIM.Add_Atom(name = "HC3", atom_type = Xponge.AtomType.types["opls_763"], x = 1., y = 1.13, z = 1.117)
NIM.Add_Atom(name = "NO", atom_type = Xponge.AtomType.types["opls_760"], x = 1.156, y = 1.199, z = 1.237)
NIM.Add_Atom(name = "ON1", atom_type = Xponge.AtomType.types["opls_761"], x = 1.177, y = 1.321, z = 1.234)
NIM.Add_Atom(name = "ON2", atom_type = Xponge.AtomType.types["opls_761"], x = 1.177, y = 1.135, z = 1.341)

#添加连接信息
NIM.Add_Connectivity(NIM.CT, NIM.HC1)
NIM.Add_Connectivity(NIM.CT, NIM.HC2)
NIM.Add_Connectivity(NIM.CT, NIM.HC3)
NIM.Add_Connectivity(NIM.CT, NIM.NO)
NIM.Add_Connectivity(NIM.ON1, NIM.NO)
NIM.Add_Connectivity(NIM.ON2, NIM.NO)

Save_PDB(NIM)
Save_SPONGE_Input(NIM)
```
代码中的"test_angle.txt"文件如下：
```plain
name b[degree] k[kJ/mol·rad^-2]
opls_763-opls_762-opls_763  107.800    276.144
opls_763-opls_762-opls_760  105.000    292.880
opls_762-opls_760-opls_761  117.500    669.440
opls_761-opls_760-opls_761  125.000    669.440
```

##### 方案2 AMBER格式读取（parmdat和frcmod）
例如[ff14SB](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/Xponge/forcefield/AMBER/ff14SB.py)和[ff19SB](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/Xponge/forcefield/AMBER/ff19SB.py)的实现，先使用`loadparmdat`
和`loadfrcmod`获得字符串，然后再使用各类型的`New_From_String`即可
```python
import Xponge
import Xponge.forcefield.AMBER
atoms, bonds, angles, propers, impropers, LJs = loadparmdat("parm19.dat")

Xponge.AtomType.New_From_String(atoms)
Xponge.forcefield.BASE.BOND.BondType.New_From_String(bonds)
Xponge.forcefield.BASE.ANGLE.AngleType.New_From_String(angles)
Xponge.forcefield.BASE.DIHEDRAL.ProperType.New_From_String(propers)
Xponge.forcefield.BASE.DIHEDRAL.ImproperType.New_From_String(impropers)
Xponge.forcefield.BASE.LJ.LJType.New_From_String(LJs)


atoms, bonds, angles, propers, impropers, LJs, cmap = loadfrcmod("ff19SB.frcmod")

Xponge.AtomType.New_From_String(atoms)
Xponge.forcefield.BASE.BOND.BondType.New_From_String(bonds)
Xponge.forcefield.BASE.ANGLE.AngleType.New_From_String(angles)
Xponge.forcefield.BASE.DIHEDRAL.ProperType.New_From_String(propers)
Xponge.forcefield.BASE.DIHEDRAL.ImproperType.New_From_String(impropers)
Xponge.forcefield.BASE.LJ.LJType.New_From_String(LJs)

from Xponge.forcefield.BASE import RCMAP
Xponge.forcefield.BASE.RCMAP.CMAP.Residue_Map.update(cmap)
```
##### 方案3 GROMACS格式读取（forcefield.itp）
例如[CHARMM27的蛋白质力场](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/Xponge/forcefield/CHARMM27/protein.py)的实现，先使用`loadffitp`获得参数字典，然后再使用各类型的`New_From_String`或`New_From_Dict`即可
```python
import Xponge
import Xponge.forcefield.CHARMM27
output = loadffitp("forcefield.itp")

Xponge.AtomType.New_From_String(output["atomtypes"])
Xponge.forcefield.BASE.BOND.BondType.New_From_String(output["bonds"])
output["dihedrals"] += "X-X-X-X 0 0 1 0\n"
Xponge.forcefield.BASE.DIHEDRAL.ProperType.New_From_String(output["dihedrals"])
Xponge.forcefield.BASE.LJ.LJType.New_From_String(output["LJ"])
Xponge.forcefield.BASE.UREY_BRADLEY.UreyBradleyType.New_From_String(output["Urey-Bradley"])
Xponge.forcefield.BASE.IMPROPER.ImproperType.New_From_String(output["impropers"])
Xponge.forcefield.BASE.NB14_EXTRA.NB14Type.New_From_String(output["nb14_extra"])
Xponge.forcefield.BASE.NB14.NB14Type.New_From_String(output["nb14"])
Xponge.forcefield.BASE.ACMAP.CMAP.New_From_Dict(output["cmaps"])
```
> `loadffitp`目前功能还很不全，并不一定完全支持函数形式，对于不支持的形式待更新或给出错误提示

> `loadffitp`只能加载力场项的itp文件，也即不包含“[ moleculetype ]”的itp

> `loadffitp`中对于宏的处理不会去寻找gromacs的路径，也即文件中的`#include “xxx.itp”`只会在当前目录下寻找

#### B2. 新力场形式
##### 类型1 非键
非键可以看作本身是原子性质，以目前实现的静电和LJ作用为例
###### 静电
本身实现于Xponge.forcefield.BASE.CHARGE中
```python
import Xponge
#给原子类别添加性质
Xponge.AtomType.Add_Property({"charge":float})

#给电荷添加单位
#三个参数的意义是：性质"charge"的单位是"charge"，程序内部单位是"SPONGE"
Xponge.AtomType.Set_Property_Unit("charge", "charge", "e")

#通过修饰器@Xponge.Molecule.Set_Save_SPONGE_Input使得Save_SPONGE_Input时调用该函数
#被调用的函数接受三个参数，分别是构建的分子、前缀名和保存的路径
@Xponge.Molecule.Set_Save_SPONGE_Input      
def write_charge(self, prefix, dirname):
    towrite = "%d\n"%(len(self.atoms)) #通过self.atoms获取所有原子
    towrite += "\n".join(["%.6f"%(atom.charge * 18.2223) for atom in self.atoms])  #通过atom.charge可以调用上面Add_Property的性质
    f = open(os.path.join(dirname, prefix + "_charge.txt"),"w")
    f.write(towrite)
    f.close()
```
###### LJ作用
本身实现于Xponge.forcefield.BASE.LJ中
```python
import Xponge
#给原子类别添加性质
Xponge.AtomType.Add_Property({"LJtype":str})

#生成一个新的类型(Xponge.Type)的子类，参数分别是名字、性质
LJType = Xponge.Generate_New_Pairwise_Force_Type("LJ", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float})

#设置单位（性质、量纲、程序内部单位）
LJType.Set_Property_Unit("rmin", "distance", "A")
LJType.Set_Property_Unit("sigma", "distance", "A")
LJType.Set_Property_Unit("epsilon", "energy", "kcal/mol")
LJType.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
LJType.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")

#通过修饰器Xponge.GlobalSetting.Add_Unit_Transfer_Function(LJType)使得LJType初始化进行单位转化时调用
@Xponge.GlobalSetting.Add_Unit_Transfer_Function(LJType)
def LJ_Unit_Transfer(self):
    if self.A != None and self.B != None:
        self.sigma = (self.A / self.B) ** (1/6)
        self.epsilon = 0.25 * B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1/12) / 2)
        self.sigma = None

#定义默认的组合规则供内部使用，自己定义的
def Lorentz_Berthelot_For_A(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)

def Lorents_Berthelot_For_B(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)
    
#通过修饰器@Xponge.Molecule.Set_Save_SPONGE_Input使得Save_SPONGE_Input时调用该函数
#被调用的函数接受三个参数，分别是构建的分子、前缀名和保存的路径
@Xponge.Molecule.Set_Save_SPONGE_Input      
def write_LJ(self, prefix, dirname):
    LJtypes = []
    LJtypemap = {}
    for atom in self.atoms:
        if atom.LJtype not in LJtypemap.keys():
            LJtypemap[atom.LJtype] = len(LJtypes)
            LJtypes.append(atom.LJtype)
             
    As = []
    Bs = []
    for i in range(len(LJtypes)):
        LJ_i = LJType.types[LJtypes[i] + "-" + LJtypes[i]]
        for j in range(len(LJtypes)):
            LJ_j = LJType.types[LJtypes[j] + "-" + LJtypes[j]]
            finded = False
            findnames = [LJtypes[i] + "-" + LJtypes[j], LJtypes[j] + "-" + LJtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    LJ_ij = LJType.types[findname]
                    As.append(LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    Bs.append(LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    break
            if not finded:
                    As.append(LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
                    Bs.append(LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
    
    checks = {}
    count = 0
    for i in range(len(LJtypes)):
        check_string_A = ""
        check_string_B = ""
        for j in range(len(LJtypes)):
            check_string_A += "%16.7e"%As[count] + " "
            check_string_B += "%16.7e"%Bs[count] + " "
            count += 1
            
        checks[i] = check_string_A+check_string_B
    
    same_type = { i: i for i in range(len(LJtypes))}
    for i in range(len(LJtypes)-1, -1, -1):
        for j in range(i+1, len(LJtypes)):
            if checks[i] == checks[j]:
                same_type[j] = i

    real_LJtypes = []
    real_As = []
    real_Bs = []
    tosub = 0
    for i in range(len(LJtypes)):
        
        if same_type[i] == i:
            real_LJtypes.append(LJtypes[i])
            same_type[i] -= tosub
        else:
            same_type[i] = same_type[same_type[i]]
            tosub += 1

    for i in range(len(real_LJtypes)):
        LJ_i = LJType.types[real_LJtypes[i] + "-" + real_LJtypes[i]]
        for j in range(i+1):
            LJ_j = LJType.types[real_LJtypes[j] + "-" + real_LJtypes[j]]
            finded = False
            findnames = [real_LJtypes[i] + "-" + real_LJtypes[j], real_LJtypes[j] + "-" + real_LJtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    LJ_ij = LJType.types[findname]
                    real_As.append(LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    real_Bs.append(LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    break
            if not finded:
                    real_As.append(LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
                    real_Bs.append(LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
    
                
    
    towrite = "%d %d\n\n"%(len(self.atoms), len(real_LJtypes))
    count = 0
    for i in range(len(real_LJtypes)):
        for j in range(i+1):
            towrite += "%16.7e"%real_As[count] + " "
            count += 1
        towrite +="\n"
    towrite += "\n"
    
    count = 0
    for i in range(len(real_LJtypes)):
        for j in range(i+1):
            towrite += "%16.7e"%real_Bs[count] + " "
            count += 1
        towrite +="\n"
    towrite += "\n"
    towrite += "\n".join(["%d"%(same_type[LJtypemap[atom.LJtype]]) for atom in self.atoms])
    f = open(os.path.join(dirname, prefix + "_LJ.txt"),"w")
    f.write(towrite)
    f.close()

```

##### 类型2 成键
###### 基于原子判断的(atom-specific)
最简单的例子是AMBER力场的简谐bond
```python
import Xponge

#生成一个新的类型(Xponge.Type)的子类，参数分别是名字、成键关系、性质和是否是必须的
BondType = Xponge.Generate_New_Bonded_Force_Type("bond", "1-2", {"k":float, "b":float}, True)


#设置性质单位
BondType.Set_Property_Unit("k", "energy·distance^-2", "kcal/mol·A^-2")
BondType.Set_Property_Unit("b", "distance", "A")

#保存操作，通过self.bonded_forces["bond"]（"bond"对应于上面的产生类型时输出的字符串）获取相关信息
#而对应的self.bonded_forces["bond"]的成员的属性中，atoms是成键的原子，k、b是上面自己定义的
#self.atom_index是一个字典，原子为key，原子序号为value
@Xponge.Molecule.Set_Save_SPONGE_Input
def write_bond(self, prefix, dirname):
    bonds = []
    for bond in self.bonded_forces["bond"]:
        order = list(range(2))
        if bond.k != 0:
            if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order
            bonds.append("%d %d %f %f"%(self.atom_index[bond.atoms[temp_order[0]]]
            , self.atom_index[bond.atoms[temp_order[1]]], bond.k, bond.b))
    
    if (bonds):
        towrite = "%d\n"%len(bonds)
        bonds.sort(key = lambda x: list(map(int, x.split()[:2])))
        towrite += "\n".join(bonds)
        
        f = open(os.path.join(dirname, prefix + "_bond.txt"),"w")
        f.write(towrite)
        f.close()
```
复杂的则以AMBER力场的二面角为例
```python
import Xponge

#AMBER力场的二面角包括恰当和非恰当两类
#恰当二面角参数分别是名字、成键关系、性质、是否是必须的和可重复的参数（恰当二面角部分参数可重复给定并叠加计算）
ProperType = Generate_New_Bonded_Force_Type("dihedral", "1-2-3-4", {"k":float, "phi0": float, "periodicity":int}, True, ["k", "phi0", "periodicity"])

ProperType.Set_Property_Unit("k", "energy", "kcal/mol")
ProperType.Set_Property_Unit("phi0", "angle", "rad")


ImproperType = Generate_New_Bonded_Force_Type("improper", "1-3-2-3", {"k":float, "phi0": float, "periodicity":int}, False)
#非恰当二面角是复杂的拓扑，需要手动指定
#该矩阵的意思是成键的第i个原子和第j个原子之间的关系，只有上三角的值有效，例如ImproperType.topology_matrix[1][2]指成键的第2个原子和第3个原子是1-2关系
ImproperType.topology_matrix = [[1, 3, 2, 3],
                                [1, 1, 2, 3],
                                [1, 1, 1, 2],
                                [1, 1, 1, 1]]

ImproperType.Set_Property_Unit("k", "energy", "kcal/mol")
ImproperType.Set_Property_Unit("phi0", "angle", "rad")

#确认相同的力。默认是正倒序，例如ProperType中的A-B-C-D和D-C-B-A是同一种力
#如果不是默认情况，需要手动设置，接受两个参数cls（自己的类型，此处即ImproperType）和atom_list（可能是一个"A-B-C-D"的字符串，也可能是[A, B, C, D]的列表）
#输出是一个列表，包含了所有的相同的atom_list
#例如，对于非恰当二面角，A-B-C-D中A、B、D可以随意互换位置
@ImproperType.Set_Same_Force_Function
def Improper_Same_Force(cls, atom_list):
    
    temp = []
    if type(atom_list) == str:
        atom_list_temp = [ atom.strip() for atom in atom_list.split("-")]
        center_atom = atom_list_temp.pop(2)
        for atom_permutation in Xponge.permutations(atom_list_temp):
            atom_permutation = list(atom_permutation)
            atom_permutation.insert(2, center_atom)
            temp.append("-".join(atom_permutation))
    else:
        atom_list_temp = [ atom for atom in atom_list]
        center_atom = atom_list_temp.pop(2)
        for atom_permutation in Xponge.permutations(atom_list_temp):
            atom_permutation = list(atom_permutation)
            atom_permutation.insert(2, center_atom)
            temp.append(atom_permutation)
    return temp

#输出
#对于self.bonded_forces["dihedral"]的元素，因为设置了可叠加的性质，叠加度属性由multiple_numbers获取，每个值由对应的性质+"s"获取
@Molecule.Set_Save_SPONGE_Input
def write_dihedral(self, prefix, dirname):
    dihedrals = []
    for dihedral in self.bonded_forces.get("dihedral", []):
        order = list(range(4))
        if self.atom_index[dihedral.atoms[order[0]]] > self.atom_index[dihedral.atoms[order[-1]]]:
            temp_order = order[::-1]
        else:
            temp_order = order
        for i in range(dihedral.multiple_numbers):            
            if dihedral.ks[i] != 0:
                dihedrals.append("%d %d %d %d %d %f %f"%(self.atom_index[dihedral.atoms[temp_order[0]]]
                , self.atom_index[dihedral.atoms[temp_order[1]]], self.atom_index[dihedral.atoms[temp_order[2]]]
                , self.atom_index[dihedral.atoms[temp_order[3]]], dihedral.periodicitys[i], dihedral.ks[i], dihedral.phi0s[i]))
 
    for dihedral in self.bonded_forces.get("improper", []):
        order = list(range(4))
        if dihedral.k != 0:
            if self.atom_index[dihedral.atoms[order[0]]] > self.atom_index[dihedral.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order
            dihedrals.append("%d %d %d %d %d %f %f"%(self.atom_index[dihedral.atoms[temp_order[0]]]
            , self.atom_index[dihedral.atoms[temp_order[1]]], self.atom_index[dihedral.atoms[temp_order[2]]]
            , self.atom_index[dihedral.atoms[temp_order[3]]], dihedral.periodicity, dihedral.k, dihedral.phi0))
        
    
    if (dihedrals):
        towrite = "%d\n"%len(dihedrals)
        dihedrals.sort(key = lambda x: list(map(float, x.split())))
        towrite += "\n".join(dihedrals)
        
        f = open(os.path.join(dirname, prefix + "_dihedral.txt"),"w")
        f.write(towrite)
        f.close()
```
###### 基于残基判断的(residue-specific)
以ff19SB力场中的cmap为例（RCMAP）
```python
import Xponge
#仍然创建一个种类
CMAP = Xponge.Generate_New_Bonded_Force_Type("residue_specific_cmap", "1-2-3-4-5", {}, False)

#所有信息相当于都自己处理
CMAP.Residue_Map = {}

@Xponge.Molecule.Set_Save_SPONGE_Input
def write_cmap(self, prefix, dirname):
    cmaps = []
    resolutions = []
    used_types = []
    used_types_map = {}
    atoms = []
    for cmap in self.bonded_forces["residue_specific_cmap"]:
        resname = cmap.atoms[2].residue.type.name
        #通过判断2号原子的resname来决定信息
        if resname in CMAP.Residue_Map.keys():
            if CMAP.Residue_Map[resname]["count"] not in used_types_map.keys():
                used_types_map[CMAP.Residue_Map[resname]["count"]] = len(used_types)
                used_types.append(CMAP.Residue_Map[resname]["parameters"])
                resolutions.append(str(CMAP.Residue_Map[resname]["resolution"]))
            cmaps.append("%d %d %d %d %d %d"%(self.atom_index[cmap.atoms[0]], self.atom_index[cmap.atoms[1]],
            self.atom_index[cmap.atoms[2]], self.atom_index[cmap.atoms[3]], 
            self.atom_index[cmap.atoms[4]], used_types_map[CMAP.Residue_Map[resname]["count"]]))
            
    
    if (cmap):
        towrite = "%d %d\n"%(len(cmaps), len(resolutions))
        towrite += " ".join(resolutions) + "\n\n"
        
        for i in range(len(used_types_map)):
            resol = int(resolutions[i])
            for j in range(resol):
                for k in range(resol):
                    towrite += "%f "%used_types[i][j * 24 + k]
                towrite += "\n"
            towrite += "\n"
        
        towrite += "\n".join(cmaps)
        
        f = open(os.path.join(dirname, prefix + "_cmap.txt"),"w")
        f.write(towrite)
        f.close()
```
##### 类型3 虚拟原子
虚拟原子是被当作一种成键作用对待的。目前添加新的虚拟原子需要修改Xponge.forcefield.BASE.VIRTUAL_ATOM
```python
import Xponge

#创建成键相互作用，其中参数包含atomN，N从0到(依赖的原子数-1)，作为判断依赖的原子
VirtualType2 = Xponge.Generate_New_Bonded_Force_Type("vatom2", "1", {"atom0":int, "atom1":int, "atom2":int, "k1":float, "k2":float}, False)
#在Xponge.GlobalSetting.VirtualAtomTypes注册：key是名字，value是依赖的原子数
Xponge.GlobalSetting.VirtualAtomTypes["vatom2"] = 3

@GlobalSetting.Molecule.Set_Save_SPONGE_Input
def write_virtual_atoms(self, prefix, dirname):
    vatoms = []
    for vatom in self.bonded_forces["vatom2"]:
        vatoms.append("2 %d %d %d %d %f %f"%(self.atom_index[vatom.atoms[0]],
            self.atom_index[vatom.atoms[0]] + vatom.atom0, 
            self.atom_index[vatom.atoms[0]] + vatom.atom1,
            self.atom_index[vatom.atoms[0]] + vatom.atom2,
            vatom.k1, vatom.k2))
    
    if (vatoms):
        towrite = ""
        towrite += "\n".join(vatoms)
        
        f = open(os.path.join(dirname, prefix + "_vatom.txt"),"w")
        f.write(towrite)
        f.close()
```

#### B3. 新的原子类型推断
```python
from Xponge import assign

#注册新的推测规则
GAFF = assign.Judge_Rule("GAFF")

#通过GAFF.Add_Judge_Rule(原子种类)进行判断
#程序会按添加顺序从上到下依次判断，直到判断True为止，返回此时的原子种类
#判断函数接受两个参数：需要判断的原子的编号i和判断任务Assign
#Assign.Atom_Judge(i, "O1"): 判断i号原子是不是氧元素且只有一根键
@GAFF.Add_Judge_Rule("o")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O1")

#Assign.Atom_Judge(i, "O2"): 判断i号原子是不是氧元素且只有两根键
#"RG3" in Assign.atom_marker[i].keys()：判断"RG3"有没有在i号原子的标志里，也即是不是在三元环上
#标志由Assign.Determine_Ring_And_Bond_Type()产生，key是标志，value是标志的数量
@GAFF.Add_Judge_Rule("op")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("oq")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2") and "RG4" in Assign.atom_marker[i].keys()

#更复杂的，寻找连接的原子是否是氢
#Assign.bonds[i]是一个字典，key是连接原子的编号，value是成键的键级
@GAFF.Add_Judge_Rule("oh")
def temp(i, Assign):
    tofind = False
    if Assign.Atom_Judge(i, "O2") or Assign.Atom_Judge(i, "O3"):
        for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "H1"):
                tofind = True
                break    
    return tofind

@GAFF.Add_Judge_Rule("os")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2")
```
定义好后，即可对Assign使用`Determine_Atom_Type("GAFF")`进行原子种类推断
#### B4. 新的力场组合
新的力场组合除了加载新的力场形式和力场参数外，还需要设置一些"描述性"的力场设置，例如
AMBER力场中，默认的LJ组合规则、14作用设置、排除到4号连接原子
```python
from Xponge import *
from Xponge.forcefield.BASE import CHARGE, MASS, LJ, BOND, ANGLE, DIHEDRAL, NB14, VIRTUAL_ATOM, EXCLUDE
import os

AMBER_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorents_Berthelot_For_B

NB14.NB14Type.New_From_String(r"""
name    kLJ     kee
X-X     0.5     0.833333
""")
EXCLUDE.Exclude(4)
```
CHARMM27力场中，默认的LJ组合规则、二面角设置、不需要的力场判断、排除到4号连接原子
```python
from Xponge import *
from Xponge.forcefield.BASE import CHARGE, MASS, LJ, BOND, DIHEDRAL, NB14, NB14_EXTRA, UREY_BRADLEY, IMPROPER, VIRTUAL_ATOM, ACMAP, EXCLUDE
import OS

CHARMM27_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorents_Berthelot_For_B

GlobalSetting.Set_Invisible_Bonded_Forces(["improper"])
DIHEDRAL.ProperType.New_From_String(r"""
name        k reset  phi0 periodicity
X-X-X-X     0 0      0    0
""")
EXCLUDE.Exclude(4)
```
## 结构处理
### A. 内坐标更改
#### A1. 更改键长
`Impose_Bond`可以更改键长
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

t = ACE + NME

Save_Mol2(t, "imposing.mol2")

Impose_Bond(t, t.residues[0].C, t.residues[1].N, 5)

Save_Mol2(t, "imposed.mol2")

#不同残基之间可以impose_bond任意两个原子，按残基分别移动
Impose_Bond(t, t.residues[0].C, t.residues[1].CH3, 5)

Save_Mol2(t, "imposed2.mol2")

#同一个残基内没有键连的原子、成环互有叠加的不能impose_bond，因为不知道怎么移动
#Impose_Bond(t, t.residues[0].C, t.residues[0].H2, 5)
#AssertionError
```
上述产生的mol2文件在vmd中观察

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/3.png)

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/4.png)

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/5.png)
#### A2. 更改键角
`Impose_Angle`可以更改键角
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

t = ACE + NME

#不同残基之间可以impose_bond任意两个原子，按残基分别移动
#Impose_Angle要求第2、3个原子之间符合Impose_Bond的要求
Impose_Angle(t, t.residues[0].C, t.residues[1].N, t.residues[1].CH3, 3.1415926 / 2)

Save_Mol2(t, "imposed.mol2")
```
产生的mol2文件在vmd中观察

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/6.png)
#### A3. 更改二面角
`Imporse_Dihedral`可以更改二面角
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

t = ACE + ALA * 10 + NME

Save_Mol2(t, "imposing.mol2")
#Impose_Dihedral要求第2、3个原子之间符合Impose_Bond的要求
for i in range(1,len(t.residues)-1):
    head = t.residues[i-1]
    res = t.residues[i]
    tail = t.residues[i+1]
    Impose_Dihedral(t, head.C, res.N, res.CA, res.C, -3.1415926/3)
    Impose_Dihedral(t, res.N, res.CA, res.C, tail.N, -3.1415926/3)

Save_Mol2(t, "imposed.mol2")
```
改变前的mol2文件在vmd中观察（Draw Method分别使用Line、Ribbon）

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/7.png)

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/8.png)

改变后的mol2文件在vmd中观察（Draw Method分别使用Line、Ribbon）

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/9.png)

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/10.png)
### B. 溶剂与离子添加
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB
import Xponge.forcefield.AMBER.tip3p

t = NALA + ARG + NME
c = round(t.charge)

#添加溶剂盒子，距离边界x<0的5埃，y<0的10埃，z<15的埃，x>0的30埃，y>0的25埃，z>0的20埃
Process_Box(t, WAT, [5,10,15,30,25,20])

#也可以简单的提供一个数
#Process_Box(t, WAT, 30)

#替代，将残基名为"WAT"的部分替代为钾离子和氯原子
Ion_Replace(t, lambda res: res.type.name == "WAT", {CL:30 + c, K:30})

#重新排序
#非必要，只为好看
t.residues.sort(key = lambda residue: {"CL":2, "K":1, "WAT":3}.get(residue.type.name, 0))

#打印出电荷，确保没错
print(t.charge)


Save_PDB(t, "test.pdb")
Save_SPONGE_Input(t, "test")
```

### C. 重复结构产生
暂未实现

### D. 主轴旋转
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

t = ACE + ALA * 100 + NME

Save_Mol2(t,"before.mol2")

Molecule_Rotate(t)

Save_Mol2(t,"after.mol2")
```
旋转前图像（从z轴看过去）

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/11.png)

旋转后图像（从z轴看过去）

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/12.png)

旋转后图像稍旋转观察角度

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/13.png)

### E. 质量重分配
```python
import Xponge
import Xponge.forcefield.AMBER.ff14SB

t = ACE + ALA + NME

import Xponge.forcefield.AMBER.tip3p
Process_Box(t, WAT, [5,10,15,30,25,20])

HMass_Repartition(t)

Save_SPONGE_Input(t, "test")
```

### F. 残基突变
暂未实现

### G. 结构堆叠
暂未实现

## 后处理
### A. 轨迹分析
调用MDAnalysis库进行计算，需自行安装MDAnalysis
```python
#从Xponge中import Universe来加载
from Xponge.analysis.MDAnalysis import Universe
u = Universe("test.pdb", "mdcrd.dat", "mdbox.txt")

#下面都是MDAnalysis的API
O = u.select_atoms("resname WAT and name O")
from MDAnalysis.analysis import rdf as RDF
rdf = RDF.InterRDF(O, O)
rdf.run()
import matplotlib.pyplot as plt
plt.plot(rdf.results.bins[1:], rdf.results.rdf[1:])
plt.show()
```
最终画出来的图像为

![imposing.mol2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/14.png)


### B. 格式转化
使用`python -m Xponge -h`可以获得更多信息

#### B1. dat2nc
将SPONGE的.dat轨迹文件转化为AMBER的.nc文件，详情见`python -m Xponge dat2nc -h`

#### B2. gro2crd
将GROMACS的.gro坐标文件转化为SPONGE的_coordinate.txt文件，详情见`python -m Xponge gro2crd -h`

#### B3. nc2rst7
将AMBER的.nc重开文件转化为AMBER的可读重开文件，详情见`python -m Xponge nc2rst7 -h`

#### B4. maskgen
调用VMD，选择原子，生成对应的原子序号，详情见`python -m Xponge maskgen -h`

#### B5. exgen
输入SPONGE的bond-like, angle-like, dihedral-like和约束文件来获得排除文件，详情见`python -m Xponge exgen -h`

#### B6. dat1frame
从SPONGE的.dat文件中抽出一帧作为SPONGE的_coordinate.txt文件，详情见`python -m Xponge dat1frame -h`


### C. FEP处理
暂未整合
