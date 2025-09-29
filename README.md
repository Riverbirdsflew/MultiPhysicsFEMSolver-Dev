# Kratos多物理场仿真框架 - 热-浓度-结构-相变耦合应用

<p align=center><img height="72.125%" width="72.125%" src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/kratos.png"></p>

[![License][license-image]][license] [![C++][c++-image]][c++standard] [![Github CI][Nightly-Build]][Nightly-link] [![DOI][DOI-image]][DOI] [![GitHub stars][stars-image]][stars] [![Twitter][twitter-image]][twitter] [![Youtube][youtube-image]][youtube]

_KRATOS Multiphysics_ ("Kratos") 是一个用于构建并行、多学科仿真软件的框架，旨在实现模块化、可扩展性和高性能。Kratos使用C++编写，并具有广泛的Python接口。

**Kratos** 在BSD-4许可下是**免费**的，甚至可以在商业软件中使用。

## 项目概述

本项目是基于Kratos Multiphysics框架的扩展，特别添加了**热-浓度-结构-相变四场耦合应用**（ThermDiffStructPhaseApplication）。该应用专门用于钢材热处理过程的多物理场耦合仿真，实现了热传导、碳扩散、结构力学和相变过程的强耦合计算。

### 核心特性

- **四场强耦合**：在单个迭代步内实现热场、浓度场、结构场和相变场的强耦合计算
- **七相钢材模型**：完整描述钢材的微观组织演化，包括奥氏体、铁素体、渗碳体、珠光体、上贝氏体、下贝氏体和马氏体
- **四种相变过程**：支持奥氏体转变、共析转变、贝氏体转变和马氏体转变
- **统一相变管理**：通过MaterialPhaseTransformationSystem统一管理所有相变过程
- **自包含模块**：无需依赖其他应用模块，实现完整的热处理仿真功能

## 新增功能 - ThermDiffStructPhaseApplication

### 多物理场耦合

- **热场**：温度分布和热传导，考虑相变潜热
- **浓度场**：碳浓度分布和扩散，考虑相间碳配分
- **结构场**：应力应变和变形，考虑相变应变、TRIP和蠕变
- **相变场**：相变动力学、相变塑性、蠕变、潜热、相变塑性、相变应变

### 七相钢材模型

1. 奥氏体 (Austenite)
2. 铁素体 (Ferrite) 
3. 渗碳体 (Cementite)
4. 珠光体 (Pearlite)
5. 上贝氏体 (Upper Bainite)
6. 下贝氏体 (Lower Bainite)
7. 马氏体 (Martensite)

### 四种相变过程

1. **奥氏体转变**：加热过程中的奥氏体化
2. **共析转变**：奥氏体向珠光体的分解，偏离共析点还形成铁素体+渗碳体
3. **贝氏体转变**：中温区的贝氏体形成
4. **马氏体转变**：快冷过程中的马氏体形成

### 相变模型

#### 相变动力学模型
- Avrami模型：适用于均匀成核和生长的相变
- ISO模型：等温相变动力学
- JMAK模型：Johnson-Mehl-Avrami-Kolmogorov模型
- KME模型：Koistinen-Marburger方程（马氏体转变）
- TTT模型：基于TTT曲线的相变动力学
- UDef模型：用户自定义相变动力学模型

#### 相变塑性TRIP模型
- Leblond模型
- Desalos模型

#### 配分模型
- Eutectoid模型：共析分解过程
- SubBn模型：贝氏体转变过程中

#### 蠕变模型
- Maxwell Linear模型
- Maxwell Power模型

## 主要组件架构

### 1. PhaseTransformationSystem (相变系统核心类)
位置: `phase_transformation/phase_transformation_system.h/cpp`

**职责**: 相变过程的总控制类，管理所有相变计算

**关键特性**:
- 七相钢材管理
- 四种相变过程协调
- 子模型管理：动力学模型、配分模型、TRIP模型、蠕变模型

### 2. ThermoDiffStructPhaseElement (多场耦合单元)
位置: `custom_elements/thermo_diff_struct_phase_element.h/cpp`

**职责**: 四场耦合的有限元单元

**继承关系**: 继承自 `BaseSolidElement`

### 3. FourFieldCoupledNewtonRaphsonStrategy (四场耦合求解策略)
位置: `custom_strategies/four_field_coupled_newton_raphson_strategy.h`

**职责**: 实现四场强耦合的迭代求解

**继承关系**: 继承自 `ImplicitSolvingStrategy`

## 使用方法

### 1. 编译应用
```bash
# 在Kratos根目录下
mkdir build && cd build
cmake .. -DKRATOS_APPLICATIONS="ThermDiffStructPhaseApplication"
make
```

### 2. Python接口使用
```python
import KratosMultiphysics as KM
import KratosMultiphysics.ThermDiffStructPhaseApplication as TDS

# 创建模型部件
model_part = KM.ModelPart("MainModelPart")

# 添加必要的变量
model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
model_part.AddNodalSolutionStepVariable(KM.CONCENTRATION)
model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
# ... 添加相变相关变量
```

## 许可证

本项目基于Kratos Multiphysics，采用BSD许可证。您可以自由使用、修改和分发本软件，但必须保留原始版权声明。

[license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/license.txt

[c++-image]: https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B
[c++standard]: https://isocpp.org/std/the-standard

[Nightly-Build]: https://github.com/KratosMultiphysics/Kratos/workflows/Nightly%20Build/badge.svg
[Nightly-link]: https://github.com/KratosMultiphysics/Kratos/actions?query=workflow%3A%22Nightly+Build%22

[DOI-image]: https://zenodo.org/badge/DOI/10.5281/zenodo.3234644.svg
[DOI]: https://doi.org/10.5281/zenodo.3234644

[stars-image]: https://img.shields.io/github/stars/KratosMultiphysics/Kratos?label=Stars&logo=github
[stars]: https://github.com/KratosMultiphysics/Kratos/stargazers

[twitter-image]: https://img.shields.io/twitter/follow/kratosmultiphys.svg?label=Follow&style=social
[twitter]: https://twitter.com/kratosmultiphys

[youtube-image]: https://badges.aleen42.com/src/youtube.svg
[youtube]:https://www.youtube.com/@kratosmultiphysics3578