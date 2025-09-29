# ThermDiffStructPhaseApplication 应用

## 概述

ThermDiffStructPhaseApplication是一个专为钢材热处理过程设计的多物理场耦合仿真应用模块。该模块实现了热-浓度-结构-相变四场强耦合求解，支持钢材常见的7种相组成和4种相变过程的完整建模与仿真。

### 核心特性
- **四场强耦合**: 在单个迭代步内实现热场、浓度场、结构场和相变场的强耦合计算
- **七相钢材模型**: 完整描述钢材的微观组织演化
- **统一相变管理**: 通过MaterialPhaseTransformationSystem统一管理所有相变过程
- **自包含模块**: 无需依赖其他应用模块，实现完整的热处理仿真功能
- **多尺度建模**: 从原子尺度的扩散到宏观尺度的变形的多尺度耦合


### 多物理场
- **热场**: 温度分布和热传导，考虑相变潜热
- **浓度场**: 碳浓度分布和扩散，考虑相间碳配分
- **结构场**: 应力应变和变形，考虑相变应变、TRIP和蠕变
- **相变场**: 相变动力学、相变塑性、蠕变、潜热、相变塑性、相变应变



## 支持模型
### 七相
- 奥氏体 (Austenite)
- 铁素体 (Ferrite) 
- 渗碳体 (Cementite)
- 珠光体 (Pearlite)
- 上贝氏体 (Upper Bainite)
- 下贝氏体 (Lower Bainite)
- 马氏体 (Martensite)
- 蠕变模型：1. Maxwell Linear 2. Maxwell power 

### 四种相变过程
- **奥氏体转变**: 加热过程中的奥氏体化
- **共析转变**: 奥氏体向珠光体的分解，偏离共析点还形成铁素体+渗碳体
- **贝氏体转变**: 中温区的贝氏体形成
- **马氏体转变**: 快冷过程中的马氏体形成

### 相变模型

#### 相变动力学模型
- **Avrami模型**: 适用于均匀成核和生长的相变
- **ISO模型**: 等温相变动力学
- **JMAK模型**: Johnson-Mehl-Avrami-Kolmogorov模型
- **KME模型**: Koistinen-Marburger方程（马氏体转变）
- **TTT模型**: 基于TTT曲线的相变动力学
- **UDef模型**: 用户自定义相变动力学模型

#### 相变塑性TRIP模型
- **Leblond模型**

- **Desalos模型**


#### 配分模型
- **Eutectoid模型**: 共析分解过程
  - 奥氏体 → 珠光体 + 铁素体+渗碳体
- **SubBn模型**: 贝氏体转变过程中
  - 奥氏体 → 上贝氏体和下贝氏体

#### 蠕变模型
- **Maxwell Linear模型**
  ```
- **Maxwell Power模型**
  ```


### 多相的材料属性

材料属性通过相分数加权平均计算：

#### 热物理属性
```
ρ_eff = Σ_{i=1}^{7} ξ_i × ρ_i        // 密度
k_eff = Σ_{i=1}^{7} ξ_i × k_i          // 热导率
Cp_eff = Σ_{i=1}^{7} ξ_i × Cp_i        // 比热
```

#### 力学属性
```
E_eff = Σ_{i=1}^{7} ξ_i × E_i          // 弹性模量
ν_eff = Σ_{i=1}^{7} ξ_i × ν_i        // 泊松比
σ_y_eff = Σ_{i=1}^{7} ξ_i × σ_{y,i}   // 屈服强度
α_eff = Σ_{i=1}^{7} ξ_i × α_i        // 热膨胀系数
```

#### 扩散属性
```
D_eff = Σ_{i=1}^{7} ξ_i × D_i          // 扩散系数
```

其中 ξ_i 为第i个的相分数

### 相变模型的参数

#### 相变动力学参数（每种相变12个参数）
- **奥氏体转变**: AUSTENIZATION_KINETICS_PARA1-12
  - 
- **共析转变**: EUTECTOID_KINETICS_PARA1-12
  - 
- **贝氏体转变**: BAINITIC_KINETICS_PARA1-12
  - 
- **马氏体转变**: MARTENSITIC_KINETICS_PARA1-12
  - 

#### TRIP参数（每种相变12个参数）
- **共析TRIP**: EUTECTOID_TRIP_PARA1-12
- **贝氏体TRIP**: BAINITIC_TRIP_PARA1-12
- **马氏体TRIP**: MARTENSITIC_TRIP_PARA1-12

#### 配分参数（每种相变12个参数）
- **共析配分**: EUTECTOID_PARTITION_PARA1-12
- **贝氏体配分**: BAINITIC_PARTITION_PARA1-12

#### 蠕变参数（每相类12个参数）
- **奥氏体蠕变**: AUSTENITE_CREEP_PARA1-12
- **铁素体蠕变**: FERRITE_CREEP_PARA1-12
- **渗碳体蠕变**: CEMENTITE_CREEP_PARA1-12
- **珠光体蠕变**: PEARLITE_CREEP_PARA1-12
- **上贝氏体蠕变**: UPPER_BAINITE_CREEP_PARA1-12
- **下贝氏体蠕变**: LOWER_BAINITE_CREEP_PARA1-12
- **马氏体蠕变**: MARTENSITE_CREEP_PARA1-12


## 代码结构

### 核心类架构

#### 1. MaterialPhaseTransformationSystem (相变系统核心类)
位置: `phase_transformation/material_phase_transformation_system.h/cpp`

**职责**: 相变过程的总控制类，管理所有相变计算

**关键特性**:
- 七相钢材管理: 奥氏体、铁素体、渗碳体、珠光体、上/下贝氏体、马氏体
- 四种相变过程协调: 奥氏体转变、共析转变、贝氏体转变、马氏体转变
- 子模型: 动力学模型、配分模型、TRIP模型、蠕变模型

**核心方法**:
```cpp
// 初始化相变系统
int Initialize(const Properties& rMaterialProperties,
               const GeometryType& rElementGeometry,
               const Vector& rShapeFunctionsValues);

// 获取有效材料属性
double GetDensity/GetConductivity/GetYoungModulus(...);

// 相变计算方法
int GetAusteniteMassFrcInc(...);
int GetEutectoidPhaseMassFrcInc(...);
int GetBainiteMassFrcInc(...);
int GetMartensiteMassFrcInc(...);

// 耦合效应计算
double GetLatent(...);                    // 潜热计算
double GetPhaseTransStrainInc(...);       // 相变应变
int GetTripInc(...);                      // TRIP应变
int GetCreepStrainInc(...);               // 蠕变应变
```

#### 2. ThermoChemoMechanicalPhaseChangeElement (多场耦合单元)
位置: `custom_elements/therm_diff_struct_phase_element.h/cpp`

**职责**: 四场耦合的有限元单元

**继承关系**: 继承自 `BaseSolidElement`

**变量结构**:
```cpp
struct CoupledElementVariables {
    // 时间积分参数
    double theta, dt, dt_inv;
    
    // 场变量值
    array_1d<double, TNumNodes> temperature_current/previous;
    array_1d<double, TNumNodes> concentration_current/previous;
    array_1d<double, TNumNodes> displacement_current/previous[TDim];
    
    // 相分数 (7相 × 节点数)
    array_1d<std::array<double, NumPhases>, TNumNodes> phase_fractions_current/previous;
    
    // 材料属性
    double effective_density, effective_conductivity, effective_young_modulus...;
    
    // 相变系统指针
    MaterialPhaseTransformationSystem::Pointer pPhaseTransSystem;
};
```

**核心方法**:
```cpp
// 基础有限元方法
void CalculateLocalSystem(...);
void EquationIdVector(...);
void GetDofList(...);

// 耦合计算方法
void InitializeCoupledElement(...);
void GetAllNodalValues(...);
void UpdateMaterialProperties(...);
PhaseTransformationResults CalculatePhaseTransformation(...);

// 各场贡献计算
void CalculateThermalContribution(...);
void CalculateDiffusionContribution(...);
void CalculateMechanicalContribution(...);
void AddCouplingTerms(...);
```

#### 3. FourFieldCoupledNewtonRaphsonStrategy (四场耦合求解策略)
位置: `custom_strategies/four_field_coupled_newton_raphson_strategy.h`

**职责**: 实现四场强耦合的迭代求解

**继承关系**: 继承自 `ImplicitSolvingStrategy`

**核心数据结构**:
```cpp
enum FieldIndex {
    TEMPERATURE_FIELD = 0,
    CONCENTRATION_FIELD = 1,
    STRUCTURE_FIELD = 2,
    PHASE_FIELD = 3
};

// 耦合系数矩阵
std::array<std::array<double, NUM_FIELDS>, NUM_FIELDS> mCouplingMatrix;

// 场收敛容差
std::array<double, NUM_FIELDS> mFieldTolerance;

// MaterialPhaseTransformationSystem指针
MaterialPhaseTransformationSystem::Pointer mpPhaseTransSystem;
```

**核心方法**:
```cpp
// 求解步验
void Initialize();
bool SolveSolutionStep();
void Clear();

// 耦合计算
void BuildCoupledSystem(...);
void ApplyCouplingTerms(...);
void UpdateCoupledFields(...);
bool CheckCoupledConvergence(...);

// 管理方法
void InitializeCouplingParameters();
void SetCouplingCoefficient(...);
```


## 主要内容

### 1. Element 单元文件 (`therm_diff_struct_phase_element.h/cpp`)





### 2. SolvingStrategy 求解策略文件 (`four_field_coupled_newton_raphson_strategy.h`)



## 文件结构

```
applications/ThermDiffStructPhaseApplication/
├── custom_elements/
│   ├── therm_diff_struct_phase_element.h        # 更新的单元头文件
│   └── therm_diff_struct_phase_element.cpp      # 更新的单元实现
├── custom_strategies/
│   └── four_field_coupled_newton_raphson_strategy.h  # 更新的策略
├── phase_transformation/
│   ├── material_phase_transformation_system.h   # 相变系统核心类
│   ├── kinetics/                                # 相变动力学模型
│   │   ├── avrami.h/cpp                        # Avrami模型
│   │   ├── iso.h/cpp                           # 等温相变模型
│   │   ├── jmak.h/cpp                          # JMAK模型
│   │   ├── kme.h/cpp                           # K-M方程模型
│   │   ├── ttt.h/cpp                           # TTT曲线模型
│   │   └── udef.h/cpp                          # 用户自定义模型
│   ├── partition/                              # 碳配分模型
│   │   ├── eutectoid.h/cpp                     # 共析分解配分
│   │   └── subbn.h/cpp                         # 贝氏体配分
│   ├── trip/                                   # TRIP模型
│   │   ├── leblond.h/cpp                       # Leblond模型
│   │   └── desalos.h/cpp                       # Desalos模型
│   ├── creep/                                  # 蠕变模型
│   │   ├── maxwell_linear.h/cpp                # 线性Maxwell模型
│   │   └── maxwell_power.h/cpp                 # 幂次Maxwell模型
│   └── phase.h/cpp                             # 单相基础类
├── therm_diff_struct_phase_application_variables.h/cpp  # 变量定义
├── test_integration.cpp                        # 集成测试
├── example_usage.py                           # Python使用示例
└── README.md                                   # 本文档
```

### 相变子模型组织架构

#### 相变动力学模型类层次
```
KineticsModel (基类)
├── Avrami          // 均匀成核成长模型
├── ISO             // 等温相变动力学
├── JMAK            // Johnson-Mehl-Avrami-Kolmogorov模型
├── KME             // Koistinen-Marburger方程
├── TTT             // TTT曲线动力学
└── UDef            // 用户自定义模型
```

#### TRIP模型类层次
```
TripModel (基类)
├── Leblond         // 适用于扩散型相变
└── Desalos         // 适用于无扩散型相变
```

#### 配分模型类层次
```
PartitionModel (基类)
├── Eutectoid       // 共析分解配分
└── SubBn           // 贝氏体转变配分
```

#### 蠕变模型类层次
```
CreepModel (基类)
├── MaxwellLinear   // 线性粘弹性蠕变
└── MaxwellPower    // 幂次粘弹性蠕变
```

### 计算流程

#### 1. 初始化阶段
```
1. 初始化MaterialPhaseTransformationSystem
2. 初始化各相的Phase对象
3. 初始化各相变的PhaseTransformationLaw对象
4. 读取材料参数和边界条件
5. 初始化节点变量（温度、浓度、位移、相分数）
```

#### 2. 时间步迭代阶段
```
for each time step:
    for each Newton-Raphson iteration:
        1. 计算相变过程
           │
           ├── GetAusteniteMassFrcInc()       // 奥氏体转变
           ├── GetEutectoidPhaseMassFrcInc()  // 共析转变
           ├── GetBainiteMassFrcInc()         // 贝氏体转变
           └── GetMartensiteMassFrcInc()      // 马氏体转变
           │
        2. 更新相分数和碳含量
           │
           ├── UpdateMassFraction()
           └── UpdateCarbonContent()
           │
        3. 计算有效材料属性
           │
           ├── GetDensity()          // 加权平均密度
           ├── GetConductivity()     // 加权平均热导率
           ├── GetYoungModulus()     // 加权平均弹性模量
           └── ...                   // 其他材料属性
           │
        4. 计算耦合效应
           │
           ├── GetLatent()               // 潜热释放
           ├── GetPhaseTransStrainInc()  // 相变应变
           ├── GetTripInc()              // TRIP应变
           └── GetCreepStrainInc()       // 蠕变应变
           │
        5. 求解各场方程
           │
           ├── CalculateThermalContribution()    // 热场
           ├── CalculateDiffusionContribution()  // 浓度场
           ├── CalculateMechanicalContribution() // 结构场
           └── AddCouplingTerms()                // 耦合项
           │
        6. 检查收敛
           │
           └── CheckCoupledConvergence()
```

### 数据流的转换

#### 输入数据
```
节点变量:
- 温度 T(x,t)
- 浓度 C(x,t) 
- 位移 u(x,t)
- 相分数 ξ_i(x,t)

材料参数:
- 各相材料属性
- 相变动力学参数
- TRIP参数
- 蠕变参数
```

#### 中间计算结果
```
相变计算:
- dξ_i/dt = 相变速率
- dC_i/dt = 碳配分速率

有效属性:
- ρ_eff, k_eff, E_eff, ...

耦合效应:
- Q_latent = 潜热源项
- ε_phase = 相变应变
- ε_trip = TRIP应变
- ε_creep = 蠕变应变
```

#### 输出结果
```
场变量更新:
- T^{n+1} = 新温度场
- C^{n+1} = 新浓度场
- u^{n+1} = 新位移场
- ξ_i^{n+1} = 新相分数场

应力应变结果:
- σ(x,t) = 应力场
- ε(x,t) = 应变场
```

## 使用指南

### 1. 材料参数设置

```python
# 在Python中设置材料属性
material_properties = Kratos.Properties(1)

# 开启多相和相变计算
material_properties[IS_MULTIPHASE_DEFINED] = True
material_properties[IS_PHASE_TRANSFORMATION_CALCULATED] = True

# 设置各相材料属性
material_properties[DENSITY_AUSTENITE] = 7900.0
material_properties[DENSITY_FERRITE] = 7870.0
# ... 其他相的属性

# 设置相变动力学参数
material_properties[AUSTENIZATION_KINETICS_PARA1] = 1000.0
# ... 其他12个参数
```

### 2. 模型创建

```python
# 创建主模型部件
main_model_part = Kratos.ModelPart("MainModelPart")

# 添加必要的变量
main_model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
main_model_part.AddNodalSolutionStepVariable(Kratos.CONCENTRATION)
main_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

# 添加相分数变量
main_model_part.AddNodalSolutionStepVariable(MASS_FRACTION_AUSTENITE)
# ... 其他6相的变量
```

### 3. 求解器设置

```cpp
// 在C++中创建四场耦合求解器
auto strategy = Kratos::make_unique<FourFieldCoupledNewtonRaphsonStrategy<...>>(
    model_part, scheme, linear_solver, convergence_criteria);

// 设置耦合系数
strategy->SetCouplingCoefficient(TEMPERATURE_FIELD, STRUCTURE_FIELD, 1e-5);
strategy->SetCouplingCoefficient(PHASE_FIELD, TEMPERATURE_FIELD, 1.0);
```

## 注意事项

1. **相分数守恒**: 确保Σ(相分数ᵢ) = 1.0
2. **参数调试**: 根据具体钢种调整相变动力学参数
3. **收敛性**: 相变计算可能影响收敛，需要合适的松弛参数
4. **时间步长**: 相变过程敏感，需要适当的时间步长
5. **初始条件**: 正确设置初始相组成和温度分布

## 技术优势

1. **统一架构**: MaterialPhaseTransformationSystem提供统一的相变管理
2. **模块化设计**: 各子模型独立实现，易于扩展和维护
3. **强耦合求解**: 四场在单个迭代步内强耦合计算
4. **物理一致性**: 保证相分数守恒和热力学一致性
5. **高度可配置**: 支持多种相变模型和参数组合
6. **性能优化**: 基于相分数的材料属性计算，避免重复计算