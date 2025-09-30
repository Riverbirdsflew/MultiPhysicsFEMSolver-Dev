# 相变效应+本构模型

## 1. 基本原理

在考虑相变效应的本构计算中，总应变需要分解为多个分量：

```
总应变 = 弹性应变 + 热应变 + 相变应变 + 相变塑性应变 + 蠕变应变
```

因此，在计算本构关系时，应该：
1. 计算各种非弹性应变分量
2. 从总应变中减去这些非弹性应变，得到纯弹性应变
3. 基于纯弹性应变计算应力

## 2. PhaseTransformationSystem提供的函数

PhaseTransformationSystem类提供了以下与应变相关的函数：

### 2.1 应变增量计算函数
- `GetThermalStrainInc()`: 返回热应变增量
- `GetPhaseTransStrainInc()`: 返回相变应变增量
- `GetCreepStrainInc()`: 返回蠕变应变增量
- `GetTripInc()`: 返回相变塑性应变(TRIP)增量

### 2.2 材料属性获取函数
- `GetElasticityModulus()`: 返回弹性模量
- `GetPoisson()`: 返回泊松比

## 3. 实现方法

### 3.1 非机械应变计算

需要在对应的单元中添加几个函数来计算非弹性应变：

```cpp
for(size_t point = 0; point < integration_points.size(); ++point) {
            // ========== 第一步：计算总应变 ==========
            Vector total_strain = CalculateTotalStrainAtIntegrationPoint(point);
            
            // ========== 第二步：计算非机械应变 ==========
            //Vector thermal_strain = CalculateThermalStrainIncrement(point, rCurrentProcessInfo);
            Vector transformation_strain = CalculateTransformationStrain(point, rCurrentProcessInfo);
            Vector trip_strain = CalculateTRIPStrain(point, rCurrentProcessInfo);
            Vector creep_strain = CalculateCreepStrainIncrement(point, rCurrentProcessInfo);
            
            // 更新累积应变
            //mThermalStrainVector[point] += thermal_strain;
            mTransformationStrainVector[point] += transformation_strain;
            mTRIPStrainVector[point] += trip_strain;
            mCreepStrainVector[point] += creep_strain;
            
            // ========== 第三步：计算力学应变 ==========
            Vector mechanical_strain = total_strain
                                     //- mThermalStrainVector[point]
                                     - mTransformationStrainVector[point]
                                     - mTRIPStrainVector[point]
                                     - mCreepStrainVector[point];
            
            // ========== 第四步：本构计算（纯弹塑性） ==========
            ConstitutiveLaw::Parameters cl_params = PrepareConstitutiveParameters(
                mechanical_strain, point, rCurrentProcessInfo
            );
            
            mConstitutiveLawVector[point]->CalculateMaterialResponseCauchy(cl_params);
            
            // ========== 第五步：组装到整体 ==========
            const Vector& stress = cl_params.GetStressVector();
            AssembleInternalForces(stress, point, rRightHandSideVector);
        }
```

### 3.2 减去非机械应变后的应力计算

```cpp
void SmallStrainJ2Elasticity3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues)
{
    // 计算弹性应变
    ConstitutiveLaw::StrainVectorType elastic_strain_vector;
    CalculateElasticStrain(rStrainVector, elastic_strain_vector, rValues);
    
    // 创建PhaseTransformationSystem实例获取材料属性
    PhaseTransformationSystem phaseTransformationSystem;
    
    // 获取节点ID（需要根据实际情况确定）
    SizeType node_id = 0;
    
    // 获取材料属性
    const double E = phaseTransformationSystem.GetElasticityModulus(
        rValues.GetMaterialProperties(),
        rValues.GetElementGeometry(),
        rValues.GetShapeFunctionsValues(),
        rValues.GetProcessInfo(),
        node_id
    );
    
    const double NU = phaseTransformationSystem.GetPoisson(
        rValues.GetMaterialProperties(),
        rValues.GetElementGeometry(),
        rValues.GetShapeFunctionsValues(),
        rValues.GetProcessInfo(),
        node_id
    );

    // 基于弹性应变计算应力
    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(
        rStressVector, elastic_strain_vector, E, NU);
}
```

## 4. 

### 4.1 节点ID
在实际实现中，需要根据具体的单元和高斯点来确定node_id。通常可以通过以下方式获取：
- 通过rValues.GetElementGeometry()获取单元几何信息
- 通过rValues.GetShapeFunctionsValues()获取形函数值
- 通过形函数值确定高斯点对应的位置

### 4.2 相分数和相分数增量的获取
相分数(mf)和相分数增量(mfinc)需要从求解器的状态变量中获取：
- 通过r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_*)获取相分数
- 通过r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_*, 1)获取上一步的相分数
