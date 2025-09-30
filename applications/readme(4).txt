本目录下的内容用于修改或调整Kratos Application部分代码，注意当前所有内容以master-10.3版本为源本，不一定适用于其他版本

以下为文件列表：
   applications
  ├── ConstitutiveLawsSmallStrainApplication    						// 新增的小应变米塞斯塑性本构准则
  │   ├── CMakeLists.txt    								// CMake文件
  │   ├── constitutive_laws_small_strain_application.cpp  					// 本构模块源文件
  │   ├── constitutive_laws_small_strain_application.h  					// 本构模块头文件
  │   ├── ConstitutiveLawsSmallStrainApplication.json
  │   ├── ConstitutiveLawsSmallStrainApplication.py 
  │   ├── constitutive_laws_small_strain_application_variables.cpp  				// 本构模块变量定义源文件
  │   ├── constitutive_laws_small_strain_application_variables.h  				// 本构模块变量定义头文件
  │   ├── custom_constitutive   							// 本构文件夹
  │   │   │   └── elasticity
  │   │   │       ├── small_strain_j2_elasticity_3d.cpp   					    	// 3维 J2准则 弹性问题
  │   │   │       ├── small_strain_j2_elasticity_3d.h     				
  │   │   │       ├── small_strain_j2_elasticity_plane_strain.cpp  				// 平面应变 J2准则 弹性问题 派生自small_strain_j2_elasticity_3d
  │   │   │       └── small_strain_j2_elasticity_plane_strain.h  
  │   │   │       ├── small_strain_j2_elasticity_plane_stress.cpp  				// 平面应力 J2准则 弹性问题 派生自small_strain_j2_elasticity_3d 未完成
  │   │   │       └── small_strain_j2_elasticity_plane_stress.h      						
  │   │   │   └── plasticity
  │   │   │       ├── small_strain_j2_plasticity_3d.cpp   					// 3维 J2准则 弹塑性问题
  │   │   │       ├── small_strain_j2_plasticity_3d.h     				
  │   │   │       ├── small_strain_j2_plasticity_plane_strain.cpp  				// 平面应变 J2准则 弹塑性问题 派生自small_strain_j2_plasticity_3d
  │   │   │       └── small_strain_j2_plasticity_plane_strain.h  
  │   │   │       ├── small_strain_j2_plasticity_plane_stress.cpp  				// 平面应力 J2准则 弹塑性问题 派生自small_strain_j2_plasticity_3d 未完成
  │   │   │       └── small_strain_j2_plasticity_plane_stress.h      
  │   │   └── thermal   								// 用于小应变的热弹塑性本构
  │   │           ├── elasticity
  │   │           │   ├── small_strain_j2_thermal_elasticity_3d.cpp					// 3维 热弹性问题
  │   │           │   ├── small_strain_j2_thermal_elasticity_3d.h
  │   │           │   ├── small_strain_j2_thermal_elasticity_plane_strain.cpp				// 平面应变 热弹性问题 派生自thermal_elastic_isotropic_3d
  │   │           │   ├── small_strain_j2_thermal_elasticity_plane_strain.h
  │   │           │   ├── small_strain_j2_thermal_elasticity_plane_stress.cpp				// 平面应力 热弹性问题 派生自thermal_elastic_isotropic_3d 未完成
  │   │           │   └── small_strain_j2_thermal_elasticity_plane_stress.h
  │   │           └── plasticity
  │   │               ├── small_strain_j2_thermal_plasticity_3d.cpp					// 3维 J2准则 热弹塑性问题
  │   │               ├── small_strain_j2_thermal_plasticity_3d.h
  │   │               ├── small_strain_j2_thermal_plasticity_plane_strain.cpp			// 平面应变 J2准则 热弹塑性问题 派生自small_strain_j2_thermal_plasticity_3d
  │   │               ├── small_strain_j2_thermal_plasticity_plane_strain.h
  │   │               ├── small_strain_j2_thermal_plasticity_plane_stress.cpp			// 平面应力 J2准则 热弹塑性问题 派生自small_strain_j2_thermal_plasticity_3d 未完成
  │   │               └── small_strain_j2_thermal_plasticity_plane_stress.h
  │   ├── custom_processes
  │   ├── custom_python   								// 链接到python
  │   │   ├── add_custom_constitutive_laws_to_python.cpp
  │   │   ├── add_custom_constitutive_laws_to_python.h
  │   │   ├── add_custom_utilities_to_python.cpp
  │   │   ├── add_custom_utilities_to_python.h
  │   │   └── constitutive_laws_small_strain_python_application.cpp
  │   ├── custom_utilities								// 本构计算使用到的工具类
  │   │   ├── advanced_constitutive_law_small_strain_utilities.cpp				// 包含一些本构计算需要用到的工具函数
  │   │   ├── advanced_constitutive_law_small_strain_utilities.h
  │   │   ├── constitutive_law_utilities.cpp    					            // 本构计算工具类           
  │   │   ├── constitutive_law_utilities.h
  │   │   └── tangent_operator_calculator_utility.h						// 用于计算切线刚度矩阵的工具类
  │   ├── README.md
  │   └── tests   									// 测试文件夹
  ├── ThermDiffStructPhaseApplication   							// 新增的相变计算模块
  │   ├── CMakeLists.txt    								// CMake文件
  │   ├── custom_conditions   								// 自定义边界条件文件夹
  │   │   ├── axisym_line_load_condition.cpp						// 轴对称 线载荷边界
  │   │   ├── axisym_line_load_condition.h			
  │   │   ├── axisym_concentration_convection_condition.cpp					// 轴对称 浓度对流边界
  │   │   ├── axisym_concentration_convection_condition.h
  │   │   ├── axisym_thermal_cconvection_condition.cpp					// 轴对称 热面对流边界
  │   │   ├── axisym_thermal_convection_condition.h
  │   │   ├── axisym_point_load_condition.cpp						// 轴对称 点载荷边界
  │   │   ├── axisym_point_load_condition.h
  │   │   ├── base_load_condition.cpp							// 基础载荷边界
  │   │   ├── base_load_condition.h
  │   │   ├── concentration_convection_condition.cpp					// 浓度对流边界
  │   │   ├── concentration_convection_condition.h
  │   │   ├── concentration_flux_condition.cpp						// 物质流
  │   │   ├── concentration_flux_condition.h
  │   │   ├── line_load_condition.cpp							// 线载荷边界
  │   │   ├── line_load_condition.h
  │   │   ├── point_load_condition.cpp							// 点载荷边界
  │   │   ├── point_load_condition.h
  │   │   ├── small_displacement_line_load_condition.cpp					// 小变形线载荷边界
  │   │   ├── small_displacement_line_load_condition.h
  │   │   ├── small_displacement_surface_load_condition_3d.cpp				// 小变形面载荷边界
  │   │   ├── small_displacement_surface_load_condition_3d.h
  │   │   ├── surface_load_condition_3d.cpp						// 面载荷边界
  │   │   ├── surface_load_condition_3d.h
  │   │   ├── thermal_convection_condition.cpp						// 热对流边界
  │   │   ├── thermal_convection_condition.h
  │   │   ├── thermal_flux_condition.cpp							// 热流
  │   │   └── thermal_flux_condition.h
  │   ├── custom_elements   								// 自定义单元文件夹
  │   │   ├── base_solid_element.cpp							// 基础固体单元，基类
  │   │   ├── base_solid_element.h
  │   │   ├── thermo_diff_struct_phase_element.cpp						// 四场耦合单元 派生自base_solid_element
  │   │   └── thermo_diff_struct_phase_element.h
  │   ├── custom_python   								// python链接文件夹
  │   │   ├── add_custom_strategies_to_python.cpp
  │   │   ├── add_custom_strategies_to_python.h
  │   │   ├── add_custom_utilities_to_python.cpp
  │   │   ├── add_custom_utilities_to_python.h
  │   │   └── therm_diff_struct_phase_python_application.cpp
  │   ├── custom_strategies   								// SolvingStrategy文件夹
  │   │   └── four_field_coupled_newton_raphson_strategy.h
  │   ├── custom_utilities    								// 工具类
  │   │   ├── constitutive_law_utilities.cpp							// 本构计算需要的工具类
  │   │   ├── constitutive_law_utilities.h							// 本构计算需要的工具类
  │   │   ├── phase_trans_registration_utilities.h    						// 用于注册相变相关模型
  │   │   ├── read_phase_trans_materials_utility.cpp  						// 用于配置文件解析
  │   │   ├── read_phase_trans_materials_utility.h    						// 用于配置文件解析，继承自utilities/read_materials_utility.h
  │   │   ├── structural_mechanics_element_utilities.cpp					//一些结构力学单元计算的通用工具函数
  │   │   ├── structural_mechanics_element_utilities.h					//一些结构力学单元计算的通用工具函数
  │   │   ├── therm_diff_struct_phase_components.cpp  					// 用于相变相关模型的component定义
  │   │   └── therm_diff_struct_phase_components.h    					// 用于相变相关模型的component定义
  │   ├── phase_transformation  							// 相变模型类文件夹
  │   │   ├── creep     								// 蠕变模型文件夹
  │   │   │   ├── creep_model.cpp       							// 蠕变模型基类
  │   │   │   ├── creep_model.h       							// 蠕变模型基类
  │   │   │   ├── maxwell_linear.cpp  							// maxwell_linear蠕变模型类
  │   │   │   ├── maxwell_linear.h    							// maxwell_linear蠕变模型类
  │   │   │   ├── maxwell_power.cpp   							// maxwell_power蠕变模型类
  │   │   │   └── maxwell_power.h     							// maxwell_power蠕变模型类
  │   │   ├── kinetics  								// 相变动力学模型文件夹
  │   │   │   ├── avrami.cpp    								// Avrami方程
  │   │   │   ├── avrami.h
  │   │   │   ├── iso.cpp  								// iso等温
  │   │   │   ├── iso.h								
  │   │   │   ├── jmak.cpp  								// jmak方程
  │   │   │   ├── jmak.h
  │   │   │   ├── kinetics_model.cpp  							// 动力学模型基类
  │   │   │   ├── kinetics_model.h
  │   │   │   ├── kme.cpp  								// kme 方程
  │   │   │   ├── kme.h
  │   │   │   ├── ttt.cpp  								// ttt 曲线
  │   │   │   ├── ttt.h
  │   │   │   ├── udef.cpp  								// 用户定义
  │   │   │   └── udef.h
  │   │   ├── partition   								// 配分模型类文件夹
  │   │   │   ├── eutectoid.cpp								// 共析分解
  │   │   │   ├── eutectoid.h											
  │   │   │   ├── partition_model.cpp							// 配分模型基类
  │   │   │   ├── partition_model.h
  │   │   │   ├── subbn.cpp								// 贝氏体分解
  │   │   │   └── subbn.h
  │   │   ├── phase.cpp     								// 定义了相 Phase 类
  │   │   ├── phase.h       								// 定义了相 Phase 类
  │   │   ├── phase_transformation_law.cpp      						// 定义了相变类
  │   │   ├── phase_transformation_law.h        						// 定义了相变类
  │   │   ├── phase_transformation_system.cpp   						// 定义了相变过程控制类
  │   │   ├── phase_transformation_system.h     						// 定义了相变过程控制类
  │   │   ├── phase_transformation_define.h       						// 定义了相的种类枚举类和相变的种类枚举类
  │   │   └── trip    									// 相变塑性模型类文件夹
  │   │       ├── desalos.cpp       							// desalos 相变塑性模型类
  │   │       ├── desalos.h         								// desalos 相变塑性模型类
  │   │       ├── leblond.cpp       							// leblond 相变塑性模型类
  │   │       ├── leblond.h         								// leblond 相变塑性模型类
  │   │       ├── trip_model.cpp    							// 相变塑性模型基类
  │   │       └── trip_model.h      							// 相变塑性模型基类
  │   ├── python_scripts
  │   │   ├── python_solvers_wrapper_therm_diff_struct_phase.py   				// 四场耦合求解的求解器包装器，用于返回不同求解类型需要的求解器
  │   │   ├── therm_diff_struct_phase_analysis.py                 					// 四场耦合求解的流程控制类
  │   │   ├── therm_diff_struct_phase_coupled_solver.py            				// 四场耦合求解的求解器
  │   │   ├── apply_thermal_convection_process.py						// 用于解析配置文件中对于热对流边界的设置，还需要修改
  │   │   ├── apply_concentration_convection_process.py					// 用于解析配置文件中对于浓度传导边界的设置，还需要修改
  │   │   └── convergence_criteria_factory.py						// 用于解析配置文件中对于收敛条件的设置，还需要修改
  │   ├── README.md
  │   ├── tests     									// 测试文件夹
  │   ├── python_registry_lists.py
  │   ├── therm_diff_struct_phase_application.cpp   						// 四场耦合模块源文件
  │   ├── therm_diff_struct_phase_application.h     						// 四场耦合模块头文件
  │   ├── ThermDiffStructPhaseApplication.py        
  │   ├── therm_diff_struct_phase_application_variables.cpp   					// 四场耦合模块变量定义
  │   └── therm_diff_struct_phase_application_variables.h     					// 四场耦合模块变量定义
  └── TrilinosApplication     								// 删除了对于部分FSIApplication代码的依赖
      └── custom_python
          ├── add_custom_utilities_to_python.cpp
          └── add_trilinos_convergence_accelerators_to_python.cpp
