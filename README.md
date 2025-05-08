**TBlabENM version 0.0.1**
TBlabENM用于自动化批量完成生态位建模及后续分析。该包主要针对MaxEnt模型进行编写，能够多物种、多线程自动完成环境变量选择、模型参数优化、模型拟合与评价、适生区投影等一系列工作。此外，同样也能够对前期数据进行处理、后期模型结果分析，极大地提高了工作效率。
在TBlabENM包中，所有函数被分为三个模块：模型之前处理函数（Pre-Modelling Functions）、模型构建函数（Modelling Functions）和模型之后处理函数（Post-Modelling Functions）。模型之前处理函数包括环境变量统一、物种分布数据的检查和地理空间过滤；模型构建函数包括为多物种自动构建调优的MaxEnt模型及相关过程函数；模型之后处理函数包括MaxEnt模型变量贡献统计、为多物种批量进行适生区重分类和面积计算、物种丰富度计算、物种气候变化避难所计算及相关栅格绘图。
可在R中使用devtools::install_github("TBlab-why/TBlabENM")下载。
