使用了类完成了 LB 模型计算 BCC 结构, 使用了两种方法: 
1. 构造了类 DGSolver, 该类包含了 DG 的一些信息, 然后在 main.cpp 中正常写程序, 用到 BCC 的信息时从 DGSolver 中获取
2. 构造了类 LBSolver, 该类继承了 DGSolver, 该类中包含了 LB 模型的计算, test.cpp 中实例化 LBSolver 类, 然后直接计算
