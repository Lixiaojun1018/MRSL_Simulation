
Rcpp::sourceCpp('STEP1.cpp')  # 加载STEP1的C++实现
Rcpp::sourceCpp('STEP2.cpp')  # 加载STEP2的C++实现
Rcpp::sourceCpp('DFS.cpp')    # 加载深度优先搜索（DFS）的C++实现
# 先有边际因果图再有边际因果图的D*D邻接矩阵：边际因果图是概念上的模型，它描述了变量之间的因果结构，而邻接矩阵是该模型的数学表示，用于进一步的统计分析和计算
# 定义MRSL函数，该函数接收汇总的遗传关联数据、标准误、阈值、调整方法等参数
# 参数说明：
# data_sum_beta 所有J个SNP对边际因果图中D个变量的汇总系数的J*D矩阵，例如，每个变量对每个SNP的线性回归系数。
# data_sum_se 对应于data_sum_beta的J*D标准误差矩阵。
# beta: J*D的指示矩阵，元素为0/1，例如，第i行第j列的元素为1表示第i个SNP是第j个节点的工具变量。
# cutoff 显著性置信水平
# adj_methods 选择调整变量的三种方法：1 移除collier；2 开放路径中的变量；3 最小充分集。
# use_eggers_step1：第一步中使用的多变量孟德尔随机化方法：1 MVMR-Egger；0 MVMR-IVW
# use_eggers_step2 第二步中使用的多变量孟德尔随机化方法：1 MVMR-Egger；0 MVMR-IVW。
# vary mvmr adj=0表示直接对mvmr中的调整变量进行调整；vary mvmr adj=1表示用adj_methods从选定的变量改变mvmr中调整变量的数量。
# R:预测的表型间的效应矩阵
MRSL <- function(data_sum_beta, data_sum_se, beta, cutoff,
                 adj_methods, use_eggers_step1, use_eggers_step2, vary_mvmr_adj, R){
  
  # 使用STEP1函数处理汇总数据
  res1 <- STEP1(data_sum_beta, data_sum_se, beta, use_eggers_step1)
  # 初始化一个列表来存储邻接矩阵
  #邻接矩阵
  amatt0 <- list() 
  amatt0[[1]] <- res1$amatt  # 将STEP1的结果存储为邻接矩阵的第一个元素
  
  # 获取邻接矩阵
  
  pp <- amatt0[[1]]
  pp<-M_amatt
  # 使用DFS算法对邻接矩阵进行拓扑排序，并获取排序后的节点索引
  v1 <- 10 - c(DFS(amattt=pp, n_nodes=ncol(pp)))  # ncol(pp)是邻接矩阵的列数，即节点数
  
  # 创建一个序列，表示节点的原始顺序
  v2 <- 1:ncol(data_sum_beta)
  
  # 计算Spearman Footrule距离，衡量两个排列之间的差异
  SpearF <- sum(sapply(seq_along(v1), function(i) abs(i - (which(v2 == v1[i])))))
  
  # 计算Kendall等级相关系数，评估两个排列之间的相关性
  kenDcor <- cor.test(v1, v2, use = "pairwise", method="kendall") 
  # 提取Kendall等级相关系数的结果
  kenDcor_res <- c(kenDcor$statistic, kenDcor$p.value, kenDcor$estimate)
  
  # 根据DFS排序重新排列邻接矩阵
  amatt0[[1]] <- amatt0[[1]][v1, v1]
  
  # 根据DFS排序重新排列相关系数矩阵R
  R1 <- R[v1, v1]
  
  # 初始化一个布尔变量，用于控制迭代过程
  just <- TRUE
  i <- 1
  while(just){
    cat("iteration--", i)  # 打印当前迭代次数
    i <- i + 1  # 迭代次数加1
    # 使用STEP2函数处理邻接矩阵，可能进行进一步的统计测试或调整
    amatt055 <- STEP2(amatt0[[i-1]], data_sum_beta, data_sum_se, beta,
                      adj_methods, use_eggers_step2, vary_mvmr_adj)
    amatt0[[i]] <- amatt055$amatt  # 更新邻接矩阵列表
    # 如果新的邻接矩阵与上一次迭代的结果相同，则停止迭代
    just <- !all(amatt0[[i]] == amatt0[[i-1]])
  }
  
  # 返回最终结果，包括最终的邻接矩阵、迭代次数、Spearman Footrule距离、Kendall等级相关系数结果和相关系数矩阵
  return(list(amatt=amatt0[[i]],
              iteration=i,
              SpearF=SpearF,
              kenDcor_res=kenDcor_res,
              R1=R1))
}

