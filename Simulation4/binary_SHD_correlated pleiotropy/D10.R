
library(openxlsx)
# D_beta 和 D_beta_M 分别是从两个表格中读取的矩阵数据(描述了DAG的结构或因果关系的假设矩阵。)
D_beta <- read.xlsx('graph_matrix.xlsx',
                    sheet='asia',rowNames = T)
D_beta_M <- read.xlsx('graph_matrix.xlsx',
                      sheet='asia_step1',rowNames = T)
D_beta_M[is.na(D_beta_M)] <- 0 #D_beta_M 中的缺失值被替换为 0。

############################## edge effect 0.25-0.5 ######################
############################## weak IVs #####################
source("Comp.R")

cutoff=5e-02 #设置参数
NN=50  #外层模拟的次数。
mc=50  #并行任务的核心数。
N=10000  
gl=0.05  #参数范围，用于筛选/剪枝。
gu=0.2
dlB=round(log(1.5),2)
duB=round(log(2),2)
each_g=50

prob0=0 #prob0当前被固定为0，似乎用作保留参数。
# prob_corp0 是一组不同的概率值，表示不同的情景（如因果关系的强度）。
prob_corp0 = c(0,0.1,0.3,0.5,0.8)  # 通过设置不同组合，实现敏感性分析

oo=1
result_all <- NULL


for(ku in 1:length(prob_corp0)){
  for(kb in 1:length(prob0)){
    prob_corp <- prob_corp0[ku]

    # 并行计算 
    registerDoMC(mc)
    tt1 <- foreach(x=1:NN,
                   .packages = c("bnlearn","MRPC","pcalg",
                                 "dagitty","Rcpp","RcppArmadillo",
                                 "MendelianRandomization","MRPRESSO",
                                 "R.utils","mr.raps","MRMix","MRCD"),
                   .export = c("STEP1.cpp", "STEP2.cpp","DFS.cpp")) %dopar% {
                     #withTimeout({
                     Comp(x,N,D_beta,each_g,gl,gu,dlB,duB,prob_corp,D_beta_M)
                     #},timeout = 960, onTimeout = "silent")
                   }
   # tt1 <- Comp(x=1,N,D_beta,each_g,gl,gu,dlB,duB,prob_corp,D_beta_M)

    # tt1 <- list()
    # for(yt in 1:NN){
    #   tt1[[yt]] <- Comp(x=yt,N,D,each_g,gl,gu,prob,dlB,duB)
    # }
    
    filename0 <- paste0("Asia_","_prob_corp_",prob_corp,"_addorder.Rdata")
    save(tt1,file=filename0)
    
    rm(tt1)
    
  }
  
}



#rm(list=ls())

