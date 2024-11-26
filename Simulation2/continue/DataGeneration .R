# N: 个体数量。10000
# D: 表型数量。
# each_g: 每个表型对应的工具变量数量。
# gl: 工具变量效应的下限。
# gu: 工具变量效应的上限。
# prob: 随机图中边存在的概率。
# dlB: 表型之间效应的下限。
# duB: 表型之间效应的上限


DataGeneration <- function(N,D,each_g,gl,gu,prob,dlB,duB){
  #总工具变量数量：D*each_g：每个表型对应的工具变量数量。each_g：额外的工具变量数量，这些工具变量影响所有表型。
  total_g <- D*each_g + each_g     
  #生成工具变量矩阵 G
  G <- NULL      
  for(j in 1:total_g){
    G <- cbind(G,rbinom(N,2,0.3))   #生成N个个体的工具变量，每个工具变量从二项分布B(2,0.3)生成。
  }
  #生成工具变量对表型的效应矩阵 G_beta
  # G-Y coefficients
  G_beta <- matrix(0,ncol=D,nrow=total_g)
  #生成效应：
  for(d in 1:D){
    ww <- (1 + each_g*(d-1)):(each_g*d)    #确定每个表型对应的工具变量索引。
    G_beta[ww,d] <- runif(length(ww),gl,gu)  #生成均匀分布的效应值。
  }
  G_beta[(D*each_g+1):total_g,] <- runif(each_g*D,gl,gu)   #生成额外工具变量对所有表型的效应值。
  colnames(G_beta) <- paste0("D",1:D)
  
  #生成表型之间的效应矩阵 D_beta
  # D-Y coefficients
  D_beta <- matrix(0,nrow=D,ncol=D)
  nmbEdges <- 0L
  for (i in seq_len(D - 2)) {
    listSize <- rbinom(1, D - i, prob)   #确定每个表型的父节点数量。
    nmbEdges <- nmbEdges + listSize      
    edgeList <- sample(seq(i + 1, D), size = listSize)   #随机选择父节点。
    weightList <- runif(length(edgeList), min = dlB, max = duB)   #生成均匀分布的效应值。
    
    D_beta[i,edgeList] <- weightList
  }
  # 生成未测量的混杂因素 U
  U <- NULL
  for(ed in 1:D){
    U <- cbind(U,rnorm(N,0,1))
  }
  # 生成表型数据 data_ind
  data_ind <- (G %*% G_beta + U) %*% solve(diag(D)-D_beta)   #生成表型数据，考虑工具变量效应、混杂因素和表型之间的相互作用。
  #添加残差项。
  for(tty in 1:ncol(data_ind)){
    data_ind[,tty] <- data_ind[,tty]+rnorm(N,0,1)
  }
  
  colnames(data_ind) <- paste0("D",1:D)
  #生成汇总数据 data_sum
  #### summary data
  data_sum <- vector("list",D)
  for(k in 1:total_g){
    for(d1 in 1:D){
      fit <- lm(data_ind[,d1]~G[,k])   #对每个表型和每个工具变量进行线性回归。
      data_sum[[d1]] <- rbind(data_sum[[d1]],c(beta=summary(fit)$coef[2,1],    
                                               se=summary(fit)$coef[2,2],
                                               pval=summary(fit)$coef[2,4]))
    }
  }
  names(data_sum) <- paste0("D",1:D)
  #生成得分 score
  data_sum_beta <- NULL
  for(f in 1:D){
    data_sum_beta <- cbind(data_sum_beta,data_sum[[f]][,1])
    data_sum[[f]] <- as.data.frame(data_sum[[f]])
  }
  
  score <- NULL
  for(d2 in 1:D){
    ed <- (1 + each_g*(d2-1)):(each_g*d2)
    score <- cbind(score,G[,ed] %*% data_sum_beta[ed,d2])
  }
  colnames(score) <- paste0("score_D",1:D)
  
  colnames(G) <- paste0("SNP",1:ncol(G))
  
  return(list(data_ind=cbind(data_ind,G),
              data_sum=data_sum,
              score=score,
              nmbEdges=nmbEdges,
              D_beta=D_beta,
              G_beta=G_beta))
}
DataGeneration(20,4,2,0.1,0.3,0.2,0.2,0.3)
