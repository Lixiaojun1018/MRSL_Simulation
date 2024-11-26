# 这是一个用于计算两个邻接矩阵m1和m2的结构汉明距离（SHD, Structural Hamming Distance）的R语言函数。
# SHD通常用于比较网络结构，例如贝叶斯网络或因果图的结构相似性。
# SHD是衡量两个有向图之间差异的指标，通常##用于比较推断的因果图和真实因果图之间的相似性##

#这个函数的核心是通过比较邻接矩阵，计算删除、添加边以及修改方向所需的最小操作数，从而得到结构汉明距离SHD。
#这在因果推断和图结构学习中非常有用，用于评估学习算法与真实结构之间的差异。

SHD  <- function(m1,m2){    #m1,m2为矩阵形式的DAG
#初始化
  shd <- 0
  ## Remove superfluous edges from g1 
  s1 <- m1 + t(m1)    #s1 s2: 将两个有向图邻接矩阵转换为无向图形式
  s2 <- m2 + t(m2)    #m1+t(m1)是矩阵与其转置的和，表示无向图的边。 
  s1[s1 == 2] <- 1    #如果某条边在两个方向上都存在，即m[i→j]=m[j→i]=1,结果会是2，此时这些值将重置为1
  s2[s2 == 2] <- 1
#删除多余边  
  ds <- s1 - s2    #ds: 计算两个无向图之间的差异。
                   #ds[i→j]=1：图1中存在的边，但图2中不存在
                   #ds[i→j]=1：两图在该位置的边一致
                   #ds[i→j]=1：图2中存在的边，但图1中不存在
  ind <- which(ds > 0)    #ind: 找到图1中多余的边的位置
  m1[ind] <- 0    #删除多余边
  shd <- shd + length(ind)/2    #更新shd，length(ind)/2为删除边的数量
#添加缺失边
  ## Add missing edges to g1
  ind <- which(ds < 0)    #ind: 找到图1中缺失的边的位置
  m1[ind] <- m2[ind]    #将这些边添加到图1中。
  shd <- shd + length(ind)/2    #更新shd，length(ind)/2为添加边的数量
#比较方向差异
  ## Compare Orientation
  d <- abs(m1-m2)
  ## return
  SHD <- shd + sum((d + t(d)) > 0)/2    #sum((d + t(d)) > 0)/2为统计方向不一致的边数
                                        #d + t(d)是无向表示的差异矩阵，将方向差异对称化。
                                        #除以2是因为对称矩阵包含重复的边。  
  return(SHD)
}

#eg.
m1 <- matrix(c(0, 1, 0,
               0, 0, 1,
               0, 0, 0), nrow=3, byrow=TRUE)

m2 <- matrix(c(0, 0, 1,
               1, 0, 0,
               0, 1, 0), nrow=3, byrow=TRUE)

SHD(m1, m2)

