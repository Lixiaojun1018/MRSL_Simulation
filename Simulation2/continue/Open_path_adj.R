# 定义 Open_path_adj 函数
Open_path_adj <- function(MM, start0, end0, adj) {
  # 获取 MM 矩阵中非零元素的位置
  locname <- which(MM != 0, arr.ind = TRUE)
  
  # 初始化 DAG 语句
  dagstatement <- "dag {"
  
  # 遍历非零元素的位置，构建 DAG 语句
  for (i in 1:nrow(locname)) {
    dagstatement <- paste0(dagstatement, "D", locname[i, 1], " -> ", "D", locname[i, 2], " ; ")
  }
  
  # 去掉最后一个多余的分号
  dagstatement <- substr(dagstatement, 1, nchar(dagstatement) - 3)
  
  # 完成 DAG 语句
  dagstatement <- paste0(dagstatement, "}")
  
  # 使用 dagitty 包创建 DAG 图
  g <- dagitty(dagstatement)
  
  # 将起始节点和结束节点转换为 "D" 开头的字符串
  start0 <- paste0("D", start0)
  end0 <- paste0("D", end0)
  
  # 根据 adj 参数的不同值进行不同的处理
  if (adj == 2) {
    # 初始化节点调整集
    node_adj <- NULL
    
    # 获取所有路径
    # pp的长度即为路径的条数
    pp <- paths(g, start0, end0, directed = FALSE)  # 列出所有路径
    # pp长度不等于0
    if (length(pp) != 0) {
      # 获取开放路径
      openpath <- pp$paths[pp$open]
      # 将路径拆分为节点列表
      nodes_all <- strsplit(openpath, " ")
      # 节点列表不为0
      if (length(nodes_all) != 0) {
        # 遍历所有路径
        for (oi in 1:length(nodes_all)) {
          # 第oi条路径
          node_once <- nodes_all[[oi]]
          # 如果这条路的节点>3
          if (length(node_once) > 3) {
            # 提取中间节点
            node_loc <- seq(3, length(node_once) - 2, 2)
            node_adj <- c(node_adj, node_once[node_loc])
          }
        }
      }
    }
  } else if (adj == 3) {
    # 获取调整集
    # g:DAG
    loc_adj0 <- adjustmentSets(g, start0, end0)
    if (length(loc_adj0) == 0) {
      node_adj <- NULL
    } else {
      node_adj <- loc_adj0[[1]]
    }
    
    # 获取所有有向路径
    pp2 <- paths(g, start0, end0, directed = TRUE)$paths
    
    if (length(pp2) != 0) {
      # 将路径拆分为节点列表
      nodes_all <- strsplit(pp2, " ")
      
      for (oi in 1:length(nodes_all)) {
        node_once <- nodes_all[[oi]]
        if (length(node_once) > 3) {
          # 提取中间节点
          node_loc <- seq(3, length(node_once) - 2, 2)
          node_adj <- c(node_adj, node_once[node_loc])
        }
      }
    }
  }
  
  # 处理节点调整集
  if (length(node_adj) != 0) {
    # 去重
    node_adj <- unique(node_adj)
    # 去掉 "D" 前缀并转换为数值
    node_adj <- as.numeric(substr(node_adj, 2, nchar(node_adj)))
    # 调整索引，使其从 0 开始
    node_adj <- node_adj - 1
  }
  
  # 确保返回值为数值类型
  node_adj <- as.numeric(node_adj)
  
  # 返回节点调整集
  return(node_adj)
}