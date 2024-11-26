#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <R.h>
#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;
using namespace arma;
using namespace std;


// combination：这个函数用于生成从n个元素中取k个元素的组合
// 定义生成下一个组合的函数
// comb:comb 是一个数组，用于存储当前生成的组合。
// 在组合生成算法中，comb 数组的每个元素表示从 0 到 n-1 中选择的一个元素的位置。具体来说，comb 数组的长度为 k，表示组合的长度，每个元素的值表示选择的元素在原集合中的索引
// comb 数组的初始值是从 0 到 k-1 的连续整数。例如，对于 k = 3，初始值为 [0, 1, 2]。
bool next_comb(int* comb, const int n, const int k) {
  // 初始化变量 i 为 k - 1，这里的 k 通常表示组合的长度。
  // i 用于从组合数组的最后一个元素开始向前遍历。
  // eg：K=2，则组合最后一个元素就是comb[1],第一个元素就是comb[0]
  int i = k - 1;
  // e 表示组合数组中每个元素的最大可能值（相对于其初始位置）。
  // eg:n=1
  const int e = n - k;
  do
    comb[i]++; // 当前索引位置的值加 1，为了尝试生成下一个组合。
  // i--：如果当前 comb[i] 超出边界，i 会递减，直到找到一个可以增加的值或 i 变为负数（即所有位置都已尝试过）
  while (comb[i] > e + i && i--); // 如果超出边界，回溯到前一个位置并增加其值
  
  if (comb[0] > e) // 如果第一个元素超出边界，表示没有更多组合;eg:
    return 0;
  
  while (++i < k) // 重置后面的元素
    comb[i] = comb[i - 1] + 1; 
  
  return 1; // 成功生成下一个组合
}
// cho[0]：组合第一个元素
int Combination(int* cho, int nl, int n, int k) {
  if (n < k || k <= 0) // 如果组合条件不满足，直接返回空列表
    return 0;
  
  int* comb = new int[k]; // 分配组合数组
  for (int i = 0; i < k; i++)
    comb[i] = i; // 初始化第一个组合 [0, 1, ..., k-1]
  
  do {
    for (int i = 0; i < k; ++i) {
      cho[nl] = comb[i] + 1; // 存储组合，加1是为了从1开始计数
      nl++;
    }
  } while (next_comb(comb, n, k)); // 生成下一个组合
  
  delete[] comb; // 释放组合数组
  return nl; // 返回组合的总数
}



// choose IVs
// 定义选择工具变量 (IVs) 的函数
// 定义 Cho_IV 函数，用于找出只对特定列 nodex 有影响的行的索引
void Cho_IV(const int nodex, const arma::mat& bet, arma::uvec& l1){
  // 初始化 ttsnp 向量，用于标记哪些行只有一个非零元素
  arma::vec ttsnp(bet.n_rows);
  
  // 遍历 bet 矩阵的每一行
  for(int rsnp=0; rsnp<bet.n_rows; rsnp++){
    // 检查当前行中非零元素的数量是否为 1
    if(sum(bet.row(rsnp) != 0) == 1){
      // 如果是，标记该行为只有一个非零元素
      ttsnp[rsnp] = 1;
    } else {
      // 如果不是，标记该行为不满足条件
      ttsnp[rsnp] = 0;
    }
  }
  
  // 使用 ttsnp 向量中的标记，提取出只有一个非零元素的行
  arma::mat bet0 = bet.rows(find(ttsnp == 1));
  
  // 在提取出的行中，找出只影响 nodex 列的行
  l1 = find(bet0.col(nodex) != 0);

  
 // cout << "l1==" << l1 << endl;
  
} 


// 定义 fastLm 函数，用于执行线性回归
// 参数:
// X1 - 样本（子节点）
// y - 响应变量（父节点）
// intercept - 布尔值，表示是否添加截距项
Rcpp::List fastLm(const arma::mat& X1, const arma::colvec& y, bool intercept) {
  
  int n = X1.n_rows; // 获取样本数量:获取X1的行数
  
  mat X; // 初始化样本矩阵（子节点）
  if(intercept) {
    // 如果需要截距项，则在X1前添加一列全1
    X = join_rows(ones(n, 1), X1);
  } else {
    // 如果不需要截距项，则直接使用X1
    X = X1;
  }
  
  int k = X.n_cols; // 获取设计矩阵的列数（包括截距项）
  
  // 使用伪逆计算系数
  // coef = arma::spsolve(X, y); // 另一种计算系数的方法，使用稀疏矩阵求解
  arma::colvec coef = arma::pinv(X) * y; // 使用伪逆计算系数
  arma::colvec res = y - X * coef; // 计算残差
  
  // 计算系数的标准误差
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0) / (n - k); // 计算残差平方和
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X) * X))); // 计算标准误差
  
  // 计算t值
  arma::colvec tval = coef / std_err;
  
  // 将t值转换为p值
  NumericVector xx = as<NumericVector>(wrap(-abs(tval))); // 计算-tval的绝对值
  double df = n - k; // 计算自由度
  
  NumericVector pval = Rcpp::pt(xx, df) * 2; // 计算双尾p值
  
  // 返回线性回归的结果
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coef, // 系数
    Rcpp::Named("stderr") = std_err,    // 标准误差
    Rcpp::Named("tvalue") = tval,       // t值
    Rcpp::Named("pvalue") = pval,       // p值
    Rcpp::Named("df") = n - k           // 自由度
  );
}

// 定义 Judge 函数，用于判断是否存在边缘效应
// betaX：子节点的对应的工具变量矩阵
// betaY：父节点的对应的工具变量矩阵
// seY：父节点对选择工具变量回归对应的标准误
int Judge(const arma::mat& betaX, const arma::vec& betaY, const arma::vec& seY,
          const int use_egger) {
  
  int edge = 0; // 初始化边缘效应标志为0（无边缘效应）
  
  // 构建标准误差的逆矩阵
  mat sigma = diagmat(1/seY);
  
  // 对betaX和betaY进行标准化处理
  mat betaX1 = sigma * betaX;
  mat betaY1 = sigma * betaY;
  
  if(use_egger == 1) {
    // 如果使用Egger回归测试，首先进行带截距的线性回归
    List lm1 = fastLm(betaX1, betaY1, true);
    arma::vec pval1 = as<arma::vec>(lm1["pvalue"]); // 获取p值
    
    // 检查截距项的p值是否大于0.05
    if(pval1[0] > 0.05) {
      // 如果截距项的p值大于0.05，进行不带截距的线性回归
      List lm2 = fastLm(betaX1, betaY1, false);
      arma::vec pval2 = as<arma::vec>(lm2["pvalue"]); // 获取p值
      // 检查第一个系数的p值是否小于0.05
      if(pval2[0] < 0.05)
        edge = 1; // 如果是，则标记为存在边缘效应
    } else if(pval1[1] < 0.05) {
      // 如果第一个系数的p值小于0.05，则标记为存在边缘效应
      edge = 1;
    }
    
  } else if(use_egger == 0) {
    // 如果不使用Egger回归测试，直接进行不带截距的线性回归
    List lm2 = fastLm(betaX1, betaY1, false);
    arma::vec pval2 = as<arma::vec>(lm2["pvalue"]); // 获取p值
    // 检查第一个系数的p值是否小于0.05
    if(pval2[0] < 0.05)
      edge = 1; // 如果是，则标记为存在边缘效应
  }
  
  // 返回边缘效应标志
  return edge;
}


// [[Rcpp::export]]

SEXP STEP1(SEXP data_sum_beta,SEXP data_sum_se,SEXP beta1,  
            SEXP use_eggers_step1){
  // data_sum_beta 所有J个SNP对边际因果图中D个变量的汇总系数的J*D矩阵，例如，每个变量对每个SNP的线性回归系数。
  // data_sum_se 对应于data_sum_beta的J*D标准误差矩阵。
 //  beta: J*D的指示矩阵，元素为0/1，例如，第i行第j列的元素为1表示第i个SNP是第j个节点的工具变量。
// use_eggers_step1：第一步中使用的多变量孟德尔随机化方法：1 MVMR-Egger；0 MVMR-IVW
  const arma::mat sum_beta = as<arma::mat>(data_sum_beta);
  const arma::mat sum_se = as<arma::mat>(data_sum_se);
  const arma::mat beta = as<arma::mat>(beta1);
  const int use_egger_step1=as<int>(use_eggers_step1);
  
  int D1 = sum_beta.n_cols, D2 = sum_se.n_cols;
  if (D1 != D2){
    perror("The dimensions of beta and se are not matched");
  } 
  
  mat amatt;
  amatt = zeros<mat>(D1,D1);// 空矩阵
  
  int* cho = new int[100000000];
  int nl=0;
  nl=Combination(cho,nl,D1,2);
  
  // step 1：第一个方向
  cout << "step1 all for the first direction " << nl << "..." << endl;
  for(int i=0; i<nl/2; i++){
    // int i=0;
    uvec l1;
    //选择与子节点有效应的工具变量
    Cho_IV(cho[2*i],beta,l1);
    
    //cout << cho[2*i] << " " << cho[2*i+1] << endl;
    mat betaX1 = sum_beta.col(cho[2*i]);
    betaX1 = betaX1.elem(l1);
    vec betaY1 = sum_beta.col(cho[2*i + 1]);
    betaY1 = betaY1.elem(l1);
    vec seY1 = sum_se.col(cho[2*i + 1]);
    seY1 = seY1.elem(l1);
    // # 判断是否存在有效的边
    amatt(cho[2*i],cho[2*i + 1]) = Judge(betaX1,betaY1,seY1,use_egger_step1);
    
  } 
  // 第二个方向
  cout << "step1 all for the second direction " << nl << "..." << endl;
  for(int i=0; i<nl/2; i++){
    // int i=0;
    uvec l2;
    Cho_IV(cho[2*i+1],beta,l2);
    
    //cout << cho[2*i+1] << " " << cho[2*i] << endl;
    mat betaX2 = sum_beta.col(cho[2*i+1]);
    betaX2 = betaX2.elem(l2);
    vec betaY2 = sum_beta.col(cho[2*i]);
    betaY2 = betaY2.elem(l2);
    vec seY2 = sum_se.col(cho[2*i]);
    seY2 = seY2.elem(l2);
    amatt(cho[2*i+1],cho[2*i]) = Judge(betaX2,betaY2,seY2,use_egger_step1);
    
  }
  
  cout << "step1 over" << endl;
  //cout << "amatt " << amatt << endl;
  
  //cout << "memory " << GetSysMemInfo() << endl;
  

  delete[] cho;
  
  return List::create(Rcpp::Named("amatt") = amatt);
}
