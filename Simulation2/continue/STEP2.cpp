#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <R.h>
#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;
using namespace arma;
using namespace std;
// 定义 Open_Path_Adj 函数，用于寻找两个节点之间的开放路径
// MM:边际因果图邻接矩阵
// start：父节点
// end：子节点
// adj：调整值。
arma::vec Open_Path_Adj(const mat& MM, const int start, const int end,
                        const int adj){
  // Function 类型通常来自 Rcpp 库，用于在 C++ 中调用 R 函数。
  // Open_path_adj：外部函数R函数
  Function f("Open_path_adj");
  NumericVector loc_adj =  f(Named("MM")=MM, 
                             Named("start")=start, 
                             Named("end")=end,
                             Named("adj")=adj);
  // as<arma::vec>(loc_adj)：将 NumericVector 类型的 loc_adj 转换为 arma::vec 类型的 loc_adj1。
  //arma::vec loc_adj1：将转换后的结果存储在 loc_adj1 中。
  arma::vec loc_adj1 = as<arma::vec>(loc_adj);
  
  
  //cout << "loc_adj1==" << loc_adj1 << endl;
  
  return loc_adj1;
} 


// 定义 Combination 函数
bool next_comb(int* comb, const int n, const int k) {
  int i = k - 1;
  const int e = n - k;
  do
    comb[i]++;
  while (comb[i] > e + i && i--);
  if (comb[0] > e)
    return 0;
  while (++i < k)
    comb[i] = comb[i - 1] + 1;
  return 1;
}

int Combination(int* cho, int nl, int n, int k) {
  
  if (n < k || k <= 0)
    return 0;
  int* comb = new int[k];
  for (int i = 0; i < k; i++)
    comb[i] = i;
  
  do
    for (int i = 0; i < k; ++i) {
      cho[nl] = comb[i];
      //cout << cho[nl] << endl;
      nl++;
    } 
    while (next_comb(comb, n, k));
  
  delete[] comb;
  return nl;
} 


// 定义 Cho_IV1 函数
//const int nodex1（子节点）：输入参数 nodex1 是一个常量整数，表示某个节点的索引或编号。
//const arma::vec& nodeadjust（调整节点）：输入参数 nodeadjust 是一个常量引用，类型为 arma::vec（Armadillo 库中的列向量），表示节点调整值。
//const arma::mat& bet1：输入参数 bet1 是一个常量引用，类型为 arma::mat（Armadillo 库中的矩阵），表示某个系数矩阵。
//arma::uvec& l19：输入参数 l19 是一个引用，类型为 arma::uvec（Armadillo 库中的无符号整数列向量），表示某个输出向量。
void Cho_IV1(const int nodex1, const arma::vec& nodeadjust, const arma::mat& bet1,arma::uvec& l19){
  //创建一个大小为 nodeadjust.n_elem + 1 的无符号整数列向量 posadj。
  // posadj：调整节点列向量
  uvec posadj(nodeadjust.n_elem+1);
  //将 nodeadjust 中的每个元素复制到 posadj 中。
  for(int hu=0; hu<(nodeadjust.n_elem); hu++){
    posadj[hu]=nodeadjust[hu];
  }  
  //将 nodex1 添加到 posadj 的最后一个位置。
  // 调整节点+子节点
  posadj[nodeadjust.n_elem]=nodex1;
  //挑选相应节点下所有工具变量的效应值：mstadj
  mat mstadj = bet1.cols(posadj);
  // 创建一个大小为 mstadj.n_rows 的列向量 ll1。
  vec ll1(mstadj.n_rows);
  //遍历 mstadj 的每一行，如果某一行中有非零元素，则将 ll1 中对应位置的值设为 1。
  for(int ko=0; ko<(mstadj.n_rows); ko++){
    if(any(mstadj.row(ko)!= 0)) ll1[ko]=1;
  }  
  // 使用 find 函数找到 ll1 中值为 1 的元素的索引，并将这些索引存储在 l19 中。
  l19=find(ll1==1);
  
}   


// # 定义 fastLm 函数
// 用于判断两个节点之间是否存在因果关系
Rcpp::List fastLm(const arma::mat& X1, const arma::colvec& y,bool intercept) {
  
  int n = X1.n_rows;
  
  mat X;
  if(intercept)
    X=join_rows(ones(n,1),X1);
  else
    X=X1;
  
  int k = X.n_cols;
  
  //arma::colvec coef = arma::spsolve(X, y);    // fit model y ~ X
  arma::mat coef = arma::pinv(X) * y;
  arma::colvec res  = y - X*coef;           // residuals
  
  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));
  
  arma::colvec tval = coef/std_err;
  
  NumericVector xx=as<NumericVector>(wrap(-abs(tval)));
  double df=n-k;
  
  NumericVector pval = Rcpp::pt(xx,df)*2;
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("tvalue")       = tval,
                            Rcpp::Named("pvalue")       = pval,
                            Rcpp::Named("df")  = n - k);
}  


// judge edge
int Judge(const arma::mat& betaX,const arma::vec& betaY,const arma::vec& seY,
          const int use_egger){
  
  int edge=0;
  mat sigma = diagmat(1/seY);
  mat betaX1 = sigma*betaX;
  mat betaY1 = sigma*betaY;
  
  if(use_egger==1){
    
    List lm1=fastLm(betaX1,betaY1,true);
    arma::vec pval1 = as<arma::vec>(lm1["pvalue"]);
    
    if(pval1[0]>0.05){
      List lm2=fastLm(betaX1,betaY1,false); 
      arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
      if(pval2[0]<0.05)
        edge=1;
    }else if(pval1[1]<0.05){ 
      edge=1; 
    }
    
  }else if(use_egger==0){ 
    
    List lm2=fastLm(betaX1,betaY1,false); 
    arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
    if(pval2[0]<0.05)
      edge=1;
    
  } 
  
  return edge;
  
}  


// [[Rcpp::export]]
// # data_sum_beta 所有J个SNP对边际因果图中D个变量的汇总系数的J*D矩阵，例如，每个变量对每个SNP的线性回归系数。
// data_sum_se 对应于data_sum_beta的J*D标准误差矩阵。
// beta: J*D的指示矩阵，元素为0/1，例如，第i行第j列的元素为1表示第i个SNP是第j个节点的工具变量。
// adj_methods 选择调整变量的三种方法：1 移除collier；2 开放路径中的变量；3 最小充分集。
// use_eggers_step2 第二步中使用的多变量孟德尔随机化方法：1 MVMR-Egger；0 MVMR-IVW。
// vary mvmr adj=0表示直接对mvmr中的调整变量进行调整；vary mvmr adj=1表示用adj_methods从选定的变量改变mvmr中调整变量的数量。
SEXP STEP2(SEXP amattt, SEXP data_sum_beta,SEXP data_sum_se,SEXP beta1, 
           SEXP adj_methods,SEXP use_eggers_step2,SEXP vary_mvmr_adj){
  
  arma::mat amatt = as<arma::mat>(amattt);
  const arma::mat sum_beta = as<arma::mat>(data_sum_beta);
  const arma::mat sum_se = as<arma::mat>(data_sum_se);
  const arma::mat beta = as<arma::mat>(beta1);
  
  const int adj_method=as<int>(adj_methods);
  const int use_egger_step2=as<int>(use_eggers_step2);
  const int vary_mvmr_adjnodes=as<int>(vary_mvmr_adj);
  
  
  int D1 = sum_beta.n_cols, D2 = sum_se.n_cols;
  if (D1 != D2){
    perror("The dimensions of beta and se are not matched");
  }  
  
  
  int* cho = new int[100000000];
  int nl=0;
  nl=Combination(cho,nl,D1,2);
  
  // step 2
  
  cout << "step2 all start the first direction--" << nl << "..." << endl;
  
  for(int i=0; i<(nl/2); i++){
    
    //cout << "step2--first " << i << endl;
    //第一步先将两个节点本省移除
    arma::vec nodesadj_all;
    // 边际因果图中存在的边
    if(amatt(cho[2*i],cho[2*i + 1])==1){
      
      if(adj_method==1){
        // 这个向量用于存储所有可能的节点。
        vec allnodes=regspace(0,D1-1);
        // cout << allnodes.n_elem << endl;
        
        // 移除节点本身
        uvec remit={0,0};;
        remit[0]=cho[2*i];
        remit[1]=cho[2*i+1];
        //shed_rows ：删除矩阵中指定的行。
        // allnodes去除了节点本身的集合
        allnodes.shed_rows(remit);
        
        //cout << "remove colliders" << endl：移除对撞点
        // 创建一个向量：-2是减掉了两个节点本身
        vec rem(D1-2);
        for(int k=0; k<(D1-2); k++){
          rem[k]=amatt(cho[2*i],allnodes[k])+amatt(cho[2*i+1],allnodes[k]);
        }
        //allnodes.elem：选择 allnodes 中指定索引的元素。
        // nodesadj_all不是对撞点的节点集合（调整集合）
        nodesadj_all=allnodes.elem(find(rem<2));
        
      }else{
        
        //cout << "problem" << endl;
        nodesadj_all = Open_Path_Adj(amatt,cho[2*i]+1,cho[2*i + 1]+1,
                                     adj_method);
        //cout << "cho[2*i]+1==" << cho[2*i]+1 << endl;
        //cout << "cho[2*i + 1]+1==" << cho[2*i + 1]+1 << endl;
      }
      
      // cout << "cpp--num of nodesadj==" << nodesadj.n_elem << endl;
      // 调整集为空集
      if(nodesadj_all.n_elem==0){
        //cout << "No variable for MVMR" << endl;
      }else{
        
        // MVMR（作多变量孟德尔）
        
        int edge_up;
        
        //cout << "MVMR... " << endl;
        
        if(vary_mvmr_adjnodes==1){
          //获取向量的元素数量，调整集中的节点数
          int nodesnum=nodesadj_all.n_elem;
          
          for(int k=0; k<nodesnum; k++){
            
            // int k=0;
            //cho2：存调整集中节点的组合
            int* cho2 = new int[100000000];
            int nl2=0;
            // 对调整集中节点做各种数量的组合
            nl2=Combination(cho2,nl2,nodesnum,k+1);
            
            //cout << "nodesnum k -- " << k << " nl2 " << nl2 << endl;
            
            // cout << cho[2*i] << " " << cho[2*i+1] << endl;
            
            // e:组合数
            for(long e=0; e<(nl2/(k+1)); e++){
              // 提取从 cho2[(k+1)*e] 行到 cho2[(k+1)*e+k]行的子矩阵或子向量，并赋值nodesadj。
              // nodesadj：每个组合的所有玄素（节点）
              vec nodesadj=nodesadj_all.rows(cho2[(k+1)*e],cho2[(k+1)*e+k]);
              //cout << "nodesadj.n_elem " << nodesadj.n_elem << endl;
              uvec l3;
              // 选择工具变量
              // 每一次nodesadj只包含一个组合的节点
              Cho_IV1(cho[2*i],nodesadj,beta,l3);
              
              
              mat betaX3 = sum_beta.col(cho[2*i]);
              betaX3 = betaX3.elem(l3);
              vec betaY3 = sum_beta.col(cho[2*i+1]);
              betaY3 = betaY3.elem(l3);
              vec seY3 = sum_se.col(cho[2*i+1]);
              seY3 = seY3.elem(l3);
              
              vec choi;
              vec betaX4;
              mat betaX5=betaX3;
              // 每个组合内还要对每个节点进行
              for(int ut=0; ut<(k+1); ut++){
                // choi：每个组合中的节点
                choi=nodesadj_all.row(cho2[(k+1)*e+ut]);
                int choi1=choi[0];
                //cout << "choi " << choi << endl;
                // betax4：调整集对应工具变量的效应
                betaX4 = sum_beta.col(choi1);
                betaX4 = betaX4.elem(l3);
                betaX5=join_rows(betaX5,betaX4);
              }
              
              edge_up=Judge(betaX5,betaY3,seY3,use_egger_step2);
              // 边都没有效应
              if(edge_up==0){
                amatt(cho[2*i],cho[2*i+1])=0;
                k=nodesnum;
                break;
              }
              // cout <<  "edge_up " << edge_up << endl;
            }
            delete[] cho2;
          }
          
        }else{
          // 直接对每个节点作多变量孟德尔
          uvec l3;
          
          Cho_IV1(cho[2*i],nodesadj_all,beta,l3);
          
          
          
          mat betaX3 = sum_beta.col(cho[2*i]);
          betaX3 = betaX3.elem(l3);
          vec betaY3 = sum_beta.col(cho[2*i+1]);
          betaY3 = betaY3.elem(l3);
          vec seY3 = sum_se.col(cho[2*i+1]);
          seY3 = seY3.elem(l3);
          
          vec betaX4;
          mat betaX5=betaX3;
          for(int ut=0; ut<(nodesadj_all.n_elem); ut++){
            
            int choi1=nodesadj_all[ut];
            betaX4 = sum_beta.col(choi1);
            betaX4 = betaX4.elem(l3);
            betaX5=join_rows(betaX5,betaX4);
          }
          
          edge_up=Judge(betaX5,betaY3,seY3,use_egger_step2);
          if(edge_up==0){
            amatt(cho[2*i],cho[2*i+1])=0;
          }
        }
        
        
        
      }
      
    }
    
  }
  
  
  //第二个方向
  cout << "step2 all start the second direction" << nl << "..." << endl;
  for(int i=0; i<(nl/2); i++){
    
    //cout << "step2--second " << i << endl;
    arma::vec nodesadj_all;
    
    if(amatt(cho[2*i+1],cho[2*i])==1){
      
      
      if(adj_method==1){
        
        vec allnodes=regspace(0,D1-1);
        // cout << allnodes.n_elem << endl;
        
        // remove itself
        uvec remit={0,0};;
        remit[0]=cho[2*i];
        remit[1]=cho[2*i+1];
        allnodes.shed_rows(remit);
        
        //cout << "remove colliders" << endl;
        vec rem(D1-2);
        for(int k=0; k<(D1-2); k++){
          rem[k]=amatt(cho[2*i],allnodes[k])+amatt(cho[2*i+1],allnodes[k]);
        } 
        nodesadj_all=allnodes.elem(find(rem<2));
        
      }else{
        
        nodesadj_all = Open_Path_Adj(amatt,cho[2*i]+1,cho[2*i + 1]+1,
                                     adj_method);
      }
      
      //cout << "cpp--num of nodesadj==" << nodesadj.n_elem << endl;
      
      if(nodesadj_all.n_elem==0){
        //cout << "No variable for MVMR" << endl;
      }else{ 
        
        //cout << "MVMR..." << endl;
        int edge_up;
        
        if(vary_mvmr_adjnodes==1){
          
          int nodesnum=nodesadj_all.n_elem;
          
          for(int k=0; k<nodesnum; k++){
            
            // int k=0;
            
            int* cho2 = new int[100000000];
            int nl2=0;
            nl2=Combination(cho2,nl2,nodesnum,k+1);
            
            //cout << "nodesnum k -- " << k << " nl2 " << nl2 << endl;
            
            // cout << cho[2*i] << " " << cho[2*i+1] << endl;
            
            
            for(long e=0; e<(nl2/(k+1)); e++){
              
              vec nodesadj=nodesadj_all.rows(cho2[(k+1)*e],cho2[(k+1)*e+k]);
              //cout << "nodesadj.n_elem " << nodesadj.n_elem << endl;
              uvec l3;
              
              Cho_IV1(cho[2*i],nodesadj,beta,l3);
              
              
              mat betaX3 = sum_beta.col(cho[2*i]);
              betaX3 = betaX3.elem(l3);
              vec betaY3 = sum_beta.col(cho[2*i+1]);
              betaY3 = betaY3.elem(l3);
              vec seY3 = sum_se.col(cho[2*i+1]);
              seY3 = seY3.elem(l3);
              
              vec choi;
              vec betaX4;
              mat betaX5=betaX3;
              for(int ut=0; ut<(k+1); ut++){
                
                choi=nodesadj_all.row(cho2[(k+1)*e+ut]);
                int choi1=choi[0];
                //cout << "choi " << choi << endl;
                
                betaX4 = sum_beta.col(choi1);
                betaX4 = betaX4.elem(l3);
                betaX5=join_rows(betaX5,betaX4);
              } 
              // 最后得到得betaX5应该是除了子节点，还包含了每个组合中的所有节点
              edge_up=Judge(betaX5,betaY3,seY3,use_egger_step2);
              // 边都没效应，即两个A都是0
              if(edge_up==0){
                // 两个节点之间不存在边
                amatt(cho[2*i],cho[2*i+1])=0;
                k=nodesnum;
                break;
              } 
              // cout <<  "edge_up " << edge_up << endl;
            } 
            // 释放cho2 以得别的元素个数组合数
            delete[] cho2;
          } 
          
        }else{
          
          // MVMR
          
          uvec l3;
          
          Cho_IV1(cho[2*i],nodesadj_all,beta,l3);
          
          
          mat betaX3 = sum_beta.col(cho[2*i]);
          betaX3 = betaX3.elem(l3);
          vec betaY3 = sum_beta.col(cho[2*i+1]);
          betaY3 = betaY3.elem(l3);
          vec seY3 = sum_se.col(cho[2*i+1]);
          seY3 = seY3.elem(l3);
          
          vec betaX4;
          mat betaX5=betaX3;
          for(int ut=0; ut<(nodesadj_all.n_elem); ut++){
            
            int choi1=nodesadj_all[ut];
            betaX4 = sum_beta.col(choi1);
            betaX4 = betaX4.elem(l3);
            betaX5=join_rows(betaX5,betaX4);
          }
          
          edge_up=Judge(betaX5,betaY3,seY3,use_egger_step2);
          if(edge_up==0){
            amatt(cho[2*i],cho[2*i+1])=0;
          }
          
        }
        
        
        
      } 
      
    }
    
  }
  
  
  delete[] cho;
  
  return List::create(Rcpp::Named("amatt") = amatt);
}