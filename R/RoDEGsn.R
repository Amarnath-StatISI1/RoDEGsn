##' @title RoDEGsn function for finding deferentially expressed genes
##' @param exp_values: The gene expression data matrix
##' @param groups: The indicator vector for the groups. Here we take 1 and 0 to
##' indicate the two groups
##' @param d: Number of top ranked genes are to be found
##' @param alpha: Value of the tuning parameter
##' @description  Find out the top ranked deferentially expressed genes from
##' the set of genes in any noisy and skewed gene expression data. The underlying
##' distribution considered here is the skew normal(SN) distribution. The 
##' function calculated the adjusted p-values of the two sample robust Wald-type
##' tests and returns 'd' top ranked  genes based on the adjusted p-values
RoDEGsn=function(exp_values,groups,d,alpha)
{
  exp_values=as.matrix(exp_values)
  groups=as.vector(groups)
  d=as.numeric(d)
  alpha=as.numeric(alpha)
  ex_gp1=exp_values[,which(groups==1)]
  ex_gp2=exp_values[,which(groups==0)]
  ngene=nrow(ex_gp1)
  p_value=rep(0,d)
  ncores=detectCores()
  cl=makeCluster(ncores-2)
  registerDoParallel(cl)
  p_value=foreach(i=1:ngene,.combine=rbind,.packages = c("MASS","sn")) %dopar%
    {
      source("Two sample MDPDE sn.R");
      if(alpha==0)
      {p_value[i]=T_mle(g1=ex_gp1[i,],g2=ex_gp2[i,])
      }
      else{
        p_value[i]=T_dpd(g1=ex_gp1[i,],g2=ex_gp2[i,],alpha)
      }
    }
  stopCluster(cl)
  adjusted_p=p.adjust(p_value,method="fdr")
  results=cbind(exp_values,adjusted_p)
  rownames(results)=rownames(exp_values)
  n=ncol(results)
  results=results[base::order(results[,n]),]
  output=cbind("p_adj"=results[1:d,n])
  rownames(output)=rownames(results[1:d,])
  #results=results[,order(results[,1])]
  return(output)
}
