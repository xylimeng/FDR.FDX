#for solving the constrained optimization
UpdateLambda=function(sol.structure,lambdas){
  p=length(sol.structure)
  i=1
  r=sum(sol.structure!=0)
  while (i<r) {
    right.ind=max(which(sol.structure==sol.structure[i]))
    if(right.ind>i){
      lambdas[i:right.ind]=mean(lambdas[i:right.ind])
    }
    i=right.ind+1
  }
  if(r<p){
    lambdas[(r+1):p]=0
  }
  return(lambdas)
}

#===helper func: solving quadratic equation ax^2+bx+c=0
quadSolver <- function(a,b,c){
  if((b^2-4*a*c) > 0){ 
    x_1 = (-b+sqrt(b^2-4*a*c))/(2*a)
    x_2 = (-b-sqrt(b^2-4*a*c))/(2*a)
    #return(c(x_1,x_2))
    if(x_1>=0 & x_2>=0){
      return(min(x_1,x_2))
    }
    else{
      return(max(x_1,x_2))
    }
  }
  else if((b^2-4*a*c) == 0){
    x = -b/(2*a)
    return(x)
  }
  else {return(NA)}
  
}


numRej_forAllC=function(betahat.absort,lambdas,n,sigma,tau,endpoint,by=0.5){
  Calpha.list=seq(from=0.01, to=endpoint,by=by)
  result=rep(0,length(Calpha.list))
  p=length(betahat.absort)
  result[1]=p
  sol.structure=1:p
  delta.tp=sqrt(Calpha.list[1]*sigma^2/(n+tau)/sum(lambdas^2))
  sol.tp=betahat.absort-delta.tp*lambdas
  delta=rep(0,length(Calpha.list))
  delta[1]=delta.tp
  for (i in 2:length(Calpha.list)) {
    Calpha.tp=Calpha.list[i]
    epsilon=quadSolver(a=sum(lambdas^2),b=2*sum((betahat.absort-sol.tp)*lambdas),
                       c=-(Calpha.tp-Calpha.list[i-1])*sigma^2/(n+tau))
    delta.tp=delta.tp+epsilon
    sol.tp=sol.tp-lambdas*epsilon
    if(sum(diff(sol.tp)>0)==0 & sum(sol.tp<0)==0){
      result[i]=sum(sol.structure>0)
      delta[i]=delta.tp
    }else{
      while(sum(diff(sol.tp)>0)>0 | sum(sol.tp<0)>0){
        if(sum(diff(sol.tp)>0)>0){
          ind=which(diff(sol.tp)>0)[1]
          ind.right=max(which(sol.structure==sol.structure[ind+1]))
          sol.structure[ind:ind.right]=mean(sol.structure[ind:ind.right])
          feasible.sol=sol.tp
          feasible.sol[ind:ind.right]=mean(feasible.sol[ind:ind.right])
        }
        if(sum(sol.tp<0)>0){
          ind=min(which(sol.tp<0))
          sol.structure[ind:p]=0
          feasible.sol=sol.tp
          feasible.sol[ind:p]=0
          
        }
        lambdas=UpdateLambda(sol.structure,lambdas)
        sdf=sum((betahat.absort-feasible.sol)^2)
        epsilon=quadSolver(a=sum(lambdas^2),b=2*sum((betahat.absort-feasible.sol)*lambdas),
                           c=-Calpha.tp*sigma^2/(n+tau)+sdf)
        
        sol.tp=feasible.sol-lambdas*epsilon
        delta[i]=delta.tp+epsilon
        
        if(sum(diff(sol.tp)>0)==0 & sum(sol.tp<0)==0){
          result[i]=sum(sol.structure>0)
        }
      }
    }
  }
  return(list("nRejs"=result,"deltas"=delta))
}