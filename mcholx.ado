***********************MCholX Program*************************
*This function is a tentative translation of the "mchol" function 
*by Michael Zibulevsky and Brian Borchers to Stata

*Xintong 1/25/2015
  
*[L, D, E, pneg]=mchol1(G)
*
*Given a symmetric matrix G, find a matrix E of "small" norm and c,
*L and D such that G+E is Positive Definite, and 
*
*                    G+E=L*D*L'
*
*Also, calculate a direction pneg, such that if G is not PD, then
*
*                    pneg'*G*pneg<0
*
* Note that if G is PD, then the routine will return pneg=[]. 
*
* Reference: Gill, Murray, and Wright, "Practical Optimization", p111.
* Author: Brian Borchers (borchers@nmt.edu)

program define mcholx, eclass
 version 10.0
 marksample touse
 set more off
 
 matrix G=`1'
 
 mata: mcholX("'G'")
 mat L=e(L)
 mat D=e(D)
 mat E=e(E)
 mat pneg=e(pneg)
 
 ereturn mat L L
 ereturn mat D D
 ereturn mat E E
 ereturn mat pneg pneg

end

mata:

function mcholX(matrix G)
{

 G=st_matrix("G")
 
 //n gives the size of the matrix.
 
 n=rows(G)
 
 //gamma, zi, nu, and beta2 are quantities used by the algorithm.  
 
 gamma=max(diagonal(G))
 zi=max(max(G-diag(diagonal(G))))
 nu=max((1\ sqrt(n^2-1)))
 beta2=max((gamma\ zi/nu \ 1.0E-15))
 
 //Initialize diag(C) to diag(G).
 
 C=diag(diagonal(G))
 
 //Loop through, calculating column j of L for j=1:n
 L=J(n, n, 0)
 D=J(n, n, 0)
 E=J(n, n, 0)
 theta=J(n,1,.)
 
 for (j=1; j<=n; j++){
 
  bb=(1..j-1) 
  ee=(j+1..n)
   
  //
  //Calculate the jth row of L.
  //
  
  if (j>1) { 
    L[j,bb]=C[j,bb]:/(diagonal(D[bb,bb]))'
  }
  //
  //Update the jth column of C.
  //
  if (j>=2){
    if (j < n){
	 C[ee,j]=G[ee,j]-(L[j,bb]*C[ee,bb]')'
	}
   }
  else {
	//change 'ee' in the original code to 'n' 
	 if (j==1){
	 C[ee,j]=G[ee,j]
     }
   }
  
  //
  //Update theta.
  //
   
  if (j==n){
   theta[j]=0
   }
  else{
   theta[j]=max(abs(C[ee,j]))
   }
 //
 //Update D
 //
 D[j,j]=max((epsilon(1)\ abs(C[j,j]) \ theta[j]^2/beta2 ))
 
 //
 //Update E
 //
 E[j,j]=D[j,j]-C[j,j]
 
 //
 //Update C again...
 //
 for (i=j+1; i<=n; i++){
  C[i,i]=C[i,i]-C[i,j]^2/D[j,j]
 }
 
 
 }
 
 //
 //Put 1's on the diagonal of L
 //
 for (v=1;v<=n; v++){
  L[v,v]=1
 }
 
 
 //
 //if needed, find a descent direction.  
 //
 
 if (min(diagonal(C))<0){
  m=min(diagonal(C))
  minindex(diagonal(C),1,k,w)
  rhs=J(n,1,0)
  krow=rows(k)
  for (i=1; i<=krow; i++){  
  rhs[k[i]]=1
  }
  pneg=lusolve(L', rhs)
 }
 else{
 pneg=.
 }
 
 
 st_matrix("e(L)", L)
 st_matrix("e(E)", E)
 st_matrix("e(D)", D)
 st_matrix("e(pneg)", pneg)

}
 
 end
  



