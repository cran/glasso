!
!
!           Lasso regularized covariance matrix estimate
!
!                     version (12/15/09)
!
! call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
!
! Input:
!    n = dimension of matrix
!    ss(n,n) = data covariance matrix
!    rho(n,n) = regularization strength parameters for each element
!              (must be symmetric: rho(i,j)=rho(j,i))
!    ia = approximation flag
!       ai =  0 => exact solution
!       ia != 0 => Meinhausen-Buhlmann approximation
!    is = initialization flag
!       is  = 0 => cold start: initialize using ss
!       is != 0 => warm start: initialize with previous solution
!                  stored in ww and wwi (see below)
!    itr = trace flag
!       itr != 0 => trace information printed
!       itr =  0 => trace information not printed
!    ipen = diagonal penalty flag
!       ipen != 0 => diagonal is penalized
!       ipen =  0 => diagonal is not penalized
!    thr = convergence threshold: iterations stop when average absolute
!          parameter change is less than thr * ave(abs(offdiag(ss)))
!          (suggested value 1.0e-4)
!    maxit = maximum number of iterations (no effect for ia ! = 0)
!
! Output:
!    ww(n,n) = solution covariance matrix estimate (ia = 0)
!               (not used for ia != 0)
!    wwi(n,n) = solution inverse covariance matrix estimate (ia = 0)
!             = off-diagonal lasso coefficients (ia != 0)
!    niter = number of iterations
!    del = average absolute parameter change at termination
!             (not used for ia != 0)
!    jerr = memory allocation error flag
!      jerr = 0 => no error
!      jerr != 0 => memory allocation error - no output returned
!
!
!              
      subroutine glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      parameter(eps=1.0e-7)                                                 
      real ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                                
      real, dimension (:,:), allocatable :: vv,xs                               
      real, dimension (:), allocatable :: s,x,z,ws,ro                           
      integer, dimension (:), allocatable :: mm                                 
      nm1=n-1                                                                   
      allocate(vv(1:nm1,1:nm1),stat=jerr)                                       
      ierr=0                                                                    
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=ierr)                             
      jerr=jerr+ierr                                                            
      allocate(s(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(x(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(z(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(mm(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      allocate(ro(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      if(ia .ne. 0)goto 10021                                               
      allocate(ws(1:n),stat=ierr)                                           
      jerr=jerr+ierr                                                        
10021 continue                                                              
      if(jerr.ne.0) return                                                  
      shr=0.0                                                               
10030 do 10031 j=1,n                                                        
10040 do 10041 k=1,n                                                        
      if(j.eq.k)goto 10041                                                  
      shr=shr+abs(ss(j,k))                                                  
10041 continue                                                              
10042 continue                                                              
10031 continue                                                              
10032 continue                                                              
      if(shr .ne. 0.0)goto 10061                                            
      ww=0.0                                                                
      wwi=0.0                                                               
10070 do 10071 j=1,n                                                        
      if(ipen .ne. 0)goto 10091                                             
      ww(j,j)=ss(j,j)                                                       
      goto 10101                                                            
10091 continue                                                              
      ww(j,j)=ss(j,j)+rho(j,j)                                              
10101 continue                                                              
10081 continue                                                              
      wwi(j,j)=1.0/max(ww(j,j),eps)                                         
10071 continue                                                              
10072 continue                                                              
      return                                                                
10061 continue                                                              
      shr=thr*shr/nm1                                                       
      if(ia .eq. 0)goto 10121                                               
      if(is.eq.0) wwi=0.0                                                   
10130 do 10131 m=1,n                                                        
      call setup(m,n,ss,rho,ss,vv,s,ro)                                     
      l=0                                                                   
10140 do 10141 j=1,n                                                        
      if(j.eq.m)goto 10141                                                  
      l=l+1                                                                 
      x(l)=wwi(j,m)                                                         
10141 continue                                                              
10142 continue                                                              
      call lasso(ro,nm1,vv,s,shr/n,x,z,mm)                                  
      l=0                                                                   
10150 do 10151 j=1,n                                                        
      if(j.eq.m)goto 10151                                                  
      l=l+1                                                                 
      wwi(j,m)=x(l)                                                         
10151 continue                                                              
10152 continue                                                              
10131 continue                                                              
10132 continue                                                              
      niter=1                                                               
      return                                                                
10121 continue                                                              
      if(is .ne. 0)goto 10171                                               
      ww=ss                                                                 
      xs=0.0                                                                
      goto 10181                                                            
10171 continue                                                              
10190 do 10191 j=1,n                                                        
      xjj=-2.0*wwi(j,j)                                                     
      l=0                                                                   
10200 do 10201 k=1,n                                                        
      if(k.eq.j)goto 10201                                                  
      l=l+1                                                                 
      xs(l,j)=wwi(k,j)/xjj                                                  
10201 continue                                                              
10202 continue                                                              
10191 continue                                                              
10192 continue                                                              
10181 continue                                                              
10161 continue                                                              
10210 do 10211 j=1,n                                                        
      if(ipen .ne. 0)goto 10231                                             
      ww(j,j)=ss(j,j)                                                       
      goto 10241                                                            
10231 continue                                                              
      ww(j,j)=ss(j,j)+rho(j,j)                                              
10241 continue                                                              
10221 continue                                                              
10211 continue                                                              
10212 continue                                                              
      niter=0                                                               
10250 continue                                                              
10251 continue                                                              
      dlx=0.0                                                               
10260 do 10261 m=1,n                                                        
      if(itr .eq. 0)goto 10281                                              
      write(6,10290)m                                                       
10290 format ('outer loop, m =',i10)                                        
10281 continue                                                              
      x=xs(:,m)                                                             
      ws=ww(:,m)                                                            
      call setup(m,n,ss,rho,ww,vv,s,ro)                                     
      call lasso(ro,nm1,vv,s,shr/sum(abs(vv)),x,z,mm)                       
      if(itr.eq.1) write(6,*) "threshold=",shr/sum(abs(vv))
      call cleanup(m,n,ww,vv,x,s,z,mm)                                      
      dlx=max(dlx,sum(abs(ww(:,m)-ws)))                                     
      xs(:,m)=x                                                             
10261 continue                                                              
10262 continue                                                              
      niter=niter+1                                                         
      if(niter.ge.maxit)goto 10252                                          
      if(dlx.lt.shr)goto 10252                                              
      goto 10251                                                            
10252 continue                                                              
      del=dlx/nm1                                                           
      call inv(n,ww,xs,wwi)                                                 
      return                                                                
      end                                                                   
      subroutine setup(m,n,ss,rho,ww,vv,s,r)                                
      real ss(n,n),rho(n,n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1)               
      l=0                                                                   
10300 do 10301 j=1,n                                                        
      if(j.eq.m)goto 10301                                                  
      l=l+1                                                                 
      r(l)=rho(j,m)                                                         
      s(l)=ss(j,m)                                                          
      i=0                                                                   
10310 do 10311 k=1,n                                                        
      if(k.eq.m)goto 10311                                                  
      i=i+1                                                                 
      vv(i,l)=2.0*ww(k,j)                                                   
10311 continue                                                              
10312 continue                                                              
10301 continue                                                              
10302 continue                                                              
      return                                                                
      end                                                                   
      subroutine lasso(rho,n,vv,s,thr,x,z,mm)                               
      real rho(n),vv(n,n),s(n),x(n),z(n)                                    
      integer mm(n)                                                         
      call fatmul(2,n,vv,x,s,z,mm)                                          
10320 continue                                                              
10321 continue                                                              
      dlx=0.0                                                               
10330 do 10331 j=1,n                                                        
      xj=x(j)                                                               
      x(j)=0.0                                                              
      t=s(j)+vv(j,j)*xj                                                     
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)          
      if(x(j).eq.xj)goto 10331                                              
      del=x(j)-xj                                                           
      dlx=max(dlx,abs(del))                                                 
      s=s-del*vv(:,j)                                                       
10331 continue                                                              
10332 continue                                                              
      if(dlx.lt.thr)goto 10322                                              
      goto 10321                                                            
10322 continue                                                              
      return                                                                
      end                                                                   
      subroutine cleanup(m,n,ww,vv,x,s,z,mm)                                
      real ww(n,n),vv(n-1,n-1),x(n-1),s(n-1),z(n-1)                         
      integer mm(n-1)                                                       
      call fatmul(1,n-1,vv,x,s,z,mm)                                        
      l=0                                                                   
10340 do 10341 j=1,n                                                        
      if(j.eq.m)goto 10341                                                  
      l=l+1                                                                 
      ww(j,m)=s(l)                                                          
      ww(m,j)=ww(j,m)                                                       
10341 continue                                                              
10342 continue                                                              
      return                                                                
      end                                                                   
      subroutine fatmul(it,n,vv,x,s,z,m)                                    
      parameter(fac=0.2)                                                    
      real vv(n,n),x(n),s(n),z(n)                                           
      integer m(n)                                                          
      l=0                                                                   
10350 do 10351 j=1,n                                                        
      if(x(j).eq.0.0)goto 10351                                             
      l=l+1                                                                 
      m(l)=j                                                                
      z(l)=x(j)                                                             
10351 continue                                                              
10352 continue                                                              
      if(l .le. int(fac*n))goto 10371                                       
      if(it .ne. 1)goto 10391                                               
      s=matmul(vv,x)                                                        
      goto 10401                                                            
10391 continue                                                              
      s=s-matmul(x,vv)                                                      
10401 continue                                                              
10381 continue                                                              
      goto 10361                                                            
10371 if(it .ne. 1)goto 10411                                               
10420 do 10421 j=1,n                                                        
      s(j)=dot_product(vv(j,m(1:l)),z(1:l))                                 
10421 continue                                                              
10422 continue                                                              
      goto 10431                                                            
10411 continue                                                              
10440 do 10441 j=1,n                                                        
      s(j)=s(j)-dot_product(vv(m(1:l),j),z(1:l))                            
10441 continue                                                              
10442 continue                                                              
10431 continue                                                              
10361 continue                                                              
      return                                                                
      end                                                                   
      subroutine inv(n,ww,xs,wwi)                                           
      real ww(n,n),xs(n-1,n),wwi(n,n)                                       
      nm1=n-1                                                               
      xs=-2.0*xs                                                            
      wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)))                 
      wwi(2:n,1)=wwi(1,1)*xs(:,1)                                           
      wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)))               
      wwi(1:nm1,n)=wwi(n,n)*xs(:,n)                                         
10450 do 10451 j=2,nm1                                                      
      jm1=j-1                                                               
      jp1=j+1                                                               
      temp1=dot_product(xs(1:jm1,j),ww(1:jm1,j))
      temp2=dot_product(xs(j:nm1,j),ww(jp1:n,j))
      wwi(j,j)=1.0/(ww(j,j)+temp1+temp2)
      wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j)                                     
      wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j)                                     
10451 continue                                                              
10452 continue                                                              
      return                                                                
      end                                                                   




! Output from Public domain Ratfor, version 1.0
!     coded by R tibshirani dec 2008
      subroutine glassopath(beta,what,jerrs,rholist, nrho,n,ss,rho,ia,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      integer nrho,n,jerrs(nrho)
      real rholist(nrho),beta(n,n,nrho),what(n,n,nrho)
      real ss(n,n),rho(n,n),ww(n,n),wwi(n,n)
      is=0
      do23000 j=1,n
      do23002 k=1,n
      rho(j,k)=rholist(nrho)
23002 continue
23003 continue
23000 continue
23001 continue
      call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      jerrs(1)=jerr
      do23004 j=1,n
      do23006 k=1,n
      beta(j,k,nrho)=wwi(j,k)
      what(j,k,nrho)=ww(j,k)
23006 continue
23007 continue
23004 continue
23005 continue
      is=1
      do23008 i =(nrho-1), 1,-1
      do23010 j=1,n
      do23012 k=1,n
      rho(j,k)=rholist(i)
23012 continue
23013 continue
23010 continue
23011 continue
      call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      jerrs(i)=jerr
      do23014 j=1,n
      do23016 k=1,n
      beta(j,k,i)=wwi(j,k)
      what(j,k,i)=ww(j,k)
23016 continue
23017 continue
23014 continue
23015 continue
23008 continue
23009 continue
      return
      end
