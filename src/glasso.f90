                                             
!
!           Lasso regularized covariance matrix estimate
!
!                         version (4/5/
!
! call lasinv(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
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
      subroutine glasso(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,niter,del,jerr)    
      real sss(nn,nn),rrho(nn,nn),www(nn,nn),wwwi(nn,nn)                    
      real, dimension (:,:), allocatable :: ss,rho,ww,wwi                       
      integer, dimension (:), allocatable :: in                                 
      allocate(in(1:nn),stat=jerr)                                          
      if(jerr.ne.0) return                                                  
      n=0                                                                   
10010 do 10011 k=1,nn                                                       
10020 do 10021 j=1,nn                                                       
      if(j.eq.k)goto 10021                                                  
      if(abs(sss(j,k)).le.rrho(j,k))goto 10021                              
      n=n+1                                                                 
      in(n)=k                                                               
      goto 10022                                                            
10021 continue                                                              
10022 continue                                                              
10011 continue                                                              
10012 continue                                                              
      allocate(ss(1:n,1:n),stat=jerr);                                          
      allocate(rho(1:n,1:n),stat=ierr); jerr=jerr+ierr                          
      allocate(ww(1:n,1:n),stat=ierr); jerr=jerr+ierr                           
      allocate(wwi(1:n,1:n),stat=ierr); jerr=jerr+ierr                          
      if(jerr.ne.0) return                                                  
10030 do 10031 k=1,n                                                        
10040 do 10041 j=1,n                                                        
      ss(j,k)=sss(in(j),in(k))                                              
      rho(j,k)=rrho(in(j),in(k))                                            
      ww(j,k)=www(in(j),in(k))                                              
      wwi(j,k)=wwwi(in(j),in(k))                                            
10041 continue                                                              
10042 continue                                                              
10031 continue                                                              
10032 continue                                                              
      call lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)    
      if(jerr.ne.0) return                                                  
      www=0.0                                                               
      wwwi=0.0                                                              
10050 do 10051 k=1,n                                                        
10060 do 10061 j=1,n                                                        
      wwwi(in(j),in(k))=wwi(j,k)                                            
10061 continue                                                              
10062 continue                                                              
10051 continue                                                              
10052 continue                                                              
      if(ia .ne. 0)goto 10081                                               
10090 do 10091 k=1,n                                                        
10100 do 10101 j=1,n                                                        
      www(in(j),in(k))=ww(j,k)                                              
10101 continue                                                              
10102 continue                                                              
10091 continue                                                              
10092 continue                                                              
10110 do 10111 j=1,nn                                                       
      if(www(j,j).ne.0.0)goto 10111                                         
      if(ipen .ne. 0)goto 10131                                             
      www(j,j)=sss(j,j)                                                     
      goto 10141                                                            
10131 continue                                                              
      www(j,j)=sss(j,j)+rrho(j,j)                                           
10141 continue                                                              
10121 continue                                                              
      wwwi(j,j)=1.0/www(j,j)                                                
10111 continue                                                              
10112 continue                                                              
10081 continue                                                              
      return                                                                
      end                                                                   
      subroutine lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      parameter(eps=1.0e-7)                                                 
      real ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                                
      real, dimension (:,:), allocatable :: vv,xs                               
      real, dimension (:), allocatable :: s,x,z,ws,ro,so                        
      integer, dimension (:), allocatable :: mm                                 
      nm1=n-1                                                                   
      allocate(vv(1:nm1,1:nm1),stat=jerr)                                       
      ierr=0                                                                    
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=ierr)                             
      jerr=jerr+ierr                                                            
      allocate(s(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(so(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      allocate(x(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(z(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(mm(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      allocate(ro(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      if(ia .ne. 0)goto 10161                                               
      allocate(ws(1:n),stat=ierr)                                           
      jerr=jerr+ierr                                                        
10161 continue                                                              
      if(jerr.ne.0) return                                                  
      shr=0.0                                                               
10170 do 10171 j=1,n                                                        
10180 do 10181 k=1,n                                                        
      if(j.eq.k)goto 10181                                                  
      shr=shr+abs(ss(j,k))                                                  
10181 continue                                                              
10182 continue                                                              
10171 continue                                                              
10172 continue                                                              
      if(shr .ne. 0.0)goto 10201                                            
      ww=0.0                                                                
      wwi=0.0                                                               
10210 do 10211 j=1,n                                                        
      if(ipen .ne. 0)goto 10231                                             
      ww(j,j)=ss(j,j)                                                       
      goto 10241                                                            
10231 continue                                                              
      ww(j,j)=ss(j,j)+rho(j,j)                                              
10241 continue                                                              
10221 continue                                                              
      wwi(j,j)=1.0/max(ww(j,j),eps)                                         
10211 continue                                                              
10212 continue                                                              
      return                                                                
10201 continue                                                              
      shr=thr*shr/nm1                                                       
      if(ia .eq. 0)goto 10261                                               
      if(is.eq.0) wwi=0.0                                                   
10270 do 10271 m=1,n                                                        
      call setup(m,n,ss,rho,ss,vv,s,ro)                                     
      l=0                                                                   
10280 do 10281 j=1,n                                                        
      if(j.eq.m)goto 10281                                                  
      l=l+1                                                                 
      x(l)=wwi(j,m)                                                         
10281 continue                                                              
10282 continue                                                              
      call lasso(ro,nm1,vv,s,shr/n,x,z,mm)                                  
      l=0                                                                   
10290 do 10291 j=1,n                                                        
      if(j.eq.m)goto 10291                                                  
      l=l+1                                                                 
      wwi(j,m)=x(l)                                                         
10291 continue                                                              
10292 continue                                                              
10271 continue                                                              
10272 continue                                                              
      niter=1                                                               
      return                                                                
10261 continue                                                              
      if(is .ne. 0)goto 10311                                               
      ww=ss                                                                 
      xs=0.0                                                                
      goto 10321                                                            
10311 continue                                                              
10330 do 10331 j=1,n                                                        
      xjj=-wwi(j,j)                                                         
      l=0                                                                   
10340 do 10341 k=1,n                                                        
      if(k.eq.j)goto 10341                                                  
      l=l+1                                                                 
      xs(l,j)=wwi(k,j)/xjj                                                  
10341 continue                                                              
10342 continue                                                              
10331 continue                                                              
10332 continue                                                              
10321 continue                                                              
10301 continue                                                              
10350 do 10351 j=1,n                                                        
      if(ipen .ne. 0)goto 10371                                             
      ww(j,j)=ss(j,j)                                                       
      goto 10381                                                            
10371 continue                                                              
      ww(j,j)=ss(j,j)+rho(j,j)                                              
10381 continue                                                              
10361 continue                                                              
10351 continue                                                              
10352 continue                                                              
      niter=0                                                               
10390 continue                                                              
10391 continue                                                              
      dlx=0.0                                                               
10400 do 10401 m=1,n                                                        
      if(itr .eq. 0)goto 10421                                              
      write(6,10430)m                                                       
10430 format ('outer loop, m =',i10)                                        
10421 continue                                                              
      x=xs(:,m)                                                             
      ws=ww(:,m)                                                            
      call setup(m,n,ss,rho,ww,vv,s,ro)                                     
      so=s                                                                  
      call lasso(ro,nm1,vv,s,shr/sum(abs(vv)),x,z,mm)                       
      l=0                                                                   
10440 do 10441 j=1,n                                                        
      if(j.eq.m)goto 10441                                                  
      l=l+1                                                                 
      ww(j,m)=so(l)-s(l)                                                    
      ww(m,j)=ww(j,m)                                                       
10441 continue                                                              
10442 continue                                                              
      dlx=max(dlx,sum(abs(ww(:,m)-ws)))                                     
      xs(:,m)=x                                                             
10401 continue                                                              
10402 continue                                                              
      niter=niter+1                                                         
      if(niter.ge.maxit)goto 10392                                          
      if(dlx.lt.shr)goto 10392                                              
      goto 10391                                                            
10392 continue                                                              
      del=dlx/nm1                                                           
      call inv(n,ww,xs,wwi)                                                 
      return                                                                
      end                                                                   
      subroutine setup(m,n,ss,rho,ww,vv,s,r)                                
      real ss(n,n),rho(n,n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1)               
      l=0                                                                   
10450 do 10451 j=1,n                                                        
      if(j.eq.m)goto 10451                                                  
      l=l+1                                                                 
      r(l)=rho(j,m)                                                         
      s(l)=ss(j,m)                                                          
      i=0                                                                   
10460 do 10461 k=1,n                                                        
      if(k.eq.m)goto 10461                                                  
      i=i+1                                                                 
      vv(i,l)=ww(k,j)                                                       
10461 continue                                                              
10462 continue                                                              
10451 continue                                                              
10452 continue                                                              
      return                                                                
      end                                                                   
      subroutine lasso(rho,n,vv,s,thr,x,z,mm)                               
      real rho(n),vv(n,n),s(n),x(n),z(n),t                                   
      integer mm(n)                                                         
      call fatmul(2,n,vv,x,s,z,mm)                                          
10470 continue                                                              
10471 continue                                                              
      dlx=0.0                                                               
10480 do 10481 j=1,n                                                        
      xj=x(j)                                                               
      x(j)=0.0                                                              
      t=s(j)+vv(j,j)*xj                                                     
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)          
      if(x(j).eq.xj)goto 10481                                              
      del=x(j)-xj                                                           
      dlx=max(dlx,abs(del))                                                 
      s=s-del*vv(:,j)                                                       
10481 continue                                                              
10482 continue                                                              
      if(dlx.lt.thr)goto 10472                                              
      goto 10471                                                            
10472 continue                                                              
      return                                                                
      end                                                                   
      subroutine fatmul(it,n,vv,x,s,z,m)                                    
      parameter(fac=0.2)                                                    
      real vv(n,n),x(n),s(n),z(n)                                           
      integer m(n)                                                          
      l=0                                                                   
10490 do 10491 j=1,n                                                        
      if(x(j).eq.0.0)goto 10491                                             
      l=l+1                                                                 
      m(l)=j                                                                
      z(l)=x(j)                                                             
10491 continue                                                              
10492 continue                                                              
      if(l .le. int(fac*n))goto 10511                                       
      if(it .ne. 1)goto 10531                                               
      s=matmul(vv,x)                                                        
      goto 10541                                                            
10531 continue                                                              
      s=s-matmul(x,vv)                                                      
10541 continue                                                              
10521 continue                                                              
      goto 10501                                                            
10511 if(it .ne. 1)goto 10551                                               
10560 do 10561 j=1,n                                                        
      s(j)=dot_product(vv(j,m(1:l)),z(1:l))                                 
10561 continue                                                              
10562 continue                                                              
      goto 10571                                                            
10551 continue                                                              
10580 do 10581 j=1,n                                                        
      s(j)=s(j)-dot_product(vv(m(1:l),j),z(1:l))                            
10581 continue                                                              
10582 continue                                                              
10571 continue                                                              
10501 continue                                                              
      return                                                                
      end                                                                   
      subroutine inv(n,ww,xs,wwi)                                           
      real ww(n,n),xs(n-1,n),wwi(n,n)                                       
      nm1=n-1                                                               
      xs=-xs                                                                
      wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)))                 
      wwi(2:n,1)=wwi(1,1)*xs(:,1)                                           
      wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)))               
      wwi(1:nm1,n)=wwi(n,n)*xs(:,n)                                         
10590 do 10591 j=2,nm1                                                      
      jm1=j-1                                                               
      jp1=j+1                                                               
      wwi(j,j)=1.0/(ww(j,j)+dot_product(xs(1:jm1,j),ww(1:jm1,j))  +dot_product(xs(j:nm1,j),ww(jp1:n,j)))    
      wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j)                                     
      wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j)                                     
10591 continue                                                              
10592 continue                                                              
      return                                                                
      end                                                                   

! Output from Public domain Ratfor, version 1.0
!     coded by R tibshirani dec 2008

      subroutine glassopath(beta,what,jerrs,rholist, nrho,n,ss,rho,ia,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
!      implicit double precision(a-h, o-z)

      integer nrho,n,jerrs(nrho)
!      double precision rholist(nrho),beta(n,n,nrho),what(n,n,nrho)
!      double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n)
      real rholist(nrho),beta(n,n,nrho),what(n,n,nrho)
      real ss(n,n),rho(n,n),ww(n,n),wwi(n,n)
      is=0
      do 23000 j=1,n
      do 23002 k=1,n
      rho(j,k)=rholist(nrho)
23002 continue
23003 continue
23000 continue
23001 continue
      call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      jerrs(1)=jerr
      do 23004 j=1,n
      do 23006 k=1,n
      beta(j,k,nrho)=wwi(j,k)
      what(j,k,nrho)=ww(j,k)
23006 continue
23007 continue
23004 continue
23005 continue
      is=1
      do 23008 i =(nrho), 1,-1
      do 23010 j=1,n
      do 23012 k=1,n
      rho(j,k)=rholist(i)
23012 continue
23013 continue
23010 continue
23011 continue
      if(itr.gt.0) write(6,*) "rho=",rholist(i)
      itr2=itr
      if(itr2.gt.0) itr2=itr-1
      call glasso(n,ss,rho,ia,is,itr2,ipen,thr,maxit,ww,wwi,niter,del,jerr)
      jerrs(i)=jerr
      do 23014 j=1,n
      do 23016 k=1,n
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
