c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine glasso(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nn
     *iter,ddel,jerr)
      implicit double precision(a-h,o-z)                                
      double precision sss(nn,nn),rrho(nn,nn),www(nn,nn),wwwi(nn,nn)    
      double precision, dimension (:), allocatable :: ss,rho,ww,wwi     
      integer, dimension (:), allocatable :: ir,ie                      
      integer, dimension (:,:), allocatable :: ic                       
      if(ia .eq. 0)goto 10021                                           
      call lasinv1(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nniter,
     *ddel,jerr)
      return                                                            
10021 continue                                                          
      allocate(ic(1:2,1:nn),stat=jerr)                                  
      if (jerr.ne.0) return                                             
      allocate(ir(1:nn),stat=jerr)                                      
      if (jerr.ne.0) return                                             
      allocate(ie(1:nn),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call connect(nn,sss,rrho,nc,ic,ir,ie)                             
      nnq=0                                                             
      do 10031 kc=1,nc                                                  
      nnq=max(ic(2,kc)-ic(1,kc)+1,nnq)                                  
10031 continue                                                          
      continue                                                          
      nnq=nnq**2                                                        
      allocate(ss(1:nnq),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(rho(1:nnq),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(ww(1:nnq),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(wwi(1:nnq),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      nniter=0                                                          
      ddel=0.0                                                          
      l=0                                                               
      do 10041 kc=1,nc                                                  
      n=ic(2,kc)-ic(1,kc)+1                                             
      if(n .gt. 1)goto 10061                                            
      k=ir(ic(1,kc))                                                    
      www(:,k)=0.0                                                      
      www(k,:)=0.0                                                      
      wwwi(:,k)=0.0                                                     
      wwwi(k,:)=0.0                                                     
      goto 10041                                                        
10061 continue                                                          
      kb=ic(1,kc)                                                       
      ke=ic(2,kc)                                                       
      l=0                                                               
      do 10071 k=kb,ke                                                  
      ik=ir(k)                                                          
      do 10081 j=kb,ke                                                  
      ij=ir(j)                                                          
      l=l+1                                                             
      ss(l)=sss(ij,ik)                                                  
      rho(l)=rrho(ij,ik)                                                
      ww(l)=www(ij,ik)                                                  
      wwi(l)=wwwi(ij,ik)                                                
10081 continue                                                          
      continue                                                          
10071 continue                                                          
      continue                                                          
      call lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,je
     *rr)
      if(jerr.ne.0) return                                              
      nniter=nniter+niter                                               
      ddel=ddel+del                                                     
      do 10091 j=kb,ke                                                  
      k=ir(j)                                                           
      www(:,k)=0.0                                                      
      www(k,:)=0.0                                                      
      wwwi(:,k)=0.0                                                     
      wwwi(k,:)=0.0                                                     
10091 continue                                                          
      continue                                                          
      l=0                                                               
      do 10101 k=kb,ke                                                  
      ik=ir(k)                                                          
      do 10111 j=kb,ke                                                  
      l=l+1                                                             
      wwwi(ir(j),ik)=wwi(l)                                             
10111 continue                                                          
      continue                                                          
10101 continue                                                          
      continue                                                          
      if(ia .ne. 0)goto 10131                                           
      l=0                                                               
      do 10141 k=kb,ke                                                  
      ik=ir(k)                                                          
      do 10151 j=kb,ke                                                  
      l=l+1                                                             
      www(ir(j),ik)=ww(l)                                               
10151 continue                                                          
      continue                                                          
10141 continue                                                          
      continue                                                          
10131 continue                                                          
10041 continue                                                          
      continue                                                          
      ddel=ddel/nc                                                      
      if(ia .ne. 0)goto 10171                                           
      do 10181 j=1,nn                                                   
      if(www(j,j).ne.0.0)goto 10181                                     
      if(ipen .ne. 0)goto 10201                                         
      www(j,j)=sss(j,j)                                                 
      goto 10211                                                        
10201 continue                                                          
      www(j,j)=sss(j,j)+rrho(j,j)                                       
10211 continue                                                          
      continue                                                          
      wwwi(j,j)=1.0/www(j,j)                                            
10181 continue                                                          
      continue                                                          
10171 continue                                                          
      return                                                            
      end                                                               
      subroutine connect(n,ss,rho,nc,ic,ir,ie)                          
      implicit double precision(a-h,o-z)                                
      double precision ss(n,n),rho(n,n)                                 
      integer ic(2,n),ir(n),ie(n)                                       
      ie=0                                                              
      nc=0                                                              
      is=1                                                              
      do 10221 k=1,n                                                    
      if(ie(k).gt.0)goto 10221                                          
      ir(is)=k                                                          
      nc=nc+1                                                           
      ie(k)=nc                                                          
      ic(1,nc)=is                                                       
      is=is+1                                                           
      call row(nc,1,ir((is-1):n),n,ss,rho,ie,na,ir(is:n))               
      if(na .ne. 0)goto 10241                                           
      ic(2,nc)=is-1                                                     
      goto 10221                                                        
10241 continue                                                          
      continue                                                          
10251 continue                                                          
      nas=na                                                            
      iss=is                                                            
      il=iss+nas-1                                                      
      if(il.ge.n)goto 10252                                             
      is=is+na                                                          
      call row(nc,nas,ir(iss:n),n,ss,rho,ie,na,ir(is:n))                
      if(na.eq.0)goto 10252                                             
      goto 10251                                                        
10252 continue                                                          
      ic(2,nc)=il                                                       
10221 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine row(nc,nr,jr,n,ss,rho,ie,na,kr)                        
      implicit double precision(a-h,o-z)                                
      double precision ss(n,n),rho(n,n)                                 
      integer jr(nr),ie(n),kr(*)                                        
      na=0                                                              
      do 10261 l=1,nr                                                   
      k=jr(l)                                                           
      do 10271 j=1,n                                                    
      if(ie(j).gt.0)goto 10271                                          
      if(j.eq.k)goto 10271                                              
      if(abs(ss(j,k)).le.rho(j,k))goto 10271                            
      na=na+1                                                           
      kr(na)=j                                                          
      ie(j)=nc                                                          
10271 continue                                                          
      continue                                                          
10261 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,
     *del,jerr)
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0e-7)                                             
      double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                
      double precision, dimension (:,:), allocatable :: vv,xs           
      double precision, dimension (:), allocatable :: s,x,z,ws,ro,so    
      integer, dimension (:), allocatable :: mm                         
      nm1=n-1                                                           
      allocate(vv(1:nm1,1:nm1),stat=jerr)                               
      if(jerr.ne.0) return                                              
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=jerr)                     
      if(jerr.ne.0) return                                              
      allocate(s(1:nm1),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(so(1:nm1),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(x(1:nm1),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(z(1:nm1),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:nm1),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(ro(1:nm1),stat=jerr)                                     
      if(ia .ne. 0)goto 10291                                           
      allocate(ws(1:n),stat=jerr)                                       
10291 continue                                                          
      if(jerr.ne.0) return                                              
      shr=0.0                                                           
      do 10301 j=1,n                                                    
      do 10311 k=1,n                                                    
      if(j.eq.k)goto 10311                                              
      shr=shr+abs(ss(j,k))                                              
10311 continue                                                          
      continue                                                          
10301 continue                                                          
      continue                                                          
      if(shr .ne. 0.0)goto 10331                                        
      ww=0.0                                                            
      wwi=0.0                                                           
      do 10341 j=1,n                                                    
      if(ipen .ne. 0)goto 10361                                         
      ww(j,j)=ss(j,j)                                                   
      goto 10371                                                        
10361 continue                                                          
      ww(j,j)=ss(j,j)+rho(j,j)                                          
10371 continue                                                          
      continue                                                          
      wwi(j,j)=1.0/max(ww(j,j),eps)                                     
10341 continue                                                          
      continue                                                          
      return                                                            
10331 continue                                                          
      shr=thr*shr/nm1                                                   
      if(ia .eq. 0)goto 10391                                           
      if(is.eq.0) wwi=0.0                                               
      do 10401 m=1,n                                                    
      call setup(m,n,ss,rho,ss,vv,s,ro)                                 
      l=0                                                               
      do 10411 j=1,n                                                    
      if(j.eq.m)goto 10411                                              
      l=l+1                                                             
      x(l)=wwi(j,m)                                                     
10411 continue                                                          
      continue                                                          
      call lasso(ro,nm1,vv,s,shr/n,x,z,mm)                              
      l=0                                                               
      do 10421 j=1,n                                                    
      if(j.eq.m)goto 10421                                              
      l=l+1                                                             
      wwi(j,m)=x(l)                                                     
10421 continue                                                          
      continue                                                          
10401 continue                                                          
      continue                                                          
      niter=1                                                           
      return                                                            
10391 continue                                                          
      if(is .ne. 0)goto 10441                                           
      ww=ss                                                             
      xs=0.0                                                            
      goto 10451                                                        
10441 continue                                                          
      do 10461 j=1,n                                                    
      xjj=-wwi(j,j)                                                     
      l=0                                                               
      do 10471 k=1,n                                                    
      if(k.eq.j)goto 10471                                              
      l=l+1                                                             
      xs(l,j)=wwi(k,j)/xjj                                              
10471 continue                                                          
      continue                                                          
10461 continue                                                          
      continue                                                          
10451 continue                                                          
      continue                                                          
      do 10481 j=1,n                                                    
      if(ipen .ne. 0)goto 10501                                         
      ww(j,j)=ss(j,j)                                                   
      goto 10511                                                        
10501 continue                                                          
      ww(j,j)=ss(j,j)+rho(j,j)                                          
10511 continue                                                          
      continue                                                          
10481 continue                                                          
      continue                                                          
      niter=0                                                           
      continue                                                          
10521 continue                                                          
      dlx=0.0                                                           
      do 10531 m=1,n                                                    
      if(itr .eq. 0)goto 10551                                          
      call intpr('m',1,m,1)                                             
10551 continue                                                          
      x=xs(:,m)                                                         
      ws=ww(:,m)                                                        
      call setup(m,n,ss,rho,ww,vv,s,ro)                                 
      so=s                                                              
      call lasso(ro,nm1,vv,s,shr/sum(abs(vv)),x,z,mm)                   
      l=0                                                               
      do 10561 j=1,n                                                    
      if(j.eq.m)goto 10561                                              
      l=l+1                                                             
      ww(j,m)=so(l)-s(l)                                                
      ww(m,j)=ww(j,m)                                                   
10561 continue                                                          
      continue                                                          
      dlx=max(dlx,sum(abs(ww(:,m)-ws)))                                 
      xs(:,m)=x                                                         
10531 continue                                                          
      continue                                                          
      niter=niter+1                                                     
      if(niter.ge.maxit)goto 10522                                      
      if(dlx.lt.shr)goto 10522                                          
      goto 10521                                                        
10522 continue                                                          
      del=dlx/nm1                                                       
      call inv(n,ww,xs,wwi)                                             
      return                                                            
      end                                                               
      subroutine setup(m,n,ss,rho,ww,vv,s,r)                            
      implicit double precision(a-h,o-z)                                
      double precision ss(n,n),rho(n,n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1
     *)
      l=0                                                               
      do 10571 j=1,n                                                    
      if(j.eq.m)goto 10571                                              
      l=l+1                                                             
      r(l)=rho(j,m)                                                     
      s(l)=ss(j,m)                                                      
      i=0                                                               
      do 10581 k=1,n                                                    
      if(k.eq.m)goto 10581                                              
      i=i+1                                                             
      vv(i,l)=ww(k,j)                                                   
10581 continue                                                          
      continue                                                          
10571 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lasso(rho,n,vv,s,thr,x,z,mm)                           
      implicit double precision(a-h,o-z)                                
      double precision rho(n),vv(n,n),s(n),x(n),z(n)                    
      integer mm(n)                                                     
      call fatmul(2,n,vv,x,s,z,mm)                                      
      continue                                                          
10591 continue                                                          
      dlx=0.0                                                           
      do 10601 j=1,n                                                    
      xj=x(j)                                                           
      x(j)=0.0                                                          
      t=s(j)+vv(j,j)*xj                                                 
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)      
      if(x(j).eq.xj)goto 10601                                          
      del=x(j)-xj                                                       
      dlx=max(dlx,abs(del))                                             
      s=s-del*vv(:,j)                                                   
10601 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10592                                          
      goto 10591                                                        
10592 continue                                                          
      return                                                            
      end                                                               
      subroutine fatmul(it,n,vv,x,s,z,m)                                
      implicit double precision(a-h,o-z)                                
      parameter(fac=0.2)                                                
      double precision vv(n,n),x(n),s(n),z(n)                           
      integer m(n)                                                      
      l=0                                                               
      do 10611 j=1,n                                                    
      if(x(j).eq.0.0)goto 10611                                         
      l=l+1                                                             
      m(l)=j                                                            
      z(l)=x(j)                                                         
10611 continue                                                          
      continue                                                          
      if(l .le. int(fac*n))goto 10631                                   
      if(it .ne. 1)goto 10651                                           
      s=matmul(vv,x)                                                    
      goto 10661                                                        
10651 continue                                                          
      s=s-matmul(x,vv)                                                  
10661 continue                                                          
      continue                                                          
      goto 10621                                                        
10631 if(it .ne. 1)goto 10671                                           
      do 10681 j=1,n                                                    
      s(j)=dot_product(vv(j,m(1:l)),z(1:l))                             
10681 continue                                                          
      continue                                                          
      goto 10691                                                        
10671 continue                                                          
      do 10701 j=1,n                                                    
      s(j)=s(j)-dot_product(vv(m(1:l),j),z(1:l))                        
10701 continue                                                          
      continue                                                          
10691 continue                                                          
10621 continue                                                          
      return                                                            
      end                                                               
      subroutine inv(n,ww,xs,wwi)                                       
      implicit double precision(a-h,o-z)                                
      double precision ww(n,n),xs(n-1,n),wwi(n,n)                       
      nm1=n-1                                                           
      xs=-xs                                                            
      wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)))             
      wwi(2:n,1)=wwi(1,1)*xs(:,1)                                       
      wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)))           
      wwi(1:nm1,n)=wwi(n,n)*xs(:,n)                                     
      do 10711 j=2,nm1                                                  
      jm1=j-1                                                           
      jp1=j+1                                                           
      wwi(j,j)=1.0/(ww(j,j)+dot_product(xs(1:jm1,j),ww(1:jm1,j))  +dot_p
     *roduct(xs(j:nm1,j),ww(jp1:n,j)))
      wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j)                                 
      wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j)                                 
10711 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine glassopath(beta,what,jerrs,rholist,nrho,n,ss,rho,ia,itr
     *, ipen,thr,maxit,ww,wwi,niter,del,jerr)
      implicit double precision(a-h,o-z)                                
      integer nrho,n,jerrs(nrho)                                        
      double precision rholist(nrho),beta(n,n,nrho),what(n,n,nrho)      
      double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                
      is=0                                                              
      do 10721 j=1,n                                                    
      do 10731 k=1,n                                                    
      rho(j,k) = rholist(nrho)                                          
10731 continue                                                          
      continue                                                          
10721 continue                                                          
      continue                                                          
      call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jer
     *r)
      jerrs(1)=jerr                                                     
      do 10741 j=1,n                                                    
      do 10751 k=1,n                                                    
      beta(j,k,nrho)=wwi(j,k)                                           
      what(j,k,nrho)=ww(j,k)                                            
10751 continue                                                          
      continue                                                          
10741 continue                                                          
      continue                                                          
      is=1                                                              
      do 10761 i=nrho,1,-1                                              
      do 10771 j=1,n                                                    
      do 10781 k=1,n                                                    
      rho(j,k)=rholist(i)                                               
10781 continue                                                          
      continue                                                          
10771 continue                                                          
      continue                                                          
      if(itr.gt.0) call dblepr('rho=', -1, rholist(i), 1)               
      itr2=itr                                                          
      if(itr2.gt.0) itr2=itr-1                                          
      call glasso(n,ss,rho,ia,is,itr2,ipen,thr,maxit,ww,wwi,niter,del,je
     *rr)
      jerrs(i)=jerr                                                     
      do 10791 j=1,n                                                    
      do 10801 k=1,n                                                    
      beta(j,k,i)=wwi(j,k)                                              
      what(j,k,i)=ww(j,k)                                               
10801 continue                                                          
      continue                                                          
10791 continue                                                          
      continue                                                          
10761 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
