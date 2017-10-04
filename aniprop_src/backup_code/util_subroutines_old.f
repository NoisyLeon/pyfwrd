      subroutine cbcktr(n1,z,zi,dr,di,er,ei,ip)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension z(n1,1),ip(1),er(1),ei(1),dr(1),di(1),zi(n1,1)
c  back transform with unit lower triangular matrix
      do 300 i=1,n1
      dr(i)=er(ip(i))
  300 di(i)=ei(ip(i))
      if(n1.eq.1) go to 400
      do 310 i=2,n1
      i1=i-1
      do 310 j=1,i1
      zkk1=z(ip(i),j)
      zkk2=zi(ip(i),j)
      dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
  310 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
  400 continue
c  back transform with upper triangular matrix
      do 320 ii=1,n1
      i=n1+1-ii
      ip1=i+1
      uii=z(ip(i),i)**2+zi(ip(i),i)**2
      if(i.eq.n1) go to 340
      do 330 j=ip1,n1
      zkk1=z(ip(i),j)
      zkk2=zi(ip(i),j)
      dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
  330 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
  340 dri=dr(i)
      dii=di(i)
      zkk1=z(ip(i),i)
      zkk2=zi(ip(i),i)
      dr(i)=(dri*zkk1+dii*zkk2)/uii
  320 di(i)=(dii*zkk1-dri*zkk2)/uii
      return
      end
 
      subroutine bcktr(n1,z,dr,er,ip)
c  performs backtransform on input vector er - 'y'
c  to find solution dr - 'x'
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension z(n1,1),ip(1),er(1),dr(1)
c  back transform with unit lower triangular matrix
      do 300 i=1,n1
  300 dr(i)=er(ip(i))
      if(n1.eq.1) go to 400
      do 310 i=2,n1
      i1=i-1
      do 310 j=1,i1
  310 dr(i)=dr(i)-z(ip(i),j)*dr(j)
  400 continue
c  back transform with upper triangular matrix
      do 320 ii=1,n1
      i=n1+1-ii
      ip1=i+1
      if(i.eq.n1) go to 320
      do 330 j=ip1,n1
  330 dr(i)=dr(i)-z(ip(i),j)*dr(j)
  320 dr(i)=dr(i)/z(ip(i),i)
      return
      end

      subroutine lup(n,a,ip)
c  finds lu decomp of a using partial pivoting         
c  output in c (upper triangle/unit lower triangle) and                      
c  pivoting sequence returned in ip(n)                                          
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension a(n,1),ip(1)                             
      tol=1.d-14                                                                
      do 50 i=1,n
   50 ip(i)=i                                                                   
      nm1=n-1                                                                   
      if(n.eq.1) go to 700                                                      
      do 100 i=1,nm1                                                            
      aam=0.d0                                                                  
      do 200 j=i,n                                                              
      aai=a(ip(j),i)**2                                          
      if(aam.gt.aai) go to 200                                                  
      aam=aai                                                                    
      ipi=ip(j)                                                                 
      jm=j                                                                      
  200 continue                                                                  
      if(aam.lt.tol) go to 400                                                  
      ip(jm)=ip(i)                                                              
      ip(i)=ipi                                                                 
      i1=i+1                                                                    
      do 100 j=i1,n                                                             
      ipj=ip(j)                                                                 
c  if victim index is already zero, dont bother to rub it out                   
      tem=dabs(a(ipj,i))
      if(tem.lt.tol) go to 100                                                  
      b=(a(ipj,i)*a(ipi,i))/aam                             
      a(ipj,i)=b                                                                
      do 500 k=i1,n                                                             
      a(ipj,k)=a(ipj,k)-b*a(ipi,k)
  500 continue
  100 continue                                                                  
  700 continue              
      return                                                                    
c  400 print 101,aam,i 
  400 nooutput=1
  101 format('near-zero pivot ',e12.5,'  on column',i3)                         
      stop
      end               

        subroutine clup(n,a,ai,ip)                                    
c  finds lu decomp of a+i*ai using partial pivoting 
c  pivoting sequence returned in ip(n)
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension a(n,1),ai(n,1),ip(1)                             
      tol=1.d-14                                                                
c  initialize permutation vector                                                
      do 50 i=1,n
   50 ip(i)=i                                                                   
      nm1=n-1                                                                   
      if(n.eq.1) go to 700                                                      
      do 100 i=1,nm1                                                            
      aam=0.d0                                                                  
      do 200 j=i,n                                                              
      aai=a(ip(j),i)**2+ai(ip(j),i)**2                                          
      if(aam.gt.aai) go to 200                                                  
      aam=aai                                                                   
      ipi=ip(j)                                                                 
      jm=j                                                                      
  200 continue                                                                  
      if(aam.lt.tol) go to 400                                                  
      ip(jm)=ip(i)                                                              
      ip(i)=ipi                                                                 
      i1=i+1                                                                    
      do 100 j=i1,n                                                             
      ipj=ip(j)                                                                 
c  if victim index is already zero, dont bother to rub it out                   
      tem=dabs(a(ipj,i))+dabs(ai(ipj,i))                                        
      if(tem.lt.tol) go to 100                                                  
      b=(a(ipj,i)*a(ipi,i)+ai(ipj,i)*ai(ipi,i))/aam                             
      bi=(ai(ipj,i)*a(ipi,i)-a(ipj,i)*ai(ipi,i))/aam                            
      a(ipj,i)=b                                                                
      ai(ipj,i)=bi                                                              
      do 500 k=i1,n                                                             
      a(ipj,k)=a(ipj,k)-b*a(ipi,k)+bi*ai(ipi,k)                                 
  500 ai(ipj,k)=ai(ipj,k)-b*ai(ipi,k)-bi*a(ipi,k)                               
  100 continue                                                                  
  700 continue                                                                  
      return                                                                    
c  400 print 101,aam,i    
  400 nooutput=1
  101 format('near-zero pivot ',e12.5,'  on column',i3)                         
      stop                                                                 
      end      

      subroutine solve(nn,a,x,y)
c  solves the nxn system of equations a*x=y using gaussian elimination 
c  and partial pivoting
c  if n<0 the lu decomposition is already done
c  note that the matrix a is modified
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      common/solve_/ip(1000)
      dimension a(1),x(1),y(1)
      n=nn
      if(n.gt.0)then
        call lup(n,a,ip)
      else
        n=-n
      endif
c      print 101,((i,j,a(i,j),j=1,n),i=1,n)
c  101 format(' a(',2i2,')=',e15.5)
c      type 102,(ip(i),i=1,n)
c  102 format(10i5)
      call bcktr(n,a,x,y,ip)
      return
      end
                                                        
      subroutine csolve(nn,a,ai,x,xi,y,yi)
c  solves the complex nxn system of equations a*x=y using gaussian elimination 
c  and partial pivoting
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension ip(1000),a(1),x(1),y(1),ai(1),xi(1),yi(1)
      n=nn
      if(n.gt.0) then 
        call clup(n,a,ai,ip)
      else
        n=-n
      endif
      call cbcktr(n,a,ai,x,xi,y,yi,ip)
      return
      end
 
      function dsind(dgr_argument)
      real*8 dsind,dgr_argument
      real*8,parameter:: quadpi=3.1415926 
      real*8,parameter:: dgr_to_rad=(quadpi/180)

      dsind=sin(dgr_to_rad*dgr_argument)
      end function 

      function dcosd(dgr_argument)
      real*8 dcosd,dgr_argument
      real*8,parameter:: quadpi=3.1415926 
      real*8,parameter:: dgr_to_rad=(quadpi/180)

      dcosd=cos(dgr_to_rad*dgr_argument)
      end function 


