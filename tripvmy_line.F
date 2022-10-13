
c======================================================================
      subroutine tripvmy_line(ami,aci,api,rrr,n1i,n1f,m1)
      implicit none
      integer :: n1i,n1f,m1
      real, dimension(m1) :: ami,aci,api,rrr
      real, dimension(m1) :: q,s,fei
      real    :: fn,p
      integer :: ia,ii,i,l
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                                                                
cm      dimension ami(m1),aci(m1),api(m1),rrr(m1)
cm      dimension q(m1),s(m1),fei(m1)
      
      ia = n1i + 1
      ii = n1i + n1f 
c
c   COEFFICIENTS FOR TRIDIAGONAL INVERSION
c

c
c  THE INVERSION STARTS
c
      q(n1i) = -api(n1i)/aci(n1i) 
      s(n1i) = -ami(n1i)/aci(n1i)
      fn = rrr(n1f)
      rrr(n1i) = rrr(n1i)/aci(n1i)                   
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 i=ia,n1f                                                
        p =1./( aci(i) + ami(i)*q(i-1))
        q(i) = - api(i)*p                    
        s(i) = - ami(i)*s(i-1)*p
        rrr(i) = ( rrr(i) - ami(i)*rrr(i-1))*p        
   10 continue             
c                                                                       
c     backward pass                                                     
c                                                                           
      s(n1f) = 1.   
      fei(n1f) = 0.
               
      do 11 l=ia,n1f       
        i = ii - l         
        s(i) = s(i) + q(i)*s(i+1)        
        fei(i) = rrr(i) + q(i)*fei(i+1)                     
   11 continue                     
              
      rrr(n1f)=(fn-api(i)*fei(n1i) - 
     %      ami(i)*fei(n1f-1))/(api(i)*s(n1i) +
     %      ami(i)*s(n1f-1)+aci(i))                  
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 l=ia,n1f         
        i = ii -l               
        rrr(i) = rrr(n1f)*s(i) + fei(i)                                
   12 continue
                                   
      return                                 
      end 
c======================================================================
