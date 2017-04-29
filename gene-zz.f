      program pgen
c----------------------------------------------------------
c      generate 2D Phosphorus structures
c      Authors: Zhen Zhu and Jie Guan
c      Date: 11 April 2014
c      use Zhen Zhu's critieria
c----------------------------------------------------------
c234567
      implicit none
C
      integer*4 iat, nat, itau, ntau, add_zat, countj
      integer*4 ntot
      integer*4 ia1, ia2, ia3, na1, na2, na3, i,k, l, m,n,j
      integer*4 sum_j
c      integer*2 atm_type(100)
      integer*4  neib(4,5000), zat(5000), z_neib, atm_count
      real*4  coor(3,1000), coor_unit(3,4), car(3,1000),lat(3)
c      real*4 d, scal12, scal3, dxy, dz
      character pho 
      write (6,200) 
 200  format (2x,'Enter number of unit cells ',
     &  'in the a1, a2 direction'/)
      read (5,*) na1, na2
      write (6,*) 'sum of neiber, 1(3) or 2?'
      read (5,*) sum_j
c     a1 is longer than a2
c       
      ntau = 4 
c     number of atom in a unit cell
      open(1,file='neib.txt')
c 
      open(2,file='position.txt')
c
      open(3,file='coord.txt')  
      
      open(4,file='xyz.ANI')
      pho='P'
      nat = ntau * na1 * na2
      lat(1)=5.0*na1
      lat(2)=3.3*na2
      lat(3)=1.0
      coor_unit(1,1)=0.0
      coor_unit(2,1)=0.0
      coor_unit(1,2)=0.16667
      coor_unit(2,2)=0.5
      coor_unit(1,3)=0.5
      coor_unit(2,3)=0.5
      coor_unit(1,4)=0.66667
      coor_unit(2,4)=0.0
      do i=1, na1
         do j=1, na2
            atm_count=4*(i-1)*na2+j
            coor(1,atm_count)=coor_unit(1,1)/na1+(i-1)/1.0/na1
            coor(2,atm_count)=coor_unit(2,1)/na2+(j-1)/1.0/na2
            coor(1,atm_count+na2)=coor_unit(1,2)/na1+(i-1)/1.0/na1
            coor(2,atm_count+na2)=coor_unit(2,2)/na2+(j-1)/1.0/na2
            coor(1,atm_count+2*na2)=coor_unit(1,3)/na1+(i-1)/1.0/na1
            coor(2,atm_count+2*na2)=coor_unit(2,3)/na2+(j-1)/1.0/na2
            coor(1,atm_count+3*na2)=coor_unit(1,4)/na1+(i-1)/1.0/na1
            coor(2,atm_count+3*na2)=coor_unit(2,4)/na2+(j-1)/1.0/na2
         end do
      end do
c
c     test composition
c      open(3,file='decom.txt')
c
c
c     generate neighbour matrix
c
      do i=1, nat
         neib(1,i)=i
         if (i-na2<=0) then
            neib(2,i)=i-na2+nat
         else
            neib(2,i)=i-na2
         end if
         if (i+na2>nat) then
            neib(3,i)=i+na2-nat
         else
            neib(3,i)=i+na2
         end if
         if (MOD((i-1)/na2,4)==0) then
             neib(4,i)=i+na2-1
             if((i-1)/na2==(neib(4,i)-1)/na2) then
                neib(4,i)=neib(4,i)+na2
             end if
             if (neib(4,i)>nat) then
                neib(4,i)=neib(4,i)-nat
             end if
         else if (MOD((i-1)/na2,4)==1) then
             neib(4,i)=i-na2+1
                if((i-1)/na2==(neib(4,i)-1)/na2) then
                   neib(4,i)=neib(4,i)-na2
                end if
             if (neib(4,i)<=0) then
                neib(4,i)=neib(4,i)+nat
             end if
         else if (MOD((i-1)/na2,4)==2) then
             neib(4,i)=i+na2+1
                if((neib(4,i)-2)/na2==i/na2) then
                   neib(4,i)=neib(4,i)-na2
                end if
             if (neib(4,i)>nat) then
                 neib(4,i)=neib(4,i)-nat
             end if
         else if (MOD((i-1)/na2,4)==3) then
             neib(4,i)=i-na2-1
                if(neib(4,i)/na2==(i-2)/na2) then
                   neib(4,i)=neib(4,i)+na2
                end if
             if (neib(4,i)<=0) then
                neib(4,i)=neib(4,i)+nat
             end if
         end if
c             if (neib(4,i)<=0) then
c                neib(4,i)=neib(4,i)+nat
c             end if
c           else
c             neib(4,i)=i-2*na2+1
c             if (neib(4,i)<=0) then
c                neib(4,i)=neib(4,i)+nat
c             end if
c           end if 
c          end if
      end do
c write the map of nearest neigbour      
      do i=1, nat
         write(1,*) (neib(j,i),j=1,4)
      end do
c end write
c
c generate the (1,0) map to the lattice
c
      ntot=1
      do i=1, nat-1
         ntot=ntot*2
      end do
      write(*,*) ntot
c
      do i=1,ntot
         zat(nat)=0
         k=i
         do j=1,nat-1
           zat(j)=MOD(k,2)
           k=k/2
c           write(*,*) k
         end do
c         do j=1,nat
c            write(*,*) j, zat(j)
c         end do
         add_zat=0
         do j=1, nat
           add_zat=add_zat+zat(j)
         end do
         if (add_zat==nat/2) then
            countj=0
c            do j=1, nat
c               write(*,*) j, zat(j)
c            end do
            do j=1, nat
            z_neib=0
            m=neib(2,j)
            n=neib(3,j)
            l=neib(4,j)
            z_neib=zat(j)+zat(m)+zat(n)+zat(l)
c            write(*,*) z_neib
            if (sum_j==2) then 
              if (z_neib==0.or.z_neib==3.or.z_neib==1.or.z_neib==4) then
                GOTO 100
              end if
            else if (sum_j==1 .or. sum_j==3) then
              if (z_neib==0 .or. z_neib==2 .or. z_neib==4) then
                GOTO 100
              end if
            end if
c            write(*,*) z_neib
               countj=j
c            write(*,*) countj 
            end do
            if (countj==nat) then
               write(2,*) nat
               write(3,*) nat
               write(4,*) nat
               write(2,*)
               write(3,*)
               write(4,*)
               do j=1, nat
                  write(2,*) j, zat(j)
                  write(3,*) (coor(k,j),k=1,2),zat(j)
               write(4,*) pho,coor(1,j)*lat(1),coor(2,j)*lat(2),zat(j)
c test neib
c
c
c                 write(4,*) (neib(k,j),k=1,4), zat(neib(1)),&
c                   zat(neib(2,j)), zat(neib(3,j)), zat(neib(4,j))
               end do
               write(2,*)
            end if
c         end if   
c     finish sum_j first
c test neib
c
c
c                 write(4,*) (neib(k,j),k=1,4), zat(neib(1)),&
c                   zat(neib(2,j)), zat(neib(3,j)), zat(neib(4,j))
c         end if
      end if
100   continue
      end do
c
c     generate the fractional coordinate
c
c      do i=1, nat
c         coor(3,i)=zat(i)
c         write(3,*) (coor(j,i),j=1,3)     
c      end do
c
c      
      close(1)
      close(2)
      close(3)
      close(4)
      stop
      end program pgen
