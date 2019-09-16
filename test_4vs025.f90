program test
   use omp_lib
   implicit none

   integer :: i,k,n,l
   real(8) :: x(2),y,t0,t1,dt,t4,t025

   n = int(1e4)
   do k=1,14
     dt = 0.
     do l=1,10
       t0 = omp_get_wtime()
       do i = 1,n
         call random_number(x)
         y = (x(1)*x(1) + x(1) - 2.)/4.*x(2)
       enddo
       t1 = omp_get_wtime()
       t4 = t1 - t0

       t0 = omp_get_wtime()
       do i = 1,n
         call random_number(x)
         y = (x(1)*x(1) + x(1) - 2.)*0.25*x(2)
       enddo
       t1 = omp_get_wtime()
       t025 = t1 - t0

       dt = dt + (t4 - t025)
       print*,'---',t4,dt/10.
     enddo
     n = n*2
     print*,'==='
   enddo
end program test
