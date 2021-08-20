module rmsd
implicit none
contains
  subroutine get_cov_matrix(coord_c, coord_r, coord_c_c, coord_r_c, cov_matrix, m, n)

    integer::m,n
    integer::i,j,k
    real(kind=4)::sum_value
    real(kind=4),dimension(n)::coord_c_c,coord_r_c
    real(kind=4),dimension(m,n)::coord_c,coord_r
    real(kind=4),dimension(m,n)::coord_c_d,coord_r_d
    real(kind=4),dimension(n,n)::cov_matrix

    !f2py intent(in)::m,n
    !f2py intent(in)::coord_c_c,coord_r_c
    !f2py intent(in)::coord_c,coord_r
    !f2py intent(out)::cov_matrix

    do i=1,m
      do j=1,n
        coord_c_d(i,j)=coord_c(i,j)-coord_c_c(j)
        coord_r_d(i,j)=coord_r(i,j)-coord_r_c(j)
      end do
    end do

    do i=1,n
      do j=1,n
        sum_value=0.0
        do k=1,m
          sum_value=sum_value+coord_c_d(k,i)*coord_r_d(k,j)
        end do
      cov_matrix(i,j)=sum_value
      end do
    end do

    return

  end subroutine get_cov_matrix

  subroutine quarternion_rotate(cov_matrix,quart_matrix,n)

    integer::n
    real(kind=4),dimension(n,n)::cov_matrix
    real(kind=4),dimension(4,4)::quart_matrix

    !f2py intent(in)::n
    !f2py intent(in)::cov_matrix
    !f2py intent(out)::quart_matrix

    quart_matrix(1,1)=cov_matrix(1,1)+cov_matrix(2,2)+cov_matrix(3,3)
    quart_matrix(1,2)=cov_matrix(2,3)-cov_matrix(3,2)
    quart_matrix(1,3)=cov_matrix(3,1)-cov_matrix(1,3)
    quart_matrix(1,4)=cov_matrix(1,2)-cov_matrix(2,1)
    quart_matrix(2,1)=cov_matrix(2,3)-cov_matrix(3,2)
    quart_matrix(2,2)=cov_matrix(1,1)-cov_matrix(2,2)-cov_matrix(3,3)
    quart_matrix(2,3)=cov_matrix(1,2)+cov_matrix(2,1)
    quart_matrix(2,4)=cov_matrix(1,3)+cov_matrix(3,1)
    quart_matrix(3,1)=cov_matrix(3,1)-cov_matrix(1,3)
    quart_matrix(3,2)=cov_matrix(1,2)+cov_matrix(2,1)
    quart_matrix(3,3)=-cov_matrix(1,1)+cov_matrix(2,2)-cov_matrix(3,3)
    quart_matrix(3,4)=cov_matrix(2,3)+cov_matrix(3,2)
    quart_matrix(4,1)=cov_matrix(1,2)-cov_matrix(2,1)
    quart_matrix(4,2)=cov_matrix(1,3)+cov_matrix(3,1)
    quart_matrix(4,3)=cov_matrix(2,3)+cov_matrix(3,2)
    quart_matrix(4,4)=-cov_matrix(1,1)-cov_matrix(2,2)+cov_matrix(3,3)

    return

  end subroutine quarternion_rotate

  subroutine quart_to_rot(quart_vec,rot_matrix,m)

    integer::m
    real(kind=4),dimension(m)::quart_vec
    real(kind=4),dimension(3,3)::rot_matrix

    !f2py intent(in)::m
    !f2py intent(in)::quart_vec
    !f2py intent(out)::rot_matrix

    rot_matrix(1,1)=quart_vec(1)**2+quart_vec(2)**2-quart_vec(3)**2-quart_vec(4)**2
    rot_matrix(1,2)=2.0*(quart_vec(2)*quart_vec(3)-quart_vec(1)*quart_vec(4))
    rot_matrix(1,3)=2.0*(quart_vec(2)*quart_vec(4)+quart_vec(1)*quart_vec(3))
    rot_matrix(2,1)=2.0*(quart_vec(2)*quart_vec(3)+quart_vec(1)*quart_vec(4))
    rot_matrix(2,2)=quart_vec(1)**2-quart_vec(2)**2+quart_vec(3)**2-quart_vec(4)**2
    rot_matrix(2,3)=2.0*(quart_vec(3)*quart_vec(4)-quart_vec(1)*quart_vec(2))
    rot_matrix(3,1)=2.0*(quart_vec(2)*quart_vec(4)-quart_vec(1)*quart_vec(3))
    rot_matrix(3,2)=2.0*(quart_vec(3)*quart_vec(4)+quart_vec(1)*quart_vec(2))
    rot_matrix(3,3)=quart_vec(1)**2-quart_vec(2)**2-quart_vec(3)**2+quart_vec(4)**2

    return

  end subroutine quart_to_rot

  subroutine get_rmsd(coord_c,coord_r,coord_c_c,coord_r_c,eig_max,rmsd_value,m,n)

    integer::m,n
    integer::i,j
    real(kind=4)::sum_value
    real(kind=4)::eig_max
    real(kind=4)::rmsd_value
    real(kind=4),dimension(n)::coord_c_c,coord_r_c
    real(kind=4),dimension(m,n)::coord_c,coord_r
    real(kind=4),dimension(m,n)::coord_c_d,coord_r_d

    !f2py intent(in)::m,n
    !f2py intent(in)::coord_c_c,coord_r_c
    !f2py intent(in)::coord_c,coord_r
    !f2py intent(out)::rmsd_value

    do i=1,m
      do j=1,n
        coord_c_d(i,j)=coord_c(i,j)-coord_c_c(j)
        coord_r_d(i,j)=coord_r(i,j)-coord_r_c(j)
      end do
    end do

    sum_value=0.0
    do i=1,m
      do j=1,n
        sum_value=sum_value+coord_c_d(i,j)**2+coord_r_d(i,j)**2
      end do
    end do

!    write(*,*)sum_value,eig_max
    rmsd_value=sqrt((sum_value-2*eig_max)/m)

    return

  end subroutine get_rmsd

end module rmsd
