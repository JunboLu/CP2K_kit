include 'vector.f90'

module statistic
use vector
implicit none
contains
  !This subroutine need revision, we will use end, start, each, rather than stat_num.
  subroutine ensemble_average(data_array,energy_array,temp,stat_num,average_value,m,n)

    integer::i,j
    integer::ii,jj
    integer::m,n
    integer::stat_num
    integer::min_num
    real(kind=4)::sum_value
    real(kind=4)::average_value
    real(kind=4),dimension(stat_num)::part_prob
    real(kind=4),dimension(m)::data_array
    real(kind=4),dimension(m)::energy_array
    real(kind=4),dimension(n)::temp
    real(kind=4)::kb

    !f2py intent(in)::stat_num,m,n
    !f2py intent(in)::data_array,energy_array,temp
    !f2py intent(out)::average_value

    kb=3.1668073E-6 !unit in Hartree/K
    if ( m>n ) then
      min_num = n
    elseif ( m<n ) then
      min_num = m
    elseif ( m==n ) then
      min_num = m
    end if
    do i=1,stat_num
      ii=min_num-stat_num+i
      sum_value=0.0
      do j=1,stat_num
        jj=min_num-stat_num+j
        sum_value=sum_value+exp(0.0-(energy_array(jj)/temp(jj)-energy_array(ii)/temp(ii))/kb/27.2114)
      end do
      part_prob(i)=1.0/sum_value
    end do
    average_value=0.0
    do i=1,stat_num
      ii=min_num-stat_num+i
      average_value=average_value+part_prob(i)*data_array(ii)
    end do

    return

  end subroutine

  subroutine numerical_average(data_array, average_value, sigma, m)

    integer::i
    integer::m
    real(kind=4)::s_value, sq_value
    real(kind=4)::average_value, sigma
    real(kind=4),dimension(m)::data_array

    !f2py intent(in)::m
    !f2py intent(in)::data_array
    !f2py intent(out)::average_value, sigma

    s_value = 0.0
    sq_value = 0.0
    do i=1,m
      s_value=s_value+data_array(i)
      sq_value=sq_value+data_array(i)*data_array(i)
    end do

    average_value=s_value/m
    sigma = sqrt(sq_value/m-(s_value/m)*(s_value/m))

    return

  end subroutine

  subroutine fourier_transform(data_tcf,time_array,freq_array,time_step,intensity_array,m,n)

    integer::m,n
    integer::i,j
    real(kind=4)::pi
    real(kind=4)::sum_value
    real(kind=4)::time_step
    real(kind=4)::freq_in_hz
    real(kind=4)::wave_num_to_hz
    real(kind=4),dimension(m)::data_tcf
    real(kind=4),dimension(m)::time_array
    real(kind=4),dimension(m)::time_in_s
    real(kind=4),dimension(n)::freq_array
    real(kind=4),dimension(n)::intensity_array

    !f2py intent(in)::m,n
    !f2py intent(in)::time_step
    !f2py intent(in)::freq_array
    !f2py intent(in)::time_array,data_tcf
    !f2py intent(out)::intensity_array

    pi=3.1415926
    wave_num_to_hz=29979245800.0
    do i=1,m
      time_in_s(i)=time_array(i)*1.0E-15
    end do
    do i=1,n
      freq_in_hz=freq_array(i)*wave_num_to_hz
      sum_value=0.0
      do j=1,m
        sum_value=sum_value+data_tcf(j)*cos(2*pi*freq_in_hz*time_in_s(j))*time_step
      end do
      intensity_array(i)=sum_value
    end do

    return

  end subroutine fourier_transform

  subroutine matrix_element_sum(matrix, index_1, index_2, sum_dim, matrix_sum, u, v, w)

    integer::i,j,k
    integer::u,v,w
    integer::sum_dim
    integer::index_1,index_2
    real(kind=4)::sum_value
    real(kind=4),dimension(sum_dim)::matrix_sum
    real(kind=4),dimension(u,v,w)::matrix

    !f2py intent(in)::u,v,w
    !f2py intent(in)::sum_dim
    !f2py intent(in)::index_1, index_2
    !f2py intent(in)::matrix
    !f2py intent(out)::matrix_sum

    if ( index_1 == 1 .and. index_2 == 2 ) then
      do i=1,w
        sum_value = 0.0
        do j=1,u
          do k=1,v
            sum_value = sum_value + matrix(j,k,i)
          end do
        end do
        matrix_sum(i) = sum_value/(u*v)
     end do

    elseif ( index_1 == 1 .and. index_2 == 3 ) then
      do i=1,v
        sum_value = 0.0
        do j=1,u
          do k=1,w
            sum_value = sum_value + matrix(j,i,k)
          end do
        end do
        matrix_sum(i) = sum_value/(u*w)
     end do

    elseif ( index_1 == 2 .and. index_2 == 3 ) then
      do i=1,u
        sum_value = 0.0
        do j=1,v
          do k=1,w
            sum_value = sum_value + matrix(i,j,k)
          end do
        end do
        matrix_sum(i) = sum_value/(v*w)
      end do
    end if

    return

  end subroutine matrix_element_sum

  subroutine data_mode(cluster_pos,cluster_vel,a_vec,b_vec,c_vec,mass,Q_mode_1,Q_mode_2,Q_mode_3,m,n)

    ! Reference: P. Bopp, Chem. Phys. 1986, 106, 205-212

    integer::i,j
    integer::m,n
    real(kind=4)::sign_value_a,sign_value_b,sign_value_c
    real(kind=4)::a_vec_len,b_vec_len,c_vec_len
    real(kind=4)::proj_a_vec_len,proj_b_vec_len,proj_c_vec_len
    real(kind=4),dimension(n)::v1,v2
    real(kind=4),dimension(n)::distance_xyz
    real(kind=4),dimension(n)::a_vec,b_vec,c_vec
    real(kind=4),dimension(n)::a_vec_basis,b_vec_basis,c_vec_basis
    real(kind=4),dimension(n)::proj_a_vec,proj_b_vec,proj_c_vec
    real(kind=4),dimension(n)::v1_norm,v2_norm
    real(kind=4),dimension(n)::cross_1,cross_1_norm
    real(kind=4),dimension(n)::vel_com
    real(kind=4),dimension(n)::cluster_vel_2,cluster_vel_3
    real(kind=4),dimension(n)::vel2_sub_com,vel3_sub_com
    real(kind=4),dimension(n)::vel2_per,vel3_per
    real(kind=4),dimension(n)::vel2_plane,vel3_plane
    real(kind=4),dimension(n)::vel2_proj_u1,vel3_proj_u2
    real(kind=4),dimension(n)::vel2_proj_w1,vel3_proj_w2
    real(kind=4),dimension(m)::mass
    real(kind=4),dimension(2,n)::vec
    real(kind=4),dimension(m,n)::cluster_pos
    real(kind=4),dimension(m,n)::cluster_vel
    real(kind=4),dimension(n)::Q_mode_1,Q_mode_2,Q_mode_3

    !f2py intent(in)::m,n
    !f2py intent(in)::mass
    !f2py intent(in)::a_vec,b_vec,c_vec
    !f2py intent(in)::cluster_pos,cluster_vel
    !f2py intent(out)::Q_mode_1,Q_mode_2,Q_mode_3

    call get_vec_len(a_vec, a_vec_len, n)
    call get_vec_len(b_vec, b_vec_len, n)
    call get_vec_len(c_vec, c_vec_len, n)

    call norm_vec(a_vec, a_vec_basis, n)
    call norm_vec(b_vec, b_vec_basis, n)
    call norm_vec(c_vec, c_vec_basis, n)

    do i=1,2
      do j=1,n
        distance_xyz(j) = cluster_pos(i+1,j)-cluster_pos(1,j)
      end do

      call project_vec(distance_xyz, a_vec_basis, n, proj_a_vec, sign_value_a)
      call project_vec(distance_xyz, b_vec_basis, n, proj_b_vec, sign_value_b)
      call project_vec(distance_xyz, c_vec_basis, n, proj_c_vec, sign_value_c)

      call get_vec_len(proj_a_vec, proj_a_vec_len, n)
      call get_vec_len(proj_b_vec, proj_b_vec_len, n)
      call get_vec_len(proj_c_vec, proj_c_vec_len, n)

      if ( proj_a_vec_len > a_vec_len/2 ) then
        proj_a_vec_len = a_vec_len - proj_a_vec_len
        sign_value_a = -sign_value_a
      end if
      if ( proj_b_vec_len > b_vec_len/2 ) then
        proj_b_vec_len = b_vec_len - proj_b_vec_len
        sign_value_b = -sign_value_b
      end if
      if ( proj_c_vec_len > c_vec_len/2 ) then
        proj_c_vec_len = c_vec_len - proj_c_vec_len
        sign_value_c = -sign_value_c
      end if

      distance_xyz = proj_a_vec_len*a_vec_basis*sign_value_a + &
                     proj_b_vec_len*b_vec_basis*sign_value_b + &
                     proj_c_vec_len*c_vec_basis*sign_value_c

      do j=1,n
        vec(i,j) = distance_xyz(j)
      end do
    end do

    do i=1,n
      v1(i)=vec(1,i)
      v2(i)=vec(2,i)
    end do

!    do i=1,n
!      v1(i) = cluster_pos(2,i)-cluster_pos(1,i)
!      v2(i) = cluster_pos(3,i)-cluster_pos(1,i)
!    end do
!    write(*,*)v1,v2

    call norm_vec(v1, v1_norm, n)
    call norm_vec(v2, v2_norm, n)
    call get_vec_cross(v1, v2, n, cross_1)
    call norm_vec(cross_1, cross_1_norm, n)

    call center_of_mass(mass,cluster_vel,vel_com,m,n)
   
    do i=1,n
      cluster_vel_2(i)=cluster_vel(2,i)
      cluster_vel_3(i)=cluster_vel(3,i)
    end do

    vel2_sub_com=cluster_vel_2-vel_com
    vel3_sub_com=cluster_vel_3-vel_com

    call project_vec(vel2_sub_com,cross_1_norm,n,vel2_per,sign_value_a)
    call project_vec(vel3_sub_com,cross_1_norm,n,vel3_per,sign_value_b)

    vel2_plane=vel2_sub_com-vel2_per
    vel3_plane=vel3_sub_com-vel3_per

    call project_vec(vel2_plane,v1_norm,n,vel2_proj_u1,sign_value_a)
    vel2_proj_w1=vel2_plane-vel2_proj_u1
    call project_vec(vel3_plane,v2_norm,n,vel3_proj_u2,sign_value_b)
    vel3_proj_w2=vel3_plane-vel3_proj_u2

!    Q_mode_1 = vel2_proj_u1*sign_value_a+vel3_proj_u2*sign_value_b
    Q_mode_1 = vel2_proj_u1+vel3_proj_u2
    Q_mode_2 = vel2_proj_w1+vel3_proj_w2
!    Q_mode_3 = vel2_proj_u1*sign_value_a-vel3_proj_u2*sign_value_b
    Q_mode_3 = vel2_proj_u1-vel3_proj_u2

    return

  end subroutine data_mode 

end module statistic
