module vector
implicit none
contains
  subroutine get_vec_cross(v1, v2, m, cross)

    integer,intent(in)::m
    real(kind=4),intent(in),dimension(m)::v1,v2
    real(kind=4),intent(out),dimension(m)::cross

    cross(1)=v1(2)*v2(3)-v1(3)*v2(2)
    cross(2)=v1(3)*v2(1)-v1(1)*v2(3)
    cross(3)=v1(1)*v2(2)-v1(2)*v2(1)

    return

  end subroutine get_vec_cross

  subroutine get_vec_len(vec,vec_len,m)

    integer::i
    integer,intent(in)::m
    real(kind=4)::sum_value
    real(kind=4),intent(in),dimension(m)::vec
    real(kind=4),intent(out)::vec_len

    sum_value=0.0
    do i=1,m
      sum_value=sum_value+vec(i)**2
    end do
    vec_len = sqrt(sum_value)

    return

  end subroutine get_vec_len

  subroutine norm_vec(vec,new_vec,m)

    integer::i
    integer,intent(in)::m
    real(kind=4)::sum_value
    real(kind=4),intent(in),dimension(m)::vec
    real(kind=4),intent(out),dimension(m)::new_vec

    sum_value = 0.0
    do i=1,m
      sum_value=sum_value+vec(i)**2
    end do
    do i=1,m
      new_vec(i)=vec(i)/sqrt(sum_value)
    end do

    return

  end subroutine norm_vec

  subroutine project_vec(vec_1,vec_2,m,vec_1_proj_vec_2,sign_value)
  ! project vec_1 on the basis of vec_2

    integer,intent(in)::m
    real(kind=4)::a,b
    real(kind=4),intent(out)::sign_value
    real(kind=4),intent(in),dimension(m)::vec_1,vec_2
    real(kind=4),intent(out),dimension(m)::vec_1_proj_vec_2

    a=DOT_PRODUCT(vec_1,vec_2)
    b=DOT_PRODUCT(vec_2,vec_2)
    vec_1_proj_vec_2=a*vec_2/b

    if ( a<0.0 ) then
      sign_value=-1.0
    else
      sign_value=1.0
    end if
 
    return

  end subroutine project_vec

  subroutine center_of_mass(mass,data_array,data_center,n,m)

    integer::i,j
    integer,intent(in)::n,m
    real(kind=4)::sum_value
    real(kind=4)::sum_value_mass
    real(kind=4),intent(in),dimension(n)::mass
    real(kind=4),intent(in),dimension(n,m)::data_array
    real(kind=4),intent(out),dimension(m)::data_center

    sum_value_mass=0.0

    do i=1,n
      sum_value_mass=sum_value_mass+mass(i)
    end do

    do i=1,m
      sum_value=0.0
      do j=1,n
        sum_value=sum_value+mass(j)*data_array(j,i)
      end do
      data_center(i)=sum_value/sum_value_mass
    end do

    return

  end subroutine center_of_mass

end module vector
