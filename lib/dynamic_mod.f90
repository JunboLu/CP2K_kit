include 'vector.f90'

module dynamic
use vector
implicit none
contains

  subroutine remove_coord_com(coord,atom_mass,coord_sub_com,l,m,n)

    integer::l,m,n
    integer::i,j,k
    real(kind=4),dimension(n)::coord_com
    real(kind=4),dimension(m)::atom_mass
    real(kind=4),dimension(m,n)::coord_temp
    real(kind=4),dimension(l,m,n)::coord
    real(kind=4),dimension(l,m,n)::coord_sub_com

    !f2py intent(in)::l,m,n
    !f2py intent(in)::atom_mass
    !f2py intent(in)::coord
    !f2py intent(out)::coord_sub_com

    do i=1,l
      do j=1,m
        do k=1,n
          coord_temp(j,k) = coord(i,j,k)
        end do
      end do
      call center_of_mass(atom_mass, coord_temp, coord_com, m, n)
      do j=1,m
        do k=1,n
          coord_sub_com(i,j,k) = coord(i,j,k) - coord_com(k)
        end do
      end do
    end do

    return

  end subroutine

  subroutine diffusion_einstein_sum(coord,u,sum_array,l,m,n)

    integer::l,m,n,u
    integer::i,j,k
    real(kind=4)::sum_value_1,sum_value_2
    real(kind=4),dimension(u)::sum_array_x
    real(kind=4),dimension(u)::sum_array_y
    real(kind=4),dimension(u)::sum_array_z
    real(kind=4),dimension(u)::sum_array
    real(kind=4),dimension(l,m,n)::coord

    !f2py intent(in)::l,m,n,u
    !f2py intent(in)::coord
    !f2py intent(out)::sum_array

    do i=0,u-1
      sum_value_1=0.0
      do j=1,m
        sum_value_2=0.0
        do k=1,(l-i)
          sum_value_2=sum_value_2+(coord(k,j,1)-coord(k+i,j,1))**2
        end do
        sum_value_1=sum_value_1+sum_value_2/(l-i)
      end do
      sum_array_x(i+1)=sum_value_1/m
    end do

    do i=0,u-1
      sum_value_1=0.0
      do j=1,m
        sum_value_2=0.0
        do k=1,(l-i)
          sum_value_2=sum_value_2+(coord(k,j,2)-coord(k+i,j,2))**2
        end do
        sum_value_1=sum_value_1+sum_value_2/(l-i)
      end do
      sum_array_y(i+1)=sum_value_1/m
    end do

    do i=0,u-1
      sum_value_1=0.0
      do j=1,m
        sum_value_2=0.0
        do k=1,(l-i)
          sum_value_2=sum_value_2+(coord(k,j,3)-coord(k+i,j,3))**2
        end do
        sum_value_1=sum_value_1+sum_value_2/(l-i)
      end do
      sum_array_z(i+1)=sum_value_1/m
    end do

    sum_array=sum_array_x+sum_array_y+sum_array_z

    return

  end subroutine diffusion_einstein_sum

  subroutine time_correlation(data_array,data_tcf_norm,u,l,m,n,normalize)

    integer::l,m,n,u
    integer::i,j,k
    integer::normalize
    real(kind=4)::sum_nume,sum_deno
    real(kind=4)::sum_nume_x,sum_nume_y,sum_nume_z
    real(kind=4)::sum_deno_x,sum_deno_y,sum_deno_z
    real(kind=4),dimension(u)::data_tcf_norm
    real(kind=4),dimension(l,m,n)::data_array
    real(kind=4),dimension(l,m,n)::new_data_array

    !f2py intent(in)::l,m,n,u
    !f2py intent(in)::data_array
    !f2py intent(out)::data_tcf_norm
    !f2py intent(in)::normalize

    !The velocity in CP2K is bohr/fs
    do i=1,l
      do j=1,m
        do k=1,n
          if (normalize == 1) then
            new_data_array(i,j,k)=data_array(i,j,k)
          else
            new_data_array(i,j,k)=data_array(i,j,k)*100.0*0.5291772489940979*(1.0E-10)/(2.4188843265857*1.0E-17)
          end if
        end do
      end do
    end do

    do i=1,l
      do j=1,m
        do k=1,n
          if (isnan(new_data_array(i,j,k))) then
            write(*,*)i,j,k,new_data_array(i,j,k)
          end if
        end do
      end do
    end do

    do i=0,u-1
      sum_nume=0.0
      sum_deno=0.0
      do j=1,m
        sum_nume_x=0.0
        sum_nume_y=0.0
        sum_nume_z=0.0
        sum_deno_x=0.0
        sum_deno_y=0.0
        sum_deno_z=0.0
        do k=1,(l-i)
          sum_nume_x=sum_nume_x+new_data_array(k,j,1)*new_data_array(k+i,j,1)
          sum_nume_y=sum_nume_y+new_data_array(k,j,2)*new_data_array(k+i,j,2)
          sum_nume_z=sum_nume_z+new_data_array(k,j,3)*new_data_array(k+i,j,3)
          sum_deno_x=sum_deno_x+new_data_array(k,j,1)*new_data_array(k,j,1)
          sum_deno_y=sum_deno_y+new_data_array(k,j,2)*new_data_array(k,j,2)
          sum_deno_z=sum_deno_z+new_data_array(k,j,3)*new_data_array(k,j,3)
        end do
        sum_nume=sum_nume+sum_nume_x+sum_nume_y+sum_nume_z
        sum_deno=sum_deno+sum_deno_x+sum_deno_y+sum_deno_z
      end do
      if (normalize == 1) then
        data_tcf_norm(i+1)=sum_nume/sum_deno
      else
        data_tcf_norm(i+1)=sum_nume/((l-i)*m)
      end if
    end do

    return

  end subroutine time_correlation 

end module dynamic
