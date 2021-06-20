include 'vector.f90'

module geometry
use vector
contains
  subroutine calculate_distance(coord_1, coord_2, a_vec, b_vec, c_vec, calc_distance, m, n)
  !a_vec, b_vec, c_vec are triclinic cell vectors

    integer::m, n
    integer::i, j
    real(kind=4)::a_vec_len, b_vec_len, c_vec_len
    real(kind=4)::proj_a_vec_len, proj_b_vec_len, proj_c_vec_len
    real(kind=4),dimension(n)::distance_xyz
    real(kind=4),dimension(n)::a_vec, b_vec, c_vec
    real(kind=4),dimension(n)::a_vec_basis, b_vec_basis, c_vec_basis
    real(kind=4),dimension(n)::proj_a_vec, proj_b_vec, proj_c_vec
    real(kind=4),dimension(m)::calc_distance
    real(kind=4),dimension(m,n)::coord_1, coord_2

    !f2py intent(in)::m, n
    !f2py intent(in)::a_vec, b_vec, c_vec
    !f2py intent(in)::coord_1, coord_2
    !f2py intent(out)::calc_distance

    i = 0
    j = 0

    a_vec_len = 0.0
    b_vec_len = 0.0
    c_vec_len = 0.0
    proj_a_vec_len = 0.0
    proj_b_vec_len = 0.0
    proj_c_vec_len = 0.0

    do i=1,n
      a_vec_basis(i) = 0.0
      b_vec_basis(i) = 0.0
      c_vec_basis(i) = 0.0
      proj_a_vec(i) = 0.0
      proj_b_vec(i) = 0.0
      proj_c_vec(i) = 0.0
    end do

    do i=1,m
      calc_distance(i) = 0.0
    end do

    call get_vec_len(a_vec, a_vec_len, n)
    call get_vec_len(b_vec, b_vec_len, n)
    call get_vec_len(c_vec, c_vec_len, n)

    call norm_vec(a_vec, a_vec_basis, n)
    call norm_vec(b_vec, b_vec_basis, n)
    call norm_vec(c_vec, c_vec_basis, n)

    do i=1,m
      do j=1,n
        distance_xyz(j) = coord_1(i,j) - coord_2(i,j)
      end do

      call project_vec(distance_xyz, a_vec_basis, n, proj_a_vec)
      call project_vec(distance_xyz, b_vec_basis, n, proj_b_vec)
      call project_vec(distance_xyz, c_vec_basis, n, proj_c_vec)

      call get_vec_len(proj_a_vec, proj_a_vec_len, n)
      call get_vec_len(proj_b_vec, proj_b_vec_len, n)
      call get_vec_len(proj_c_vec, proj_c_vec_len, n)

      if ( proj_a_vec_len > a_vec_len/2 ) then
        proj_a_vec_len = a_vec_len - proj_a_vec_len
      end if
      if ( proj_b_vec_len > b_vec_len/2 ) then
        proj_b_vec_len = b_vec_len - proj_b_vec_len
      end if
      if ( proj_c_vec_len > c_vec_len/2 ) then
        proj_c_vec_len = c_vec_len - proj_c_vec_len
      end if

      distance_xyz = proj_a_vec_len*a_vec_basis + &
                     proj_b_vec_len*b_vec_basis + &
                     proj_c_vec_len*c_vec_basis

      call get_vec_len(distance_xyz, calc_distance(i), n)

    end do

    return

  end subroutine calculate_distance

  subroutine calculate_angle(coord_1,coord_2,coord_3,calc_angle,m,n)

    integer::m,n
    integer::i,j
    real(kind=4)::a,b,c
    real(kind=4),dimension(m)::calc_angle
    real(kind=4),dimension(n)::vector_1,vector_2
    real(kind=4),dimension(m,n)::coord_1,coord_2,coord_3

    !f2py intent(in)::m,n
    !f2py intent(in)::coord_1,coord_2,coord_3
    !f2py intent(out)::calc_angle

    do i=1,m
      do j=1,n
        vector_1(j)=coord_1(i,j)-coord_2(i,j)
        vector_2(j)=coord_3(i,j)-coord_2(i,j)
      end do
      a=vector_1(1)*vector_2(1)+vector_1(2)*vector_2(2)+vector_1(3)*vector_2(3)
      b=sqrt(vector_1(1)*vector_1(1)+vector_1(2)*vector_1(2)+vector_1(3)*vector_1(3))
      c=sqrt(vector_2(1)*vector_2(1)+vector_2(2)*vector_2(2)+vector_2(3)*vector_2(3))
      calc_angle(i)=acos(a/(b*c))
    end do

    return

  end subroutine calculate_angle

  subroutine trans_box_center(coord, atoms_mass, a_vec, b_vec, c_vec, new_coord, m, n)

    integer::m, n
    integer::i, j
    real(kind=4)::shift_a_len, shift_b_len, shift_c_len
    real(kind=4)::proj_a_vec_len, proj_b_vec_len, proj_c_vec_len
    real(kind=4),dimension(n)::coord_temp
    real(kind=4),dimension(n)::box_center_vec
    real(kind=4),dimension(n)::coord_center_mass
    real(kind=4),dimension(n)::a_vec, b_vec, c_vec
    real(kind=4),dimension(n)::a_vec_basis, b_vec_basis, c_vec_basis
    real(kind=4),dimension(n)::proj_a_vec, proj_b_vec, proj_c_vec
    real(kind=4),dimension(n)::shift
    real(kind=4),dimension(m)::atoms_mass
    real(kind=4),dimension(m,n)::coord
    real(kind=4),dimension(m,n)::new_coord

    !f2py intent(in)::m, n
    !f2py intent(in)::coord
    !f2py intent(in)::a_vec, b_vec, c_vec
    !f2py intent(in)::atoms_mass
    !f2py intent(out)::new_coord

    i = 0
    j = 0

    proj_a_vec_len = 0.0
    proj_b_vec_len = 0.0
    proj_c_vec_len = 0.0
    shift_a_len = 0.0
    shift_b_len = 0.0
    shift_c_len = 0.0

    do i=1,n
      box_center_vec(i) = 0.0
      coord_center_mass(i) = 0.0
      a_vec_basis(i) = 0.0
      b_vec_basis(i) = 0.0
      c_vec_basis(i) = 0.0
      proj_a_vec(i) = 0.0
      proj_b_vec(i) = 0.0
      proj_c_vec(i) = 0.0
    end do

    call get_vec_len(a_vec, a_vec_len, n)
    call get_vec_len(b_vec, b_vec_len, n)
    call get_vec_len(c_vec, c_vec_len, n)

    call norm_vec(a_vec, a_vec_basis, n)
    call norm_vec(b_vec, b_vec_basis, n)
    call norm_vec(c_vec, c_vec_basis, n)

    box_center_vec = 0.5*(a_vec + b_vec + c_vec)
    call center_of_mass(atoms_mass, coord, coord_center_mass, m, n) 
    shift = box_center_vec - coord_center_mass

    call project_vec(shift, a_vec_basis, n, proj_a_vec)
    call project_vec(shift, b_vec_basis, n, proj_b_vec)
    call project_vec(shift, c_vec_basis, n, proj_c_vec)

    call get_vec_len(proj_a_vec, proj_a_vec_len, n)
    call get_vec_len(proj_b_vec, proj_b_vec_len, n)
    call get_vec_len(proj_c_vec, proj_c_vec_len, n)

    if ( DOT_PRODUCT(a_vec_basis,shift) > 0 ) then
      shift_a_len = proj_a_vec_len 
    else
      shift_a_len = -proj_a_vec_len
    end if
    if ( DOT_PRODUCT(b_vec_basis,shift) > 0 ) then
      shift_b_len = proj_b_vec_len 
    else
      shift_b_len = -proj_b_vec_len
    end if
    if ( DOT_PRODUCT(c_vec_basis,shift) > 0 ) then
      shift_c_len = proj_c_vec_len
    else
      shift_c_len = -proj_c_vec_len
    end if

    do i=1,m
      do j=1,n
        coord_temp(j) = coord(i,j)
      end do
      coord_temp = coord_temp + shift_a_len*a_vec_basis + &
                   shift_b_len*b_vec_basis + shift_c_len*c_vec_basis
      do j=1,n
        new_coord(i,j) = coord_temp(j)
      end do
    end do

    return

  end subroutine trans_box_center

  subroutine periodic_center_box(coord, a_vec, b_vec, c_vec, trans_type, exclude_group_id, &
                                 group_atom_1_id, group_atoms_mass, new_coord, m, n, l, x, y, z)

    integer::m, n, l, x, y, z
    integer::i, j, k, q, a
    integer::trans_type
    real(kind=4)::shift_a_len, shift_b_len, shift_c_len
    real(kind=4)::a_vec_len, b_vec_len, c_vec_len
    real(kind=4)::proj_a_vec_len, proj_b_vec_len, proj_c_vec_len
    real(kind=4),dimension(n)::coord_temp, coord_center_mass
    real(kind=4),dimension(n)::a_vec, b_vec, c_vec
    real(kind=4),dimension(n)::a_vec_basis, b_vec_basis, c_vec_basis
    real(kind=4),dimension(n)::proj_a_vec, proj_b_vec, proj_c_vec
    integer,dimension(l)::exclude_group_id
    real(kind=4),dimension(x,z)::group_atoms_mass
    integer,dimension(x,y)::group_atom_1_id
    real(kind=4),dimension(m,n)::coord
    real(kind=4),dimension(m,n)::new_coord

    !f2py intent(in)::trans_type
    !f2py intent(in)::m, n, l, x, y
    !f2py intent(in)::coord
    !f2py intent(in)::a_vec, b_vec, c_vec
    !f2py intent(in)::exclude_group_id
    !f2py intent(in)::group_atom_1_id
    !f2py intent(in)::group_atoms_mass
    !f2py intent(out)::new_coord

    i = 0
    j = 0
    k = 0
    q = 0
    a = 0

    shift_a_len = 0.0
    shift_b_len = 0.0
    shift_c_len = 0.0
    proj_a_vec_len = 0.0
    proj_b_vec_len = 0.0
    proj_c_vec_len = 0.0

    do i=1,n
      coord_temp(i) = 0.0
      coord_center_mass(i) = 0.0
      a_vec_basis(i) = 0.0
      b_vec_basis(i) = 0.0
      c_vec_basis(i) = 0.0
      proj_a_vec(i) = 0.0
      proj_b_vec(i) = 0.0
      proj_c_vec(i) = 0.0
    end do

    call get_vec_len(a_vec, a_vec_len, n)
    call get_vec_len(b_vec, b_vec_len, n)
    call get_vec_len(c_vec, c_vec_len, n)

    call norm_vec(a_vec, a_vec_basis, n)
    call norm_vec(b_vec, b_vec_basis, n)
    call norm_vec(c_vec, c_vec_basis, n)

    do i=1,m
      do j=1,n
        new_coord(i,j) = 0.0
      end do
    end do

    if ( l /= 0 ) then
      do i=1,l
        a=exclude_group_id(i)
        do j=1,n
          coord_temp(j) = coord(a,j)
        end do

        call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
        call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
        call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

        call get_vec_len(proj_a_vec, proj_a_vec_len, n)
        call get_vec_len(proj_b_vec, proj_b_vec_len, n)
        call get_vec_len(proj_c_vec, proj_c_vec_len, n)

        if ( DOT_PRODUCT(a_vec_basis,coord_temp) > 0 ) then
          shift_a_len = a_vec_len*anint(proj_a_vec_len/a_vec_len)
        else
          shift_a_len = a_vec_len*anint(-proj_a_vec_len/a_vec_len)
        end if
        if ( DOT_PRODUCT(b_vec_basis,coord_temp) > 0 ) then
          shift_b_len = b_vec_len*anint(proj_b_vec_len/b_vec_len)
        else
          shift_b_len = b_vec_len*anint(-proj_b_vec_len/b_vec_len)
        end if
        if ( DOT_PRODUCT(c_vec_basis,coord_temp) > 0 ) then
          shift_c_len = c_vec_len*anint(proj_c_vec_len/c_vec_len)
        else
          shift_c_len = c_vec_len*anint(-proj_c_vec_len/c_vec_len)
        end if

        coord_temp = coord_temp - shift_a_len*a_vec_basis - &
                     shift_b_len*b_vec_basis - shift_c_len*c_vec_basis

        do j=1,n
          new_coord(a,j) = coord_temp(j)
        end do

      end do
    end if

    if (trans_type == 0) then !center the system based on the center of mass of group
      do i =1,x
        do j=1,y
          a=group_atom_1_id(i,j)
          call center_of_mass(group_atoms_mass(i,:), coord(a:(a+z-1),:), coord_center_mass, z, n)
          
          call project_vec(coord_center_mass, a_vec_basis, n, proj_a_vec)
          call project_vec(coord_center_mass, b_vec_basis, n, proj_b_vec)
          call project_vec(coord_center_mass, c_vec_basis, n, proj_c_vec)

          call get_vec_len(proj_a_vec, proj_a_vec_len, n)
          call get_vec_len(proj_b_vec, proj_b_vec_len, n)
          call get_vec_len(proj_c_vec, proj_c_vec_len, n)

          if ( DOT_PRODUCT(a_vec_basis,coord_center_mass) > 0 ) then
            shift_a_len = a_vec_len*anint(proj_a_vec_len/a_vec_len)
          else 
            shift_a_len = a_vec_len*anint(-proj_a_vec_len/a_vec_len)
          end if
          if ( DOT_PRODUCT(b_vec_basis,coord_center_mass) > 0 ) then
            shift_b_len = b_vec_len*anint(proj_b_vec_len/b_vec_len)
          else
            shift_b_len = b_vec_len*anint(-proj_b_vec_len/b_vec_len)
          end if
          if ( DOT_PRODUCT(c_vec_basis,coord_center_mass) > 0 ) then
            shift_c_len = c_vec_len*anint(proj_c_vec_len/c_vec_len)
          else
            shift_c_len = c_vec_len*anint(-proj_c_vec_len/c_vec_len)
          end if

          do k=0,z-1
            do q=1,n
              coord_temp(q) = coord(a+k,q)
            end do

            coord_temp = coord_temp - shift_a_len*a_vec_basis - &
                         shift_b_len*b_vec_basis - shift_c_len*c_vec_basis
            do q=1,n
              new_coord(a+k,q)=coord_temp(q)
            end do
          end do
        end do
      end do

    elseif (trans_type == 1) then
      do i=1,m
        do j=1,n
          coord_temp(j) = coord(i,j)
        end do
        call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
        call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
        call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

        call get_vec_len(proj_a_vec, proj_a_vec_len, n)
        call get_vec_len(proj_b_vec, proj_b_vec_len, n)
        call get_vec_len(proj_c_vec, proj_c_vec_len, n)

        if ( DOT_PRODUCT(a_vec_basis,coord_temp) > 0 ) then
          shift_a_len = a_vec_len*anint(proj_a_vec_len/a_vec_len)
        else
          shift_a_len = a_vec_len*anint(-proj_a_vec_len/a_vec_len)
        end if
        if ( DOT_PRODUCT(b_vec_basis,coord_temp) > 0 ) then
          shift_b_len = b_vec_len*anint(proj_b_vec_len/b_vec_len)
        else
          shift_b_len = b_vec_len*anint(-proj_b_vec_len/b_vec_len)
        end if
        if ( DOT_PRODUCT(c_vec_basis,coord_temp) > 0 ) then
          shift_c_len = c_vec_len*anint(proj_c_vec_len/c_vec_len)
        else
          shift_c_len = c_vec_len*anint(-proj_c_vec_len/c_vec_len)
        end if
        coord_temp = coord_temp - shift_a_len*a_vec_basis - &
                     shift_b_len*b_vec_basis - shift_c_len*c_vec_basis

        do j=1,n
          new_coord(i,j) = coord_temp(j)
        end do
      end do
    endif
 
    return

  end subroutine

  subroutine periodic_center_image(coord, a_vec, b_vec, c_vec, center_coord, trans_type, exclude_group_id, &
                                   group_atom_1_id, group_atoms_mass, new_coord, m, n, l, x, y, z)

    integer::m, n, l, x, y, z
    integer::i, j, k, q, a
    integer::trans_type
    real(kind=4)::a_vec_len, b_vec_len, c_vec_len
    real(kind=4)::shift_a_len, shift_b_len, shift_c_len
    real(kind=4)::proj_a_vec_len, proj_b_vec_len, proj_c_vec_len
    integer,dimension(l)::exclude_group_id
    integer,dimension(x,y)::group_atom_1_id
    real(kind=4),dimension(x,z)::group_atoms_mass
    real(kind=4),dimension(n)::center_coord, group_center_coord
    real(kind=4),dimension(n)::coord_temp
    real(kind=4),dimension(n)::a_vec, b_vec, c_vec
    real(kind=4),dimension(n)::a_vec_basis, b_vec_basis, c_vec_basis
    real(kind=4),dimension(n)::proj_a_vec, proj_b_vec, proj_c_vec
    real(kind=4),dimension(m,n)::coord
    real(kind=4),dimension(m,n)::new_coord

    !f2py intent(in)::trans_type
    !f2py intent(in)::m, n, l, x, y, z
    !f2py intent(in)::center_coord
    !f2py intent(in)::coord
    !f2py intent(in)::exclude_group_id
    !f2py intent(in)::group_atom_1_id
    !f2py intent(in)::group_atoms_mass
    !f2py intent(in)::a_vec, b_vec, c_vec
    !f2py intent(out)::new_coord

    i = 0
    j = 0
    k = 0
    q = 0
    a = 0
    
    a_vec_len = 0.0
    b_vec_len = 0.0
    c_vec_len = 0.0
    shift_a_len = 0.0
    shift_b_len = 0.0
    shift_c_len = 0.0
    proj_a_vec_len = 0.0
    proj_b_vec_len = 0.0
    proj_c_vec_len = 0.0

    do i=1,n
      group_center_coord(i) = 0.0
      a_vec_basis(i) = 0.0
      b_vec_basis(i) = 0.0
      c_vec_basis(i) = 0.0
      proj_a_vec(i) = 0.0
      proj_b_vec(i) = 0.0
      proj_c_vec(i) = 0.0
    end do

    do i=1,m
      do j=1,n
        new_coord(i,j) = 0.0
      end do
    end do

    call get_vec_len(a_vec, a_vec_len, n)
    call get_vec_len(b_vec, b_vec_len, n)
    call get_vec_len(c_vec, c_vec_len, n)

    call norm_vec(a_vec, a_vec_basis, n)
    call norm_vec(b_vec, b_vec_basis, n)
    call norm_vec(c_vec, c_vec_basis, n)

    if ( l /= 0 ) then
      do i=1,l
        a=exclude_group_id(i)
        do j=1,n
          coord_temp(j) = coord(a,j) - center_coord(j)
        end do

        call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
        call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
        call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

        call get_vec_len(proj_a_vec, proj_a_vec_len, n)
        call get_vec_len(proj_b_vec, proj_b_vec_len, n)
        call get_vec_len(proj_c_vec, proj_c_vec_len, n)

        if ( DOT_PRODUCT(coord_temp,a_vec_basis) < 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
          shift_a_len = a_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,a_vec_basis) > 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
          shift_a_len = -a_vec_len
        end if

        if ( DOT_PRODUCT(coord_temp,b_vec_basis) < 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
          shift_b_len = b_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,b_vec_basis) > 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
          shift_b_len = -b_vec_len
        end if

        if ( DOT_PRODUCT(coord_temp,c_vec_basis) < 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
          shift_c_len = c_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,c_vec_basis) > 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
          shift_c_len = -c_vec_len
        end if

        do j=1,n
          coord_temp(j) = coord(a,j)
        end do

        coord_temp = coord_temp + shift_a_len*a_vec_basis + &
                     shift_b_len*b_vec_basis + shift_c_len*c_vec_basis

        shift_a_len = 0.0
        shift_b_len = 0.0
        shift_c_len = 0.0
        do j=1,n
          new_coord(a,j) = coord_temp(j)
        end do
      end do
    end if

    if (trans_type == 0) then !center the system based on the center of mass of group
      do i =1,x
        do j=1,y
          a=group_atom_1_id(i,j)
          call center_of_mass(group_atoms_mass(i,:), coord(a:(a+z-1),:), group_center_coord, z, n)
          do k=1,n
            coord_temp(k) = group_center_coord(k) - center_coord(k)
          end do

          call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
          call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
          call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

          call get_vec_len(proj_a_vec, proj_a_vec_len, n)
          call get_vec_len(proj_b_vec, proj_b_vec_len, n)
          call get_vec_len(proj_c_vec, proj_c_vec_len, n)

          if ( DOT_PRODUCT(coord_temp,a_vec_basis) < 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
            shift_a_len = a_vec_len
          end if
          if ( DOT_PRODUCT(coord_temp,a_vec_basis) > 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
            shift_a_len = -a_vec_len
          end if

          if ( DOT_PRODUCT(coord_temp,b_vec_basis) < 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
            shift_b_len = b_vec_len
          end if
          if ( DOT_PRODUCT(coord_temp,b_vec_basis) > 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
            shift_b_len = -b_vec_len
          end if

          if ( DOT_PRODUCT(coord_temp,c_vec_basis) < 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
            shift_c_len = c_vec_len
          end if
          if ( DOT_PRODUCT(coord_temp,c_vec_basis) > 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
            shift_c_len = -c_vec_len
          end if

          do k=0,z-1
            do q=1,n
              coord_temp(q) = coord(a+k,q)
            end do

            !call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
            !call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
            !call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

            !call get_vec_len(proj_a_vec, proj_a_vec_len, n)
            !call get_vec_len(proj_b_vec, proj_b_vec_len, n)
            !call get_vec_len(proj_c_vec, proj_c_vec_len, n)

            coord_temp = coord_temp + shift_a_len*a_vec_basis + &
                         shift_b_len*b_vec_basis + shift_c_len*c_vec_basis

            do q=1,n
              new_coord(a+k,q) = coord_temp(q)
            end do
          end do
          shift_a_len = 0.0
          shift_b_len = 0.0
          shift_c_len = 0.0

        end do
      end do

    else if (trans_type == 1) then
      do i=1,m
        do j=1,n
          coord_temp(j) = coord(i,j) - center_coord(j)
        end do

        call project_vec(coord_temp, a_vec_basis, n, proj_a_vec)
        call project_vec(coord_temp, b_vec_basis, n, proj_b_vec)
        call project_vec(coord_temp, c_vec_basis, n, proj_c_vec)

        call get_vec_len(proj_a_vec, proj_a_vec_len, n)
        call get_vec_len(proj_b_vec, proj_b_vec_len, n)
        call get_vec_len(proj_c_vec, proj_c_vec_len, n)

        if ( DOT_PRODUCT(coord_temp,a_vec_basis) < 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
          shift_a_len = a_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,a_vec_basis) > 0 .and. proj_a_vec_len > a_vec_len/2.0 ) then
          shift_a_len = -a_vec_len
        end if

        if ( DOT_PRODUCT(coord_temp,b_vec_basis) < 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
          shift_b_len = b_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,b_vec_basis) > 0 .and. proj_b_vec_len > b_vec_len/2.0 ) then
          shift_b_len = -b_vec_len
        end if

        if ( DOT_PRODUCT(coord_temp,c_vec_basis) < 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
          shift_c_len = c_vec_len
        end if
        if ( DOT_PRODUCT(coord_temp,c_vec_basis) > 0 .and. proj_c_vec_len > c_vec_len/2.0 ) then
          shift_c_len = -c_vec_len
        end if

        coord_temp = coord_temp + shift_a_len*a_vec_basis + &
                     shift_b_len*b_vec_basis + shift_c_len*c_vec_basis

        do k=1,n
          new_coord(a,k) = coord_temp(k)
        end do
        shift_a_len = 0.0
        shift_b_len = 0.0
        shift_c_len = 0.0
      end do

    end if

    return

  end subroutine

  subroutine rdf(distance, increment, vol, total_num_1, total_num_2, data_num, rdf_value_final, integral_value_final, u, v, w)

    integer::data_num
    integer::u, v, w
    integer::i, j, k, l
    integer::total_num_1, total_num_2
    integer::num_2_r
    real(kind=4)::density
    real(kind=4)::increment
    real(kind=4)::sum_value_1, sum_value_2
    real(kind=4)::pi
    real(kind=4)::vol
    real(kind=4),dimension(data_num)::rdf_value_final
    real(kind=4),dimension(data_num)::integral_value_final
    real(kind=4),dimension(u,v,w)::distance
    real(kind=4),dimension(u,v,data_num)::rdf_value
    real(kind=4),dimension(u,v,data_num)::integral_value

    !f2py intent(in)::u, v, w, data_num
    !f2py intent(in)::distance
    !f2py intent(in)::vol
    !f2py intent(in)::increment
    !f2py intent(in)::total_num_1, total_num_2
    !f2py intent(out)::rdf_value_final, integral_value_final

    i = 0
    j = 0
    k = 0
    l = 0
    sum_value_1 = 0.0
    sum_value_2 = 0.0
    pi=3.1415926

    density = total_num_2/vol

    do i=1,u
      do j=1,v
        do k=1,data_num
          num_2_r = 0
          do l=1,w
            if ( distance(i,j,l) > k*increment .and. distance(i,j,l) <= (k+1)*increment ) then
              num_2_r = num_2_r+1
            end if
          end do
          rdf_value(i,j,k) = num_2_r/(4*pi*density*((k+1)*increment)**2*increment*total_num_1)
        end do
      end do
    end do

    do i=1,u
      do j=1,v
        do k=1,data_num
          sum_value_1 = 0.0
          do l=1,k
            sum_value_1 = sum_value_1+4*pi*density*(increment*(l+1))**2*rdf_value(i,j,l)*increment
          end do
          integral_value(i,j,k) = sum_value_1
        end do
      end do
    end do

    do i=1,data_num
      sum_value_1 = 0.0
      do j=1,u
        sum_value_2 = 0.0
        do k=1,v
          sum_value_2 = sum_value_2+rdf_value(j,k,i)
        end do
        sum_value_1 = sum_value_1+sum_value_2/v
      end do
      rdf_value_final(i) = sum_value_1/u
    end do

    do i=1,data_num
      sum_value_1 = 0.0
      do j=1,u
        sum_value_2 = 0.0
        do k=1,v
          sum_value_2 = sum_value_2+integral_value(j,k,i)
        end do
        sum_value_1 = sum_value_1+sum_value_2/v
      end do
      integral_value_final(i) = sum_value_1/u
    end do

    return

  end subroutine rdf

end module geometry
