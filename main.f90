program main
  ! gfortran mod_parameter.f90 mod_neighbor.f90 main.f90 -o target -O3
  use iso_fortran_env, only: LL => int64
  use parameter, only: WP, STDOUT
  use neighbor_list, only: model, sort
  implicit none
  type(model) :: Al
  character(len=*), parameter :: filename = './Al3.lmp'
  integer(LL) :: start_tick, end_tick, system_rate

  call Al%init(filename)
  call Al%introduction()
  call system_clock(count=start_tick, count_rate=system_rate)
  call Al%build_neighbor_list(cutoff=0.8536_WP*4.046_WP, style='cell') !0.8536_WP*4.046_WP
  call system_clock(count=end_tick, count_rate=system_rate)
  write(STDOUT, '(A, G0, A)') "Time cost: ", (end_tick-start_tick)/real(system_rate, 8), ' s.'
  write(STDOUT, '(A)') 'Neighbor atom IDs of the first atoms:'
  write(STDOUT, '(*(I0,  ", "))') sort(Al%find_neighbor(1))

  !call Al%write_neighbor_list('neighbor_list_fortran.txt')  ! 将邻居列表写出到txt文件中与ovito计算结果对比
end program main
