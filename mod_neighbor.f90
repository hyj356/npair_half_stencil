module neighbor_list
  use parameter, only: BIAS_FULL, WP, STDOUT
  implicit none

  private
  public :: model, WP, STDOUT, sort

  integer, parameter, public :: MAXNEIGHBOR = 100    !< 定义截断半径以内, 每个原子最多拥有多少个邻居原子

  type model
    real(WP), allocatable, private :: x(:)   !< 所有粒子的x坐标
    real(WP), allocatable, private :: y(:)   !< 所有粒子的y坐标
    real(WP), allocatable, private :: z(:)   !< 所有粒子的z坐标
    real(WP), allocatable, private :: atom_mass(:) !< 每种粒子的质量
    integer, allocatable, private :: neighbor_count(:)   !< 一维整数数组, 表示每个原子在截断半径内的邻居原子数量
    integer, allocatable, private :: neighbor_list(:)    !< 存储了每个原子的近邻ID的一维数组
    integer, allocatable, private :: type_list(:)        !< 每种原子的类型
    integer, private :: atom_types    !< 当前data文件含有多少种原子
    integer, private :: atom_count    !< 当前data文件含有多少个原子
    real(WP), private :: xlo, xhi     !< 盒子在x方向上的坐标上下限
    real(WP), private :: ylo, yhi     !< 盒子在y方向上的坐标上下限
    real(WP), private :: zlo, zhi     !< 盒子在z方向上的坐标上下限
    real(WP), private :: lx     !< 盒子在x方向上的长度
    real(WP), private :: ly     !< 盒子在y方向上的长度
    real(WP), private :: lz     !< 盒子在z方向上的长度
    real(WP), private :: lxi     !< 盒子在x方向上的长度的倒数
    real(WP), private :: lyi     !< 盒子在y方向上的长度的倒数
    real(WP), private :: lzi     !< 盒子在z方向上的长度的倒数
  contains
    procedure, pass(self), public :: init
    procedure, pass(self), public :: introduction
    procedure, pass(self), public :: build_neighbor_list
    procedure, pass(self), public :: find_neighbor
    procedure, pass(self), public :: write_neighbor_list
  end type

contains

  subroutine introduction(self)
    !! 对读取到的模型输出总结
    class(model), intent(in) :: self  !< 包含盒子和原子信息的派生变量

    if (.not. allocated(self%x)) then
      write(STDOUT, '(A)') 'Derived variables are not initialized!'
      stop
    end if

    write(STDOUT, '(T1, A, I0, A)') 'There are total ', self%atom_count, " atoms in current model."
    write(STDOUT, '(T1, A, I0, A)') 'There are total ', self%atom_types, " atom types in current model."
    write(STDOUT, '(T1, A)') 'Message about box:'
    write(STDOUT, '(T1, 4(A, F20.12, 1x))') 'xlo: ', self%xlo, 'xhi: ', self%xhi, &
                                            'lx: ', self%lx, 'lxi: ', self%lxi
    write(STDOUT, '(T1, 4(A, F20.12, 1x))') 'ylo: ', self%ylo, 'yhi: ', self%yhi, &
                                            'ly: ', self%ly, 'lyi: ', self%lyi
    write(STDOUT, '(T1, 4(A, F20.12, 1x))') 'zlo: ', self%zlo, 'zhi: ', self%zhi, &
                                            'lz: ', self%lz, 'lzi: ', self%lzi
    write(STDOUT, '(T1, A, F20.12)') 'The volume of the box: ', self%lx*self%ly*self%lz
    write(STDOUT, '(T1, A)') 'Coordinates of the first 10 atoms:'
    write(STDOUT, '(T8, A,T29, A, T50, A)') ' Position_x ', ' Position_y ', ' Position_z '
    write(STDOUT, '(3(F20.12, 1x))') self%x(1), self%y(1), self%z(1)
    write(STDOUT, '(3(F20.12, 1x))') self%x(2), self%y(2), self%z(2)
    write(STDOUT, '(3(F20.12, 1x))') self%x(3), self%y(3), self%z(3)
    write(STDOUT, '(3(F20.12, 1x))') self%x(4), self%y(5), self%z(5)
    write(STDOUT, '(3(F20.12, 1x))') self%x(5), self%y(5), self%z(5)
    write(STDOUT, '(3(F20.12, 1x))') self%x(6), self%y(6), self%z(6)
    write(STDOUT, '(3(F20.12, 1x))') self%x(7), self%y(7), self%z(7)
    write(STDOUT, '(3(F20.12, 1x))') self%x(8), self%y(8), self%z(8)
    write(STDOUT, '(3(F20.12, 1x))') self%x(9), self%y(9), self%z(9)
    write(STDOUT, '(3(F20.12, 1x))') self%x(10), self%y(10), self%z(10)
  end subroutine introduction

  subroutine init(self, filename)
    !! 从filename中初始化模型
    class(model), intent(inout) :: self       !< 包含盒子和原子信息的派生变量
    character(len=*), intent(in) :: filename  !< lmp文件名字
    logical :: is_exist   !< 判断文件是否存在
    integer :: file_ID    !< 文件的读写通道ID
    integer :: tmp, i

    ! 查询文件是否存
    inquire(file=filename, exist=is_exist)
    if (.not. is_exist) then
      write(STDOUT, '(A)') filename, " doesn't found in current path!"
      stop
    end if
    ! 只读方式打开文件
    open(newunit=file_ID, file=filename, action='read', status='old')
    ! 跳过前两行
    read(file_ID, *)
    read(file_ID, *)

    read(file_ID, *) self%atom_count  ! 读取原子数量
    read(file_ID, *) self%atom_types  ! 读取原子种类数量
    read(file_ID, *)

    ! 读取盒子信息
    read(file_ID, *) self%xlo, self%xhi
    read(file_ID, *) self%ylo, self%yhi
    read(file_ID, *) self%zlo, self%zhi
    read(file_ID, *)
    read(file_ID, *)
    read(file_ID, *)
    self%lx = self%xhi - self%xlo
    self%ly = self%yhi - self%ylo
    self%lz = self%zhi - self%zlo
    self%lxi = 1.0_WP / self%lx
    self%lyi = 1.0_WP / self%ly
    self%lzi = 1.0_WP / self%lz

    ! 读取到关键信息之后开始分配内存
    allocate(self%x(self%atom_count), self%y(self%atom_count),         &
             self%z(self%atom_count), self%atom_mass(self%atom_types), &
             source = 0.0_WP)
    allocate(self%neighbor_count(self%atom_count),                      &
             self%neighbor_list(MAXNEIGHBOR*self%atom_count),          &
             self%type_list(self%atom_count), source = 0)

    ! 读取原子质量
    do i = 1, self%atom_types
      read(file_ID, *) tmp, self%atom_mass(i)
    end do
    ! 跳过三行
    read(file_ID, *)
    read(file_ID, *)
    read(file_ID, *)
    ! 读取原子坐标
    do i = 1, self%atom_count
      read(file_ID, *) tmp, self%type_list(i), self%x(i), self%y(i), self%z(i)
    end do

    close(file_ID)
  end subroutine init

  pure function distance(self, i, j) result(res)
    !! 此函数用于返回第i个和第j个原子之间的距离, 默认xyz均为周期性边界条件
    type(model), intent(in) :: self !< 包含盒子和原子信息的派生变量
    integer, intent(in) :: i, j     !< 两个原子之间的ID
    real(WP) :: res                 !< 两个原子之间的距离的平方
    real(WP) :: dx, dy, dz

    dx = self%x(i) - self%x(j)
    dx = dx - nint(dx * self%lxi) * self%lx
    dy = self%y(i) - self%y(j)
    dy = dy - nint(dy * self%lyi) * self%ly
    dz = self%z(i) - self%z(j)
    dz = dz - nint(dz * self%lzi) * self%lz
    res = dx*dx + dy*dy + dz*dz

  end function distance

  pure function find_cell(self, i, ibin_size, num_cell) result(cell)
    !! 计算第i个原子属于哪一个cell
    type(model), intent(in) :: self   !< 包含盒子和原子信息的派生变量
    integer, intent(in) :: i          !< 原子序号
    real(WP), intent(in) :: ibin_size !< cell边长的倒数, 转换除法为乘法, 加快速度
    integer, intent(in), dimension(4) :: num_cell !< 前三个数字分别表示盒子在xyz上划分了多少份, 第四个数字表示模型总共有几个cell
    integer, dimension(4) :: cell !< 前三个数字分别表示原子i在xyz上属于第几份, 第四个数字表示原子属于哪个cell
    real(WP) :: dx, dy, dz        !< 确定原子i距离盒子xyz下边界的距离

    ! ceiling(0.0) = 0, 而Fortran中数组下标从1开始, 所以需要进行修正
    dx = self%x(i) - self%xlo
    cell(1) = ceiling(dx * ibin_size) ! 向上取整
    if (cell(1) <= 0)          cell(1) = cell(1) + num_cell(1)
    if (cell(1) > num_cell(1)) cell(1) = cell(1) - num_cell(1)

    dy = self%y(i) - self%ylo
    cell(2) = ceiling(dy * ibin_size) ! 向上取整
    if (cell(2) <= 0)          cell(2) = cell(2) + num_cell(2)
    if (cell(2) > num_cell(2)) cell(2) = cell(2) - num_cell(2)

    dz = self%z(i) - self%zlo
    cell(3) = ceiling(dz * ibin_size) ! 向上取整
    if (cell(3) <= 0)          cell(3) = cell(3) + num_cell(3)
    if (cell(3) > num_cell(3)) cell(3) = cell(3) - num_cell(3)

    ! 计算原子i最终属于哪个Cell
    cell(4) = cell(1) + num_cell(1) *                    &
             (cell(2) - 1 + num_cell(2) * (cell(3) - 1))
  end function find_cell

  subroutine build_neighbor_list(self, cutoff, style)
    !! 对近邻列表的构建子程序的封装, 因为纯子程序速度快但是不能进行任何I/O操作
    !! 所以用这个子程序封装, 首先检查数据输入, 并对传入数据确保没问题之后再传入
    !! 给对应子程序计算
    class(model), intent(inout) :: self     !< 包含盒子和原子信息的派生变量
    real(WP), intent(in) :: cutoff          !< 用户设置的截断半径
    character(len=*), intent(in) :: style   !< 构建邻居列表的方法
    integer :: ierror !< 确定第几个原子的邻居数量超出上限
    real(WP) :: cutoff2     !< 截断半径的平方
    real(WP) :: ibin_size   !< cell边长的倒数

    ! 检查是否对邻居列表分配了内存
    if (.not. allocated(self%neighbor_list)) then
      write(STDOUT, '(A)') 'Neighbor list is not initialized in memory.'
      stop
    end if
    cutoff2 = cutoff * cutoff
    if (style == 'verlet') then
      call verlet_neighbor_kernel(self, cutoff2, ierror)
      if (ierror /= 0) then
        write(STDOUT, '(A, I0, A, I0,A)') 'The upper limit of the maximum number of neighbors: ',&
         MAXNEIGHBOR, ' is exceeded when constructing the nearest neighbor list of the ', &
         ierror,'th atom.'
        stop
      end if
    else if (style=='cell') then
      ibin_size = 1.0_WP / (0.5_WP * cutoff)
      call cell_neighbor_kernel(self, cutoff2, ibin_size)
    else
      write(STDOUT, '(A, A, A, A)') 'Unknown neighbor list construction method: ', '"',style, '"'
      write(STDOUT, '(A)') 'Please choose one style from: ["verlet", "cell"]'
      stop
    end if


  end subroutine build_neighbor_list

  pure subroutine verlet_neighbor_kernel(self, cutoff2, ierror)
    !! 使用verlet-list法搜索所有原子的近邻原子, 理论时间复杂度为O(N^2)
    type(model), intent(inout) :: self  !< 包含盒子和原子信息的派生变量
    real(WP), intent(in) :: cutoff2     !< 用户设置的截断半径的平方
    integer, intent(inout) :: ierror    !< 用于判断邻居列表构建过程是否出错的整数变量, 正常为0, 如果不正常, 返回超过上限的原子ID
    integer :: i, j

    do i = 1, self%atom_count - 1
      do j = i + 1, self%atom_count
        if (distance(self, i, j) < cutoff2) then
          self%neighbor_count(i) = self%neighbor_count(i) + 1
          self%neighbor_list((i - 1)*MAXNEIGHBOR + self%neighbor_count(i)) = j
          self%neighbor_count(j) = self%neighbor_count(j) + 1
          self%neighbor_list((j - 1)*MAXNEIGHBOR + self%neighbor_count(j)) = i
          ! 如果发现某个原子的最大邻居原子数量超出了上限, 给出报错
          if (self%neighbor_count(i) > MAXNEIGHBOR) then
            ierror = i
            return
          else if (self%neighbor_count(j) > MAXNEIGHBOR) then
            ierror = j
            return
          end if
        end if
      end do
    end do
  end subroutine verlet_neighbor_kernel

  subroutine cell_neighbor_kernel(self, cutoff2, ibin_size)
    !! 利用cell_list法寻找在截断半径cutoff以内的近邻原子, 理论时间复杂度为O(N)
    class(model), intent(inout) :: self  !< 包含盒子和原子信息的派生变量
    real(WP), intent(in) :: ibin_size   !< bin_size的倒数
    real(WP), intent(in) :: cutoff2     !< 截断半径的平方
    integer :: num_cell(4)  !< 前三个数字分别表示盒子在xyz上划分了多少份, 第四个数字表示模型总共有几个cell
    integer :: cell(4)      !< 某个原子在xyz上属于第几份, 第四个数字表示对应原子属于哪个cell
    integer :: tmp_cell(3)  !< 临时计算的变量
    integer :: i, l, j !< 循环变量
    integer :: neighbor_cell !< 邻居cell的id
    integer :: neighbor_id   !< 邻居原子的id
    integer, allocatable :: cell_count(:)     !< 记录每个cell里面一共有多少个原子以确定循环次数
    integer, allocatable :: cell_count_sum(:) !< 维数与cell_count一样, 是cell_count的前缀和数组
    integer, allocatable :: cell_content(:)   !< 记录每个原子对应的cellid

    ! 计算整个模型在xyz上可以划分几份, 以及模型中总共有几个cell
    num_cell(1) = ceiling(self%lx * ibin_size)
    num_cell(2) = ceiling(self%ly * ibin_size)
    num_cell(3) = ceiling(self%lz * ibin_size)
    num_cell(4) = num_cell(1) * num_cell(2) * num_cell(3)
    write(STDOUT, *) num_cell
    ! 分配内存并进行初始化
    allocate(cell_content(self%atom_count), cell_count(num_cell(4)), &
             cell_count_sum(num_cell(4)), source = 0)
    ! 确定每个Cell里面有几个原子
    do i = 1, self%atom_count
      cell = find_cell(self, i, ibin_size, num_cell)
      cell_count(cell(4)) = cell_count(cell(4)) + 1
    end do
    ! 构建cell_count的前缀和数组, 注意cell_count_sum第一个元素为0,
    do i = 2, num_cell(4)
      cell_count_sum(i) = cell_count_sum(i-1) + cell_count(i-1)
    end do
    ! 清零cell_count数组
    cell_count = 0
    ! 将原子id按照之前计算好的数组规律依次填入
    do i = 1, self%atom_count
      cell = find_cell(self, i, ibin_size, num_cell)  ! 计算第i个原子的cell
      cell_count(cell(4)) = cell_count(cell(4)) + 1   ! 对应cell的位置, 原子数量加一
      cell_content(cell_count_sum(cell(4)) + cell_count(cell(4))) = i ! 填入cell_content对应位置
    end do
    ! 处理完成之后, 开始正式构建邻居列表
    !$omp parallel do default(shared) private(cell, tmp_cell, neighbor_cell, neighbor_id) num_threads(2)
    do i = 1, self%atom_count
      cell = find_cell(self, i, ibin_size, num_cell)  ! 计算第i个原子属于哪个cell
      do l = 1, size(BIAS_FULL, dim=2)
        tmp_cell = cell(1:3) + BIAS_FULL(1:3, l) ! 在中心原子基础上进行偏移
        if (cell(1) + BIAS_FULL(1, l) <= 0)           tmp_cell(1) = tmp_cell(1) + num_cell(1)
        if (cell(1) + BIAS_FULL(1, l) > num_cell(1))  tmp_cell(1) = BIAS_FULL(1, l)
        if (cell(2) + BIAS_FULL(2, l) <= 0)           tmp_cell(2) = tmp_cell(2) + num_cell(2)
        if (cell(2) + BIAS_FULL(2, l) > num_cell(2))  tmp_cell(2) = BIAS_FULL(2, l)
        if (cell(3) + BIAS_FULL(3, l) <= 0)           tmp_cell(3) = tmp_cell(3) + num_cell(3)
        if (cell(3) + BIAS_FULL(3, l) > num_cell(3))  tmp_cell(3) = BIAS_FULL(3, l)
        ! 计算修正之后的邻居cell的id
        neighbor_cell = tmp_cell(1) + num_cell(1) *                    &
                        (tmp_cell(2) - 1 + num_cell(2) * (tmp_cell(3) - 1))
        ! 循环邻居cell里面的所有原子
        do j = 1, cell_count(neighbor_cell)
          neighbor_id = cell_content(cell_count_sum(neighbor_cell) + j)
          if ( distance(self, i, neighbor_id) < cutoff2) then !
            ! 对应原子的邻居数量加一
            !$omp critical
            self%neighbor_count(i) = self%neighbor_count(i) + 1
            self%neighbor_list((i - 1)*MAXNEIGHBOR + self%neighbor_count(i)) = neighbor_id
            self%neighbor_count(neighbor_id) = self%neighbor_count(neighbor_id) + 1
            self%neighbor_list((neighbor_id - 1)*MAXNEIGHBOR + &
                            self%neighbor_count(neighbor_id)) = i
            !$omp end critical
          end if
        end do
      end do
    end do
    !$omp end parallel  do
    ! 释放内存
    deallocate(cell_content, cell_count, cell_count_sum)

  end subroutine cell_neighbor_kernel

  pure function find_neighbor(self, index) result(list)
    !! 寻找第index个原子的所有的近邻原子的ID
    class(model), intent(in) :: self  !< 包含盒子和原子信息的派生变量
    integer, intent(in) :: index      !< 第几个原子的ID
    integer, allocatable :: list(:)   !< 第index个原子的所有邻居的原子ID

    allocate(list(self%neighbor_count(index)), &
             source=self%neighbor_list(        &
             (index-1)*MAXNEIGHBOR+1 :         &
             (index-1)*MAXNEIGHBOR+self%neighbor_count(index)))
  end function

  subroutine write_neighbor_list(self, filename)
    !! 将计算的结果写入到以filename命名的txt文件中
    class(model), intent(in) :: self  !< 包含盒子和原子信息的派生变量
    character(len=*), intent(in) :: filename  !< 输出的txt文件的名字
    integer :: file_ID  !< txt文件对应的读写通道ID
    integer :: i  !< 循环变量

    open(newunit=file_ID, file=filename, action='write')
    do i = 1, self%atom_count
      write(file_ID, '(*(I0, 1x))') &
      sort(self%neighbor_list(           &
      (i-1)*MAXNEIGHBOR+1 :         &
      (i-1)*MAXNEIGHBOR+self%neighbor_count(i)))
    end do
    close(file_ID)
  end subroutine write_neighbor_list

  function sort(array) result(ret)
    integer, intent(in) :: array(:)
    integer, allocatable :: ret(:)   !< 第index个原子的所有邻居的原子ID
    integer :: i, j, n
    integer :: tmp

    n = size(array)
    allocate(ret, source=array)
    do i = 1, n - 1
        do j = i + 1, n
            if(ret(i) > ret(j))then
                tmp = ret(i)
                ret(i) = ret(j)
                ret(j) = tmp
            end if
        end do
    end do
  end function

end module neighbor_list
