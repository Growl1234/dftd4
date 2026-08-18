! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_ncoord
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   implicit none
   private

   public :: get_coordination_number, add_coordination_number_derivs
   public :: add_coordination_number_hessian


   !> Steepness of counting function
   real(wp), parameter :: default_kcn = 7.5_wp


contains


!> Geometric fractional coordination number, supports error function counting.
subroutine get_coordination_number(mol, trans, cutoff, rcov, en, cn, dcndr, dcndL)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_coordination_number

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%dftd4, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov, en=en, error=error)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

end subroutine get_coordination_number


subroutine add_coordination_number_derivs(mol, trans, cutoff, rcov, en, dEdcn, gradient, sigma)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(inout) :: gradient(:, :)

   !> Derivative of the CN with respect to strain deformations
   real(wp), intent(inout) :: sigma(:, :)


   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%dftd4, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov, en=en, error=error)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%add_coordination_number_derivs(mol, trans, dEdcn, gradient, sigma)

end subroutine add_coordination_number_derivs


!> Add the second derivative of the D4 coordination number contracted with
!> the derivative of the energy w.r.t. the coordination number.
!>
!> The counting function is reproduced here because mctc-lib currently
!> provides derivatives only up to first order.
subroutine add_coordination_number_hessian(mol, trans, cutoff, rcov, en, dEdcn, hessian)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: hessian(:, :)

   integer :: ipair, npair, iat, jat, izp, jzp, itr, ic, jc, ii, jj
   real(wp) :: vec(3), r2, r, cutoff2, rc, exponent, expterm
   real(wp) :: den, dcf, d2cf, dEdcnij, block(3, 3)
   real(wp), allocatable :: diagonal_local(:, :, :)
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: k4 = 4.10451_wp
   real(wp), parameter :: k5 = 19.08857_wp
   real(wp), parameter :: k6 = 2.0_wp*11.28174_wp**2

   cutoff2 = cutoff*cutoff
   npair = mol%nat*(mol%nat - 1)/2

   ! Every off-diagonal Cartesian block belongs to exactly one unordered atom
   ! pair.  It can therefore be written directly by the thread owning that pair;
   ! only the diagonal blocks require a thread-private reduction.  This keeps the
   ! OpenMP scratch storage O(N) instead of O(N^2) per thread.
   !$omp parallel default(none) &
   !$omp shared(mol, trans, cutoff2, rcov, en, dEdcn, hessian, npair) &
   !$omp private(ipair, iat, jat, izp, jzp, itr, ic, jc, ii, jj, vec, r2, r, rc, &
   !$omp& exponent, expterm, den, dcf, d2cf, dEdcnij, block, diagonal_local)
   allocate(diagonal_local(3, 3, mol%nat), source=0.0_wp)
   !$omp do schedule(guided, 1)
   do ipair = 1, npair
      iat = int(0.5_wp*(1.0_wp + sqrt(8.0_wp*real(ipair, wp) + 1.0_wp)))
      if (iat*(iat - 1)/2 < ipair) iat = iat + 1
      jat = ipair - (iat - 1)*(iat - 2)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      rc = rcov(izp) + rcov(jzp)
      den = k4*exp(-(abs(en(izp) - en(jzp)) + k5)**2/k6)
      dEdcnij = dEdcn(iat) + dEdcn(jat)

      do itr = 1, size(trans, 2)
         vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
         r2 = sum(vec*vec)
         if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
         r = sqrt(r2)

         exponent = default_kcn*(r - rc)/rc
         expterm = exp(-exponent*exponent)
         dcf = -den*default_kcn*expterm/(sqrtpi*rc)
         d2cf = 2.0_wp*den*default_kcn**3*(r - rc)*expterm/(sqrtpi*rc**3)

         do ic = 1, 3
            do jc = 1, 3
               block(ic, jc) = dEdcnij * ( &
                  & d2cf*vec(ic)*vec(jc)/r2 &
                  & - dcf*vec(ic)*vec(jc)/(r2*r))
            end do
            block(ic, ic) = block(ic, ic) + dEdcnij*dcf/r
         end do

         diagonal_local(:, :, iat) = diagonal_local(:, :, iat) + block
         diagonal_local(:, :, jat) = diagonal_local(:, :, jat) + block
         do ic = 1, 3
            ii = 3*(iat - 1) + ic
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               ! These two off-diagonal blocks are private to (iat,jat).
               hessian(ii, jj) = hessian(ii, jj) - block(ic, jc)
               hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) = &
                  & hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) - block(ic, jc)
            end do
         end do
      end do
   end do
   !$omp end do nowait

   ! Only 9*N values are reduced per thread; overlap the short reduction with
   ! threads that are still processing atom pairs.
   !$omp critical (add_coordination_number_hessian_)
   do iat = 1, mol%nat
      do ic = 1, 3
         ii = 3*(iat - 1) + ic
         do jc = 1, 3
            hessian(ii, 3*(iat - 1) + jc) = &
               & hessian(ii, 3*(iat - 1) + jc) + diagonal_local(ic, jc, iat)
         end do
      end do
   end do
   !$omp end critical (add_coordination_number_hessian_)
   deallocate(diagonal_local)
   !$omp end parallel

end subroutine add_coordination_number_hessian


end module dftd4_ncoord
