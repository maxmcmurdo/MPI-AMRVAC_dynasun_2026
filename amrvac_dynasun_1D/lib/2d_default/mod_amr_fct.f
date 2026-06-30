module mod_amr_fct
  implicit none
  private

  type facealloc
    double precision, dimension(:,:), pointer :: face
  end type facealloc

  type fake_neighbors
    integer :: igrid
    integer :: ipe
  end type fake_neighbors

  type(facealloc), dimension(:,:,:), allocatable, public :: pface

  type(fake_neighbors), dimension(:,:,:,:), allocatable,&
      public :: fine_neighbors

  integer, dimension(:,:,:,:), allocatable, public :: old_neighbor

  integer :: itag, isend, irecv
  integer :: nrecv, nsend, ibuf_recv, ibuf_send, ibuf_send_next
  integer, dimension(2) :: isize
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus
  double precision, allocatable :: recvbuffer(:), sendbuffer(:)

  public :: store_faces
  public :: comm_faces
  public :: end_comm_faces
  public :: deallocateBfaces
  public :: old_neighbors
  public :: prolong_2nd_stg
  public :: already_fine

contains
  !> This subroutine performs a 2nd order prolongation for a staggered field F,
  !> preserving the divergence of the coarse cell.
  !> This is useful for preserving DivF=0.
  !> If DivF=f(x), a different algorithm must be used.
  subroutine prolong_2nd_stg(sCo,sFi,ixComin1in,ixComin2in,ixComax1in,&
     ixComax2in,ixFimin1in,ixFimin2in,ixFimax1in,ixFimax2in,dxCo1,dxCo2,&
     xComin1,xComin2,dxFi1,dxFi2,xFimin1,xFimin2,ghost,fine_min1in,fine_min2in,&
     fine_max1in,fine_max2in)
    use mod_global_parameters
    use mod_physics

    logical, intent(in)          :: ghost
    integer, intent(in)          :: ixComin1in,ixComin2in,ixComax1in,&
       ixComax2in, ixFimin1in,ixFimin2in,ixFimax1in,ixFimax2in
    double precision, intent(in) :: dxCo1,dxCo2, xComin1,xComin2, dxFi1,dxFi2,&
        xFimin1,xFimin2
    type(state), intent(in)      :: sCo
    type(state), intent(inout)   :: sFi

    logical, optional :: fine_min1in,fine_min2in,fine_max1in,fine_max2in
    logical           :: fine_min1,fine_min2,fine_max1,fine_max2

    double precision :: eta1,eta2, invdxCo1,invdxCo2
    integer :: ixComin1,ixComin2,ixComax1,ixComax2,ixFimin1,ixFimin2,ixFimax1,&
       ixFimax2
    integer :: idim1,idim2,ix2,idim3,ixFismin1,ixFismin2,ixFismax1,ixFismax2,&
       ixGsmin1,ixGsmin2,ixGsmax1,ixGsmax2,ixCosmin1,ixCosmin2,ixCosmax1,&
       ixCosmax2,ixFisCmin1,ixFisCmin2,ixFisCmax1,ixFisCmax2
    integer :: ixCosVmin1(1:ndim),ixCosVmin2(1:ndim),ixCosVmax1(1:ndim),&
       ixCosVmax2(1:ndim),ixFisVmin1(1:ndim),ixFisVmin2(1:ndim),&
       ixFisVmax1(1:ndim),ixFisVmax2(1:ndim)
    integer :: hxCosmin1,hxCosmin2,hxCosmax1,hxCosmax2,jxCosmin1,jxCosmin2,&
       jxCosmax1,jxCosmax2,ixCosEmin1,ixCosEmin2,ixCosEmax1,ixCosEmax2,&
       ixFisEmin1,ixFisEmin2,ixFisEmax1,ixFisEmax2,hxFisCmin1,hxFisCmin2,&
       hxFisCmax1,hxFisCmax2,jxFisCmin1,jxFisCmin2,jxFisCmax1,jxFisCmax2,&
       ipxFisCmin1,ipxFisCmin2,ipxFisCmax1,ipxFisCmax2,ixCosCmin1,ixCosCmin2,&
       ixCosCmax1,ixCosCmax2,imxFisCmin1,imxFisCmin2,imxFisCmax1,imxFisCmax2,&
       jpxFisCmin1,jpxFisCmin2,jpxFisCmax1,jpxFisCmax2,jmxFisCmin1,jmxFisCmin2,&
       jmxFisCmax1,jmxFisCmax2,hpxFisCmin1,hpxFisCmin2,hpxFisCmax1,hpxFisCmax2
    integer :: hxFimin1,hxFimin2,hxFimax1,hxFimax2,jxFimin1,jxFimin2,jxFimax1,&
       jxFimax2,hijxFimin1,hijxFimin2,hijxFimax1,hijxFimax2,hjixFimin1,&
       hjixFimin2,hjixFimax1,hjixFimax2,hjjxFimin1,hjjxFimin2,hjjxFimax1,&
       hjjxFimax2
    integer :: iihxFimin1,iihxFimin2,iihxFimax1,iihxFimax2,iijxFimin1,&
       iijxFimin2,iijxFimax1,iijxFimax2,ijhxFimin1,ijhxFimin2,ijhxFimax1,&
       ijhxFimax2,ijjxFimin1,ijjxFimin2,ijjxFimax1,ijjxFimax2,ihixFimin1,&
       ihixFimin2,ihixFimax1,ihixFimax2,ijixFimin1,ijixFimin2,ijixFimax1,&
       ijixFimax2,ihjxFimin1,ihjxFimin2,ihjxFimax1,ihjxFimax2
    integer :: jihxFimin1,jihxFimin2,jihxFimax1,jihxFimax2,jijxFimin1,&
       jijxFimin2,jijxFimax1,jijxFimax2,jjhxFimin1,jjhxFimin2,jjhxFimax1,&
       jjhxFimax2,jjjxFimin1,jjjxFimin2,jjjxFimax1,jjjxFimax2,jhixFimin1,&
       jhixFimin2,jhixFimax1,jhixFimax2,jjixFimin1,jjixFimin2,jjixFimax1,&
       jjixFimax2,jhjxFimin1,jhjxFimin2,jhjxFimax1,jhjxFimax2
    double precision :: bfluxCo(sCo%ixGsmin1:sCo%ixGsmax1,&
       sCo%ixGsmin2:sCo%ixGsmax2,nws),bfluxFi(sFi%ixGsmin1:sFi%ixGsmax1,&
       sFi%ixGsmin2:sFi%ixGsmax2,nws)
    double precision :: slopes(sCo%ixGsmin1:sCo%ixGsmax1,&
       sCo%ixGsmin2:sCo%ixGsmax2,ndim),B_energy_change(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2)
    

    

    
    ! Note on the indices:
    ! ixCo  Cells where 
    !       divergence-preserving prolongation will be applied.
    ! ixFi  Fine cells which correspond to that extent.
    ! ixCoE For 'expanded', faces that need to be used to calculate
    !       the slopes for the prolongation
    ! ixCos
    ! ixFis
    ! ixCosV And
 !ixFisV For 'variable', since the ranges depend on the direction of the faces
    ! ixCosC For 'copy',
 !ixFisC For 'copy', to fill the information from the fine grid needed in the internal faces prolongation step.
    ! At the end, ixFisC is used to copy the assign the magnetic
    ! fields components to the state structure.

    fine_min1=.false.;fine_min2=.false.;fine_max1=.false.;fine_max2=.false.;

    if(present(fine_min1in)) fine_min1=fine_min1in
    if(present(fine_min2in)) fine_min2=fine_min2in;
    if(present(fine_max1in)) fine_max1=fine_max1in
    if(present(fine_max2in)) fine_max2=fine_max2in;

    ! When filling ghost cells, ixFimin1,ixFimin2,ixFimax1,ixFimax2 are given.
    ! When refining the block, ixComin1,ixComin2,ixComax1,ixComax2 are given.
    if (ghost) then
      ixFimin1=ixFimin1in;ixFimin2=ixFimin2in;ixFimax1=ixFimax1in
      ixFimax2=ixFimax2in;
      ixComin1=int((ixFimin1+nghostcells+1)/2)
      ixComin2=int((ixFimin2+nghostcells+1)/2);
      ixComax1=int((ixFimax1+nghostcells+1)/2)
      ixComax2=int((ixFimax2+nghostcells+1)/2);
    else
      ixComin1=ixComin1in;ixComin2=ixComin2in;ixComax1=ixComax1in
      ixComax2=ixComax2in;
      ixFimin1=ixMlo1;ixFimin2=ixMlo2;ixFimax1=ixMhi1;ixFimax2=ixMhi2;
    end if

    ! Expanded range for staggered variables
    ixCosmin1=ixComin1-1;ixCosmin2=ixComin2-1;
    ixCosmax1=ixComax1;ixCosmax2=ixComax2;

    ixFismin1=ixFimin1-1;ixFismin2=ixFimin2-1;
    ixFismax1=ixFimax1;ixFismax2=ixFimax2;


    associate(wCos=>sCo%ws, wFis=>sFi%ws,wCo=>sCo%w, wFi=>sFi%w)
    ! Assemble general indices
    ixGsmin1=sFi%ixGsmin1;ixGsmin2=sFi%ixGsmin2;
    ixGsmax1=sFi%ixGsmax1;ixGsmax2=sFi%ixGsmax2;

    do idim1=1,ndim
      ixCosVmin1(idim1)=ixComin1-kr(1,idim1)
      ixCosVmin2(idim1)=ixComin2-kr(2,idim1);
      ixCosVmax1(idim1)=ixComax1;ixCosVmax2(idim1)=ixComax2;
      ixFisVmin1(idim1)=ixFimin1-kr(1,idim1)
      ixFisVmin2(idim1)=ixFimin2-kr(2,idim1);
      ixFisVmax1(idim1)=ixFimax1;ixFisVmax2(idim1)=ixFimax2;
    end do

    ! Initialize auxiliary arrays at zero
    bfluxCo = zero
    bfluxFi = zero 
    slopes  = zero 


    invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;

    ! Fill coarse magnetic flux array
    do idim1=1,ndim
      ! Fill information from parts already at the fine level
      ! First set up indices
      if (ghost) then
        ixFisEmin1=max(1-kr(1,idim1),ixFisVmin1(idim1)-2*(1-kr(1,idim1)))
        ixFisEmin2=max(1-kr(2,idim1),ixFisVmin2(idim1)-2*(1-kr(2,idim1)));
        ixFisEmax1=min(ixGsmax1,ixFisVmax1(idim1)+2*(1-kr(1,idim1)))
        ixFisEmax2=min(ixGsmax2,ixFisVmax2(idim1)+2*(1-kr(2,idim1)));
        ixCosEmin1=int((ixFisEmin1+nghostcells+1)/2)
        ixCosEmin2=int((ixFisEmin2+nghostcells+1)/2)
        ixCosEmax1=int((ixFisEmax1+nghostcells+1)/2)
        ixCosEmax2=int((ixFisEmax2+nghostcells+1)/2);
      else
        ixCosEmin1=ixCosVmin1(idim1)-(1-kr(idim1,1))
        ixCosEmin2=ixCosVmin2(idim1)-(1-kr(idim1,2))
        ixCosEmax1=ixCosVmax1(idim1)+(1-kr(idim1,1))
        ixCosEmax2=ixCosVmax2(idim1)+(1-kr(idim1,2));
        ixFisEmin1=ixFisVmin1(idim1)-2*(1-kr(idim1,1))
        ixFisEmin2=ixFisVmin2(idim1)-2*(1-kr(idim1,2))
        ixFisEmax1=ixFisVmax1(idim1)+2*(1-kr(idim1,1))
        ixFisEmax2=ixFisVmax2(idim1)+2*(1-kr(idim1,2));
      end if
      ! Convert fine fields to fluxes
      bfluxFi(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
         idim1)=wFis(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
         idim1)*sFi%surfaceC(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
         idim1)
      
      idim2=1+mod(idim1,2)
     
      
      bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,idim1) = zero
      ! Add fine fluxes sharing the same fine face
     do ix2=0,1
        ixFisCmin1=ixFisEmin1+ix2*kr(idim2,1)
        ixFisCmin2=ixFisEmin2+ix2*kr(idim2,2)
        ixFisCmax1=ixFisEmax1+ix2*kr(idim2,1)
        ixFisCmax2=ixFisEmax2+ix2*kr(idim2,2);
        
        bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
           idim1)=bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
           idim1)+bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
           idim1)
     end do
    end do

    ! Omit indices for already refined face, if any
    
    if (fine_min1) then
      ixCosVmin1(1)=ixCosVmin1(1)+1
      ixFisVmin1(1)=ixFisVmin1(1)+2
    end if
    if (fine_max1) then
      ixCosVmax1(1)=ixCosVmax1(1)-1
      ixFisVmax1(1)=ixFisVmax1(1)-2
    end if
    
    
    if (fine_min2) then
      ixCosVmin2(2)=ixCosVmin2(2)+1
      ixFisVmin2(2)=ixFisVmin2(2)+2
    end if
    if (fine_max2) then
      ixCosVmax2(2)=ixCosVmax2(2)-1
      ixFisVmax2(2)=ixFisVmax2(2)-2
    end if
    

    do idim1=1,ndim
      ixCosEmin1=ixCosVmin1(idim1);ixCosEmin2=ixCosVmin2(idim1)
      ixCosEmax1=ixCosVmax1(idim1);ixCosEmax2=ixCosVmax2(idim1);
      ! Omit part already refined
      
      if (1/=idim1) then
        if ((.not.fine_min1).or.(.not.ghost)) then
         ixCosEmin1=ixCosVmin1(idim1)-1
        end if
        if ((.not.fine_max1).or.(.not.ghost)) then
         ixCosEmax1=ixCosVmax1(idim1)+1
        end if
      end if
      
    
      if (2/=idim1) then
        if ((.not.fine_min2).or.(.not.ghost)) then
         ixCosEmin2=ixCosVmin2(idim1)-1
        end if
        if ((.not.fine_max2).or.(.not.ghost)) then
         ixCosEmax2=ixCosVmax2(idim1)+1
        end if
      end if
      
      ! Fill coarse flux array from coarse field
      bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         idim1)=wCos(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         idim1)*sCo%surfaceC(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         idim1)
    end do
    ! Finished filling coarse flux array

    ! Distribute coarse fluxes among fine fluxes
    ! There are too many loops here, perhaps optimize later
    do idim1=1,ndim
       ixCosmin1=ixCosVmin1(idim1);ixCosmin2=ixCosVmin2(idim1)
       ixCosmax1=ixCosVmax1(idim1);ixCosmax2=ixCosVmax2(idim1);
       do idim2=1,ndim
         if(idim1==idim2) cycle
         ! Calculate slope in direction idim2
         ! Set up indices
         jxCosmin1=ixCosmin1+kr(idim2,1);jxCosmin2=ixCosmin2+kr(idim2,2)
         jxCosmax1=ixCosmax1+kr(idim2,1);jxCosmax2=ixCosmax2+kr(idim2,2);
         hxCosmin1=ixCosmin1-kr(idim2,1);hxCosmin2=ixCosmin2-kr(idim2,2)
         hxCosmax1=ixCosmax1-kr(idim2,1);hxCosmax2=ixCosmax2-kr(idim2,2);
         slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
            idim2)=0.125d0*(bfluxCo(jxCosmin1:jxCosmax1,jxCosmin2:jxCosmax2,&
            idim1)-bfluxCo(hxCosmin1:hxCosmax1,hxCosmin2:hxCosmax2,idim1))
    
         do ix2=0,1
           ixFisCmin1=ixFisVmin1(idim1)+ix2*kr(1,idim2)
           ixFisCmin2=ixFisVmin2(idim1)+ix2*kr(2,idim2);
           ixFisCmax1=ixFisVmax1(idim1)+ix2*kr(1,idim2)
           ixFisCmax2=ixFisVmax2(idim1)+ix2*kr(2,idim2);
           bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
              idim1)=half*(bfluxCo(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
              idim1)+(2*ix2-1)*slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
              idim2))
         end do
   
    
       end do
    end do

    ! Calculate interior fine fluxes
    

    do idim1=1,ndim
      do idim2=1,ndim
        
        if (idim1==idim2) cycle
        ixFisCmin1=ixFismin1+1;ixFisCmin2=ixFismin2+1;
        ixFisCmax1=ixFismax1-1;ixFisCmax2=ixFismax2-1;
        jxFisCmin1=ixFisCmin1+kr(idim1,1);jxFisCmin2=ixFisCmin2+kr(idim1,2)
        jxFisCmax1=ixFisCmax1+kr(idim1,1);jxFisCmax2=ixFisCmax2+kr(idim1,2);
        hxFisCmin1=ixFisCmin1-kr(idim1,1);hxFisCmin2=ixFisCmin2-kr(idim1,2)
        hxFisCmax1=ixFisCmax1-kr(idim1,1);hxFisCmax2=ixFisCmax2-kr(idim1,2);
        ipxFisCmin1=ixFisCmin1+kr(idim2,1);ipxFisCmin2=ixFisCmin2+kr(idim2,2)
        ipxFisCmax1=ixFisCmax1+kr(idim2,1);ipxFisCmax2=ixFisCmax2+kr(idim2,2);
        imxFisCmin1=ixFisCmin1-kr(idim2,1);imxFisCmin2=ixFisCmin2-kr(idim2,2)
        imxFisCmax1=ixFisCmax1-kr(idim2,1);imxFisCmax2=ixFisCmax2-kr(idim2,2);
        jpxFisCmin1=jxFisCmin1+kr(idim2,1);jpxFisCmin2=jxFisCmin2+kr(idim2,2)
        jpxFisCmax1=jxFisCmax1+kr(idim2,1);jpxFisCmax2=jxFisCmax2+kr(idim2,2);
        jmxFisCmin1=jxFisCmin1-kr(idim2,1);jmxFisCmin2=jxFisCmin2-kr(idim2,2)
        jmxFisCmax1=jxFisCmax1-kr(idim2,1);jmxFisCmax2=jxFisCmax2-kr(idim2,2);
        hpxFisCmin1=hxFisCmin1+kr(idim2,1);hpxFisCmin2=hxFisCmin2+kr(idim2,2)
        hpxFisCmax1=hxFisCmax1+kr(idim2,1);hpxFisCmax2=hxFisCmax2+kr(idim2,2);

        bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
           idim1)=half*(bfluxFi(jxFisCmin1:jxFisCmax1:2,&
           jxFisCmin2:jxFisCmax2:2,idim1)+bfluxFi(hxFisCmin1:hxFisCmax1:2,&
           hxFisCmin2:hxFisCmax2:2,idim1))-&
           quarter*(bfluxFi(ipxFisCmin1:ipxFisCmax1:2,&
           ipxFisCmin2:ipxFisCmax2:2,idim2)-bfluxFi(jpxFisCmin1:jpxFisCmax1:2,&
           jpxFisCmin2:jpxFisCmax2:2,idim2)-bfluxFi(imxFisCmin1:imxFisCmax1:2,&
           imxFisCmin2:imxFisCmax2:2,idim2)+bfluxFi(jmxFisCmin1:jmxFisCmax1:2,&
           jmxFisCmin2:jmxFisCmax2:2,idim2))

        bfluxFi(ipxFisCmin1:ipxFisCmax1:2,ipxFisCmin2:ipxFisCmax2:2,&
           idim1)=half*(bfluxFi(jpxFisCmin1:jpxFisCmax1:2,&
           jpxFisCmin2:jpxFisCmax2:2,idim1)+bfluxFi(hpxFisCmin1:hpxFisCmax1:2,&
           hpxFisCmin2:hpxFisCmax2:2,idim1))-&
           quarter*(bfluxFi(ipxFisCmin1:ipxFisCmax1:2,&
           ipxFisCmin2:ipxFisCmax2:2,idim2)-bfluxFi(jpxFisCmin1:jpxFisCmax1:2,&
           jpxFisCmin2:jpxFisCmax2:2,idim2)-bfluxFi(imxFisCmin1:imxFisCmax1:2,&
           imxFisCmin2:imxFisCmax2:2,idim2)+bfluxFi(jmxFisCmin1:jmxFisCmax1:2,&
           jmxFisCmin2:jmxFisCmax2:2,idim2))
       
        
      end do
    end do

    ! Go back to magnetic fields
    do idim1=1,ndim
      ixFisCmax1=ixFimax1;ixFisCmax2=ixFimax2;
      ixFisCmin1=ixFimin1-kr(1,idim1);ixFisCmin2=ixFimin2-kr(2,idim1);
      where(sFi%surfaceC(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
         idim1)/=zero)
        wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
           idim1)=bfluxFi(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
           idim1)/sFi%surfaceC(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
           idim1)
      elsewhere
        wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,idim1)=zero
      end where
    end do

    if(phys_total_energy.and. .not.prolongprimitive) then
      B_energy_change(ixFimin1:ixFimax1,ixFimin2:ixFimax2)=0.5d0*sum(wFi(&
         ixFimin1:ixFimax1,ixFimin2:ixFimax2,iw_mag(:))**2,dim=ndim+1)
    end if
    call phys_face_to_center(ixFimin1,ixFimin2,ixFimax1,ixFimax2,sFi)
    if(phys_total_energy.and. .not.prolongprimitive) then
      B_energy_change(ixFimin1:ixFimax1,ixFimin2:ixFimax2)=0.5d0*sum(wFi(&
         ixFimin1:ixFimax1,ixFimin2:ixFimax2,iw_mag(:))**2,&
         dim=ndim+1)-B_energy_change(ixFimin1:ixFimax1,ixFimin2:ixFimax2)
      wFi(ixFimin1:ixFimax1,ixFimin2:ixFimax2,iw_e)=wFi(ixFimin1:ixFimax1,&
         ixFimin2:ixFimax2,iw_e)+B_energy_change(ixFimin1:ixFimax1,&
         ixFimin2:ixFimax2)
    end if

    end associate

    ! END NOONED
  end subroutine prolong_2nd_stg

  !> To achive consistency and thus conservation of divergence,
  !> when refining a block we take into account the faces of the
  !> already fine neighbours, if any. This routine stores them.
  subroutine store_faces
    use mod_forest, only: refine 
    use mod_global_parameters
    integer :: igrid, iigrid, idims, iside, ineighbor, ipe_neighbor
    integer :: nx1,nx2, i1,i2, ic1,ic2, inc1,inc2

    if (npe>1) then
      nsend_fc=0
      nrecv_fc=0
    end if

    ! Size of the block face
    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ! Check whether it is necessary to store any block face, i.e.
      ! if any coarser neighbour is going to be refined 
     do iside=1,2
        i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
        if (neighbor_pole(i1,i2,igrid)/=0) cycle
        if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
          ineighbor   =neighbor(1,i1,i2,igrid)
          ipe_neighbor=neighbor(2,i1,i2,igrid)
          if (refine(ineighbor,ipe_neighbor)) then
            allocate(pface(iside,1,igrid)%face(1,1:nx2))
            !! Store the faces
            if (iside==1) then !! left side
              pface(iside,1,igrid)%face(1,1:nx2)=ps(igrid)%ws(ixMlo1-1,&
                 ixMlo2:ixMhi2,1)
            else !! right side
              pface(iside,1,igrid)%face(1,1:nx2)=ps(igrid)%ws(ixMhi1,&
                 ixMlo2:ixMhi2,1)
            end if
            if (ipe_neighbor/=mype) nsend_fc(1)=nsend_fc(1)+1
          end if
        end if
      end do
     do iside=1,2
        i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
        if (neighbor_pole(i1,i2,igrid)/=0) cycle
        if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
          ineighbor   =neighbor(1,i1,i2,igrid)
          ipe_neighbor=neighbor(2,i1,i2,igrid)
          if (refine(ineighbor,ipe_neighbor)) then
            allocate(pface(iside,2,igrid)%face(1:nx1,1))
            !! Store the faces
            if (iside==1) then !! left side
              pface(iside,2,igrid)%face(1:nx1,1)=ps(igrid)%ws(ixMlo1:ixMhi1,&
                 ixMlo2-1,2)
            else !! right side
              pface(iside,2,igrid)%face(1:nx1,1)=ps(igrid)%ws(ixMlo1:ixMhi1,&
                 ixMhi2,2)
            end if
            if (ipe_neighbor/=mype) nsend_fc(2)=nsend_fc(2)+1
          end if
        end if
      end do

      ! If a grid is going to be refined,
      ! remember what are its neighbours.
      if (refine(igrid,mype)) then
        ! Initialize neighbour array
        fine_neighbors(:,:,:,igrid)%igrid=-1
        fine_neighbors(:,:,:,igrid)%ipe=-1
        do idims=1,ndim
          do iside=1,2
            i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
            if (neighbor_pole(i1,i2,igrid)/=0) cycle
            if (neighbor_type(i1,i2,igrid)==neighbor_fine) then
             do ic2=1+int((1+i2)/2),2-int((1-i2)/2)
                inc2=ic2+i2
             do ic1=1+int((1+i1)/2),2-int((1-i1)/2)
                inc1=ic1+i1
                ineighbor=neighbor_child(1,inc1,inc2,igrid)
                ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)

                fine_neighbors(ic1,ic2,idims,igrid)%igrid= ineighbor
                fine_neighbors(ic1,ic2,idims,igrid)%ipe=ipe_neighbor

                if (ipe_neighbor/=mype) nrecv_fc(idims)=nrecv_fc(idims)+1
             end do
             end do
            end if
          end do
        end do
      end if

    end do

  end subroutine store_faces

  !> When refining a coarse block with fine neighbours, it is necessary
  !> prolong consistently with the already fine faces.
  !> This routine takes care of the communication of such faces.
  subroutine comm_faces
    use mod_forest, only: refine
    use mod_global_parameters

    integer                   :: iigrid,igrid,ineighbor,ipe_neighbor
    integer                   :: idims,iside,i1,i2,ic1,ic2,inc1,inc2,nx1,nx2
    integer                   :: recvsize, sendsize

    ! Communicate the block faces to achieve consistency when refining
    ! Initialize communication structures

    nrecv=0
    nsend=0
    recvsize=0
    sendsize=0

    do idims=1,ndim
       select case (idims)
       case (1)
          nrecv=nrecv+nrecv_fc(1)
          nsend=nsend+nsend_fc(1)
          nx1=1;nx2=ixMhi2-ixMlo2+1;
          isize(1)=nx1*nx2
          recvsize=recvsize+nrecv_fc(1)*isize(1)
          sendsize=sendsize+nsend_fc(1)*isize(1)
       
       case (2)
          nrecv=nrecv+nrecv_fc(2)
          nsend=nsend+nsend_fc(2)
          nx1=ixMhi1-ixMlo1+1;nx2=1;
          isize(2)=nx1*nx2
          recvsize=recvsize+nrecv_fc(2)*isize(2)
          sendsize=sendsize+nsend_fc(2)*isize(2)
       
       end select
    end do

    if (nrecv>0) then
    ! Allocate receive buffer
      allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv),&
          recvrequest(nrecv))
      recvrequest=MPI_REQUEST_NULL
      ibuf_recv=1
      irecv=0

    ! Receive
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if (refine(igrid,mype)) then
         do ic2=1,2
         do ic1=1,2
            ! Only one of the sides will be necessary,
            ! so we do the loop only over dimensions, instead of
            ! over dimensions and sizes as in the routines
            ! old_neighbors and already_fine.
            do idims=1,ndim
              ipe_neighbor=fine_neighbors(ic1,ic2,idims,igrid)%ipe
              ineighbor   =fine_neighbors(ic1,ic2,idims,igrid)%igrid
              if (ineighbor>0.and.ipe_neighbor/=mype) then
               if (idims==1) iside=ic1
               if (idims==2) iside=ic2
                !!! Check indices
                i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
                if (neighbor_pole(i1,i2,igrid)/=0) cycle
                inc1=ic1+i1;inc2=ic2+i2;
                irecv=irecv+1
                itag=4**2*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)

                call MPI_IRECV(recvbuffer(ibuf_recv),isize(idims),&
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                   recvrequest(irecv),ierrmpi)
                ibuf_recv=ibuf_recv+isize(idims)
              end if
            end do
         end do
         end do
        end if
      end do

    end if

    if (nsend>0) then
    ! Allocate send buffer
      allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),&
         sendrequest(nsend))
      sendrequest=MPI_REQUEST_NULL
      isend=0
      ibuf_send=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ! Check whether it is necessary to store any block face, i.e.
        ! if any coarser neighbour is going to be refined 
       do iside=1,2
          i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
 !When there is a pole, faces are always zero and this is not necessary
          if (neighbor_pole(i1,i2,igrid)/=0) cycle
          if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
            ineighbor   =neighbor(1,i1,i2,igrid)
            ipe_neighbor=neighbor(2,i1,i2,igrid)
            if (refine(ineighbor,ipe_neighbor)) then
              if (ipe_neighbor/=mype) then
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2);
                inc1=-2*i1+ic1;inc2=ic2;
                itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                isend=isend+1
                ibuf_send_next=ibuf_send+isize(1)
                sendbuffer(ibuf_send:ibuf_send_next-1)=reshape(pface(iside,1,&
                   igrid)%face,(/isize(1)/))
                call MPI_ISEND(sendbuffer(ibuf_send),isize(1),&
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                   sendrequest(isend),ierrmpi)
                ibuf_send=ibuf_send_next
              end if
            end if
          end if
        end do
       do iside=1,2
          i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
 !When there is a pole, faces are always zero and this is not necessary
          if (neighbor_pole(i1,i2,igrid)/=0) cycle
          if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
            ineighbor   =neighbor(1,i1,i2,igrid)
            ipe_neighbor=neighbor(2,i1,i2,igrid)
            if (refine(ineighbor,ipe_neighbor)) then
              if (ipe_neighbor/=mype) then
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2);
                inc1=ic1;inc2=-2*i2+ic2;
                itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                isend=isend+1
                ibuf_send_next=ibuf_send+isize(2)
                sendbuffer(ibuf_send:ibuf_send_next-1)=reshape(pface(iside,2,&
                   igrid)%face,(/isize(2)/))
                call MPI_ISEND(sendbuffer(ibuf_send),isize(2),&
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                   sendrequest(isend),ierrmpi)
                ibuf_send=ibuf_send_next
              end if
            end if
          end if
        end do
      end do
    end if

    ! Waitalls
    if (nrecv>0) then
       call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
       deallocate(recvstatus,recvrequest)
       ibuf_recv=1
    end if

    if (nsend>0) then
       call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
       deallocate(sendbuffer,sendstatus,sendrequest)
    end if

  end subroutine comm_faces

  subroutine end_comm_faces
    use mod_global_parameters
    ! Deallocate receive buffer
    if (nrecv>0) deallocate(recvbuffer)
  end subroutine end_comm_faces

  subroutine deallocateBfaces
    use mod_global_parameters
    integer :: igrid, iigrid, iside

    do iigrid=1,igridstail; igrid=igrids(iigrid);
     do iside=1,2
        if (associated(pface(iside,1,igrid)%face)) then
          deallocate(pface(iside,1,igrid)%face)
        end if
      end do
     do iside=1,2
        if (associated(pface(iside,2,igrid)%face)) then
          deallocate(pface(iside,2,igrid)%face)
        end if
      end do
    end do

  end subroutine deallocateBfaces

  subroutine old_neighbors(child_igrid,child_ipe,igrid,ipe)
    use mod_global_parameters
    integer, dimension(2,2), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
    integer :: iside, i1,i2, ic1,ic2

    do ic2=1,2
    do ic1=1,2
      old_neighbor(:,:,:,child_igrid(ic1,ic2))=-1
     do iside=1,2
        if (ic1==iside) then
          i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
          old_neighbor(1,i1,i2,child_igrid(ic1,ic2))=fine_neighbors(ic1,ic2,1,&
             igrid)%igrid
          old_neighbor(2,i1,i2,child_igrid(ic1,ic2))=fine_neighbors(ic1,ic2,1,&
             igrid)%ipe
        end if
      end do
     do iside=1,2
        if (ic2==iside) then
          i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
          old_neighbor(1,i1,i2,child_igrid(ic1,ic2))=fine_neighbors(ic1,ic2,2,&
             igrid)%igrid
          old_neighbor(2,i1,i2,child_igrid(ic1,ic2))=fine_neighbors(ic1,ic2,2,&
             igrid)%ipe
        end if
      end do
    end do
    end do

  end subroutine old_neighbors

  !> This routine fills the fine faces before prolonging.
  !> It is the face equivalent of fix_conserve 
  subroutine already_fine(sFi,ichild,fine_min1,fine_min2,fine_max1,fine_max2)
    use mod_forest
    use mod_global_parameters
    type(tree_node_ptr) :: tree
    type(state) :: sFi
    integer, intent(in) :: ichild
    logical :: fine_min1,fine_min2,fine_max1,fine_max2

    integer :: ineighbor,ipe_neighbor,ibufnext
    integer :: iside,iotherside,i1,i2,nx1,nx2

    ! Size of the block face
    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;

    ! Initialise everything to zero and false
    fine_min1=.false.;fine_min2=.false.;
    fine_max1=.false.;fine_max2=.false.;
    sFi%ws=zero

    
    ! This face communication is not needed in 1D
   do iside=1,2
      i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
      ! This is not necessary at the pole.
      ! We are inside a loop over the children, so the grid index is ichild
      if (neighbor_pole(i1,i2,ichild)/=0) cycle

      ! Get old ipe and igrid of neighbour from array fake_neighbor
      ! Then plug it into the structure pfaces and get the faces 
      ineighbor   =old_neighbor(1,i1,i2,ichild)
      ipe_neighbor=old_neighbor(2,i1,i2,ichild)

      iotherside=3-iside
      if (ineighbor>0) then
        if (iside==1) then ! left side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMlo1-1,ixMlo2:ixMhi2,1)=pface(iotherside,1,&
               ineighbor)%face(1,1:nx2)
          else
            ibufnext=ibuf_recv+isize(1)
            sFi%ws(ixMlo1-1,ixMlo2:ixMhi2,&
               1)=reshape(source=recvbuffer(ibuf_recv:ibufnext-1),&
               shape=shape(sFi%ws(ixMlo1-1,ixMlo2:ixMhi2,1)))
            ibuf_recv=ibufnext
          end if
          fine_min1=.true.
        else ! right side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMhi1,ixMlo2:ixMhi2,1)=pface(iotherside,1,&
               ineighbor)%face(1,1:nx2)
          else
            ibufnext=ibuf_recv+isize(1)
            sFi%ws(ixMhi1,ixMlo2:ixMhi2,1)=reshape(source=recvbuffer(&
               ibuf_recv:ibufnext-1),shape=shape(sFi%ws(ixMhi1,ixMlo2:ixMhi2,&
               1)))
            ibuf_recv=ibufnext
          end if
          fine_max1=.true.
        end if
      end if
    end do
    do iside=1,2
      i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
      ! This is not necessary at the pole.
      ! We are inside a loop over the children, so the grid index is ichild
      if (neighbor_pole(i1,i2,ichild)/=0) cycle

      ! Get old ipe and igrid of neighbour from array fake_neighbor
      ! Then plug it into the structure pfaces and get the faces 
      ineighbor   =old_neighbor(1,i1,i2,ichild)
      ipe_neighbor=old_neighbor(2,i1,i2,ichild)

      iotherside=3-iside
      if (ineighbor>0) then
        if (iside==1) then ! left side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMlo1:ixMhi1,ixMlo2-1,2)=pface(iotherside,2,&
               ineighbor)%face(1:nx1,1)
          else
            ibufnext=ibuf_recv+isize(2)
            sFi%ws(ixMlo1:ixMhi1,ixMlo2-1,&
               2)=reshape(source=recvbuffer(ibuf_recv:ibufnext-1),&
               shape=shape(sFi%ws(ixMlo1:ixMhi1,ixMlo2-1,2)))
            ibuf_recv=ibufnext
          end if
          fine_min2=.true.
        else ! right side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMlo1:ixMhi1,ixMhi2,2)=pface(iotherside,2,&
               ineighbor)%face(1:nx1,1)
          else
            ibufnext=ibuf_recv+isize(2)
            sFi%ws(ixMlo1:ixMhi1,ixMhi2,2)=reshape(source=recvbuffer(&
               ibuf_recv:ibufnext-1),shape=shape(sFi%ws(ixMlo1:ixMhi1,ixMhi2,&
               2)))
            ibuf_recv=ibufnext
          end if
          fine_max2=.true.
        end if
      end if
    end do
   
  end subroutine already_fine

end module mod_amr_fct
