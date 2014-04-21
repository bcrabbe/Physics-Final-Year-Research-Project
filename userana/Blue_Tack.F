      SUBROUTINE Blue_Tack (nchar,tmom,Vxyz,vpnt,vdir)
c.
c.  Author: Alan Norton
c.  Corrects a track for the blue tube field.   !  revised 15/11/04 <<<<
c.  The track charge is input as nchar = charge[+1,-1], tmom = momentum,
c.  vpnt = (x,y) at Dch1 and vdir = (dxdz,dydz) at Dch1.
c.  Vxyz (x,y,z) is an estimate of the vertex.
C.  On output vpnt and vdir are replaced by effective values which can
c.  be linearly extrapolated back to the vertex region.
c.
c.------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER    nchar,imap
      REAL*4     tmom, Vxyz(3),vpnt(2),vdir(2)
     +,          charmom,Vtem(3),binteg(6)
c
      charmom = float(nchar) * tmom
      Vtem(1) = Vxyz(1)
      Vtem(2) = Vxyz(2)
      Vtem(3) = Vxyz(3) + 12000.
c
      Call FDCALINT_VT (Vtem,vdir,binteg,imap,2)
      vpnt(1) = vpnt(1) - 0.0003*binteg(1)/charmom
      vdir(1) = vdir(1) - 0.0003*binteg(2)/charmom
      vpnt(2) = vpnt(2) - 0.0003*binteg(4)/charmom
      vdir(2) = vdir(2) - 0.0003*binteg(5)/charmom      
      RETURN
      END
*--------1---------2---------3---------4---------5---------6---------712.......8

