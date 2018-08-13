call setmhdmagneticpressure(1.)
call setdiffisotropic(1)
call setdiffconductivity(1.)
rho1 = 1.
prs1 = 1.
gamma = 2.

brdx1 = -10.
brdx2 =  10.
if (resol == 0) then
  resol = int((brdx2-brdx1)/pspc1)
end if
pspc1 = (brdx2-brdx1)/resol
pspc2 = pspc1
brdy1 = -d2null*1.
brdy2 =  d2null*1.
brdz1 = -d3null*1.
brdz2 =  d3null*1.
bordersize = nb*pspc2
! call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store)
call createFixedBorders(store, ebc_all)
