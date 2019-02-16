;
; regrid_NH_GIA_on_top.ncl
;
; Regrid the Noth. Hemi. GIA data on top of the Peltier data in
;
begin
   ;
   ; Open the different files
   ;
   sdate    = systemfunc( "date +%y%m%d" );
   inglcfilename = "mksrf_glacier_lgm21ka.example.nc";
   print( "Open: "+inglcfilename );
   nc = addfile( inglcfilename, "r" )
   outglcfilename = "mksrf_glacier_lgm21ka.wNH_GIA.c"+sdate+".nc";
   print( "Open: "+outglcfilename );
   if ( systemfunc("test -f "+outglcfilename+"; echo $?" ) .eq. 0 )then
      system( "/bin/rm "+outglcfilename );
   end if
   nco = addfile( outglcfilename, "c" )
   ingiafilename = "/glade/u/home/mpetrini/forErik/LGM_TOPO/HUY3-SCEAN1D_21_cism4km.nc"
   print( "Open: "+ingiafilename );
   ncgia = addfile( ingiafilename, "r" )
   ;
   ; Get the elevation bins from the glacier dataset
   ;
   n_z = dimsizes(nc->BIN_CENTERS)
   z_edge = nc->BIN_EDGES
   nlat = dimsizes(nc->LAT)
   nlon = dimsizes(nc->LON)
   nlatgia = dimsizes(ncgia->thk(0,:,0))
   nlongia = dimsizes(ncgia->thk(0,0,:))
   nlat2 = nlat / 2
   print( nlat2 )
   ;
   ; Get the elevations of the 4km NH data
   ;
   elev = ncgia->thk(0,:,:) + ncgia->topg(0,:,:)
   pct_gec = new( (/ n_z, nlatgia, nlongia /), "float" );
   top_gec = pct_gec
   pct_gec = 0.0
   top_gec = 0.0
   nfill = 0
   nabov = 0
   do i = 0, nlatgia-1
   do j = 0, nlongia-1
      if ( elev(i,j) .ge. 0.0 )then
         nabov = nabov + 1
         do k = 0, n_z
            if ( (elev(i,j) .ge. z_edge(k)) .and. (elev(i,j) .lt. z_edge(k+1)) )then
               nfill = nfill + 1
               pct_gec(k,i,j) = 100.0
               top_gec(k,i,j) = elev(i,j)
               break
            end if
         end do
      end if
   end do
   end do
   if ( nfill .ne. nabov )then
       print( "nfill is NOT what it should be" )
       print( nfill )
       print( nabov )
       exit()
   end if
   do k = 0, n_z-1
      print( "k       :"+k+" z="+nc->BIN_CENTERS(k) )
      if ( sum(pct_gec(k,:,:)) .gt. 0.0 )then
         print( "avg elev:"+sum(top_gec(k,:,:)*pct_gec(k,:,:))/sum(pct_gec(k,:,:)) )
         print( "avg pct: "+avg(pct_gec(k,:,:)) )
      end if
   end do
   ;
   ; Regrid from the 4km NH to half degree NH
   ;
   print( "regrid" );
   NHpct_gec = rcm2rgrid ( ncgia->lat, ncgia->lon, pct_gec, nc->LAT(nlat2:), nc->LON, 1 );
   NHtop_gec = rcm2rgrid ( ncgia->lat, ncgia->lon, top_gec, nc->LAT(nlat2:), nc->LON, 1 );
   do k = 0, n_z-1
      print( "k       :"+k+" z="+nc->BIN_CENTERS(k) )
      if ( sum(NHpct_gec(k,:,:)) .gt. 0.0 )then
         print( "avg elev:"+sum(NHtop_gec(k,:,:)*NHpct_gec(k,:,:))/sum(NHpct_gec(k,:,:)) )
         print( "avg pct: "+avg(NHpct_gec(k,:,:)) )
      end if
   end do
   ;
   ; Add the NH data to the SH data
   ;
   print( "Add NH to SH" );
   new_pct_gec = (/ nc->PCT_GLC_ICESHEET /)
   print( dimsizes(new_pct_gec) );
   print( dimsizes(NHpct_gec(:,:,:)) )
   new_pct_gec(:,nlat2:,:) = NHpct_gec(:,:,:)
   new_top_gec = (/ nc->PCT_GLC_ICESHEET /)
   new_top_gec(:,nlat2:,:) = NHtop_gec(:,:,:)

   do i = 0, nlat2-1
   do j = 0, nlon-1
      do k = 0, n_z-1
         if ( ismissing(new_pct_gec(k,i,j)) )then
            new_top_gec(k,i,j) = (/ nc->BIN_CENTERS(k) /);
         end if
      end do
   end do
   end do
   ; Create pct_glc, landmask, bin_centers
   bin_centers = (/ nc->BIN_CENTERS /);
   do k = 0, n_z-1
      sump = doubletofloat(sum(new_pct_gec(k,:,:)))
      if ( sump .gt. 0.0 )then
         bin_centers(k) = (/ doubletofloat(sum(new_top_gec(k,:,:)*new_pct_gec(k,:,:)))/sump /)
      end if
   end do
   pct_glc  = (/ new_pct_gec(0,:,:) /)
   landmask = (/ new_pct_gec(0,:,:) /)
   do i = 0, nlat2
   do j = 0, nlon-1
      pct_glc(i,j) = (/ sum(new_pct_gec(:,i,j)) /)
      if ( new_pct_glc(i,j) .gt. 0.0d00 )then
         landmask(i,j) = 1.0
      end if
   end do
   end do
   ;
   ; Output the new data
   ;

   ;
   ; Define Dimensions
   ;
   print( "Define dims" );
   dsizes = getfiledimsizes( nc );
   ndims  = dimsizes(dsizes);
   ; Get dimension names
   dimnames = getfiledims( nc );
   ;
   is_unlim = new( (/ ndims /), "logical" );
   is_unlim = False;
   do i = 0, ndims-1
      if ( dimnames(i) .eq. "ncl_scalar" )then
         dimnames(i) = "scalar";
      end if
   end do
   filedimdef( nco, dimnames, dsizes, is_unlim );
   ;
   ; Define Variables
   ;
   print( "Define vars" );
   vars = getfilevarnames( nc );
   do i = 0, dimsizes(vars)-1
      print( "var: "+vars(i) );
      vardims = getfilevardims( nc, vars(i) );
      do j = 0, dimsizes(vardims)-1
         if ( vardims(j) .eq. "ncl_scalar" )then
            vardims(j) = "scalar";
         end if
      end do
      filevardef ( nco, (/ vars(i) /), (/ typeof(nc->$vars(i)$) /),  vardims );
      filevarattdef ( nco, (/ vars(i) /), nc->$vars(i)$ );
      delete( vardims );
   end do
   ;
   ; Global attributes
   ;
   print( "set global attributes" );
   fileattdef( nco, nc );
   date = systemfunc( "date" );
   desc = "Regrid using regrid_NH_GIA_on_top.ncl from: "+ingiafilename;
   history = ""
   if ( isvar(nc@history) )then
      history = nc@history;
   end if
   nco@history = date + " : " + desc + str_get_nl() + history;
   ;
   ; Copy data in
   ;
   print( "Set vars" );
   vars = getfilevarnames( nc );
   do i = 0, dimsizes(vars)-1
      print( "var: "+vars(i) );
      if ( vars(i) .eq. "PCT_GLACIER" )then
         nco->$vars(i)$ = (/ pct_glc /)
      else
        if ( vars(i) .eq. "PCT_GLC_GIC" )then
           nco->$vars(i)$ = 0.0d00
        else
           if ( vars(i) .eq. "PCT_GLC_ICESHEET" )then
              nco->$vars(i)$ = (/ new_pct_gec /)
           else
              if ( vars(i) .eq. "BIN_CENTERS" )then
                 nco->$vars(i)$ = (/ bin_centers /)
              else
                 if ( vars(i) .eq. "LANDMASK" )then
                    nco->$vars(i)$ = (/ landmask /)
                 ;
                 ; Otherwise just copy from the base file
                 ;
                 else
                    nco->$vars(i)$ = (/ nc->$vars(i)$ /)
                 end if
              end if
           end if
        end if
      end if

   end do
end
