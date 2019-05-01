;
; regrid_GIA_glc.ncl
;
; Regrid the GIA glacier data
;
begin
   ;
   ; Open the different files
   ;
   sdate    = systemfunc( "date -u +%y%m%d" );
   inglcfilename = "mksrf_glacier_0.5x0.5_simyr2000.c190225.nc"
   print( "Open: "+inglcfilename );
   nc = addfile( inglcfilename, "r" )
   outglcfilename = "mksrf_glacier_lgm21ka.wGIA.wLGMicedam.c"+sdate+".nc";
   print( "Open: "+outglcfilename );
   if ( systemfunc("test -f "+outglcfilename+"; echo $?" ) .eq. 0 )then
      system( "/bin/rm "+outglcfilename );
   end if
   nco = addfile( outglcfilename, "c" )
   ingiafilename = "/gpfs/fs1/work/sbradley/CESM2.0_SETUP/GIA_CAM-CISM/CAM-CLM/Masks_USGS-huy3-S1D_21_10min_2019.nc"
   print( "Open: "+ingiafilename );
   ncgia = addfile( ingiafilename, "r" )
   ingiaofilename = "/gpfs/fs1/work/sbradley/CESM2.0_SETUP/GIA_CAM-CISM/CAM-CLM/Orog_USGS-huy3-S1D_21_10min_2019.nc"
   print( "Open: "+ingiaofilename );
   ncgiao = addfile( ingiaofilename, "r" )
   basedir    = "/gpfs/fs1/work/sbradley/CESM2.0_SETUP/CLM/LAKE_DATA"
   icedamlake  = basedir+"/icedammedlakes_LGM_10min_190427.nc"
   print( "Open: "+icedamlake );
   ncicedam = addfile(icedamlake, "r" );
   ;
   ; Get the elevation bins from the glacier dataset
   ;
   n_z = dimsizes(nc->BIN_CENTERS)
   z_edge = nc->BIN_EDGES
   nlat = dimsizes(nc->LAT)
   nlon = dimsizes(nc->LON)
   nlatgia = dimsizes(ncgia->sftgif(:,0))
   nlongia = dimsizes(ncgia->sftgif(0,:))
   ;
   ; Get the elevations of the 4km data
   ;
   elev = (/ ncgiao->RELIEF(:,:) /)
   pct_gec = new( (/ n_z, nlatgia, nlongia /), "double" );
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
               pct_gec(k,i,j) = ncgia->sftgif(i,j)*100.0d00 ; Ice fraction
               top_gec(k,i,j) = doubletofloat( elev(i,j) )                ; Elevation
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
   ; Regrid from the 4km to half degree
   ;
   print( "regrid" );
   print( "pct_gec" );
   printVarSummary( pct_gec );
   printVarSummary( nc->LON );
   printVarSummary( nc->LAT );
   HFpct_gec = area_conserve_remap ( ncgia->lon, ncgia->lat, pct_gec, nc->LON, nc->LAT, False );
   print( "top_gec"  );
   HFtop_gec = area_conserve_remap ( ncgia->lon, ncgia->lat, top_gec, nc->LON, nc->LAT, False );
   do k = 0, n_z-1
      print( "k       :"+k+" z="+nc->BIN_CENTERS(k) )
      if ( sum(HFpct_gec(k,:,:)) .gt. 0.0 )then
         print( "avg elev:"+sum(HFtop_gec(k,:,:)*HFpct_gec(k,:,:))/sum(HFpct_gec(k,:,:)) )
         print( "avg pct: "+avg(HFpct_gec(k,:,:)) )
      end if
   end do
   ; Create pct_glc, landmask, bin_centers
   bin_centers = (/ nc->BIN_CENTERS /);
   do k = 0, n_z-1
      sump = doubletofloat(sum(HFpct_gec(k,:,:)))
      if ( sump .gt. 0.0 )then
         bin_centers(k) = (/ doubletofloat(sum(HFtop_gec(k,:,:)*HFpct_gec(k,:,:)))/sump /)
      end if
   end do
   pct_glc  = (/ HFpct_gec(0,:,:) /)
   landmask = area_conserve_remap ( ncgia->lon, ncgia->lat, doubletofloat( ncgia->sftlf ), nc->LON, nc->LAT, False );
   do i = 0, nlat-1
   do j = 0, nlon-1
      pct_glc(i,j) = (/ sum(HFpct_gec(:,i,j)) /)
      if ( pct_glc(i,j) .gt. 0.0d00 )then
         landmask(i,j) = 1.0
      end if
      if ( ncicedam->sftlf(i,j) .gt. 0.0d00 )then
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
      print( "var: "+vars(i)+": "+typeof(nc->$vars(i)$) );
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
   date = systemfunc( "date -u" );
   desc = "Regrid using regrid_GIA_glc.ncl from: "+ingiafilename;
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
         nco->$vars(i)$ = doubletofloat(pct_glc)
      else
        if ( vars(i) .eq. "PCT_GLC_GIC" )then
           nco->$vars(i)$ = 0.0
        else
           if ( vars(i) .eq. "PCT_GLC_ICESHEET" )then
              nco->$vars(i)$ = doubletofloat( HFpct_gec )
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
   print( "" );
   print( "Successfully, regridded file" );
end
