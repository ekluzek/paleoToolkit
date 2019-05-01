load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; ----------------------------------------------------
; Read present-day lake dataset, LGM lakes, LGM ice,
; rewrite to netcdf with new LGM dataset
; author:  Michele Petrini
; ----------------------------------------------------

begin

 print ("=========================================")
 print ("Starting:  lake_dataset_21ka.ncl")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")

 dd= systemfunc("date -u +%y%m%d")

; ==============================================
; read files
; ==============================================

 basedir    = "/gpfs/fs1/work/sbradley/CESM2.0_SETUP/CLM/LAKE_DATA"
 ifile1     = "mksrf_LakePnDepth_0.5x0.5_simyr2004_csp1k.c190225.nc"
 ifile2     = basedir+"/../with_lakes/icemask_LGM_10min_190427.nc
 ifile3      = basedir+"/icedammedlakes_LGM_10min_190427.nc"
 ofile1     = "mksrf_LakePnDepth_0.5x0.5_simyr2004_csp1k_wLGMicedam_c"+dd+".nc"

;---Read present-day lake and wetland dataset
 in1 = addfile(ifile1,"r")
 lato        = in1->LATIXY
 lono        = in1->LONGXY
 pct_lake    = in1->PCT_LAKE
 nlato       = 360
 nlono       = 720

;---Remove lakes and wetland when overlapping with ice sheets
 in2 = addfile(ifile2,"r")
 lati        = in2->lat
 loni        = in2->lon
 lgm_ice     = in2->sftgif   
 lgm_mask    = in2->sftlf

;---The ice-dammed lakes to be added -> 1, 0 elsewhere 
 in3 = addfile(ifile3,"r")
 lgm_lakes   = in3->sftlf   

; ==============================================
; remap 10 min to 1 deg
; ==============================================

;---Set up options for regridding
    Opt                = True
    Opt@SrcGridLon     = loni   
    Opt@SrcGridLat     = lati
    Opt@SrcRegional    = False     
    Opt@DstGridLon     = lono
    Opt@DstGridLat     = lato
    Opt@DstRegional    = False            
    Opt@ForceOverwrite = True          
    Opt@PrintTimings   = True         
    Opt@InterpMethod   = "bilinear"  
    Opt@WgtFileName    = "./remap_weights.nc" 
;---Regrid lakes and ice sheet masks
    lgm_ice1d          = ESMF_regrid(lgm_ice,Opt)
    lgm_lakes1d        = ESMF_regrid_with_weights(lgm_lakes,basedir+"/remap_weights.nc",Opt)
    lgm_mask1d         = ESMF_regrid_with_weights(lgm_mask,basedir+"/remap_weights.nc",Opt) 
 
; ==============================================
;   Some loops/if-then cycles to update dataset
; ==============================================

;---Set landmask = land over the ice-dammed lakes
    do j=0,nlato-1
       do i=0,nlono-1

          if(lgm_lakes1d(j,i) .gt. 0) then
             lgm_mask1d(j,i) = 1.
          end if

       end do
    end do

;---Add the  ice-dammed lakes

    do j=0,nlato-1
       do i=0,nlono-1
       
          if(lgm_lakes1d(j,i) .gt. 0) then
             pct_lake(j,i) = 100.
          end if
          
       end do
    end do 

; ==============================================
;  write out new lake dataset 
; ==============================================

 system("/bin/rm -f " + ofile1)
 fout1 = addfile(ofile1,"c")

 fileAtt             = True
 fileAtt@title       = "New LGM lake dataset"
 fileAtt@create_date = systemfunc("date -u")
 fileAtt@src_code = "lake_dataset_lgm21ka.ncl"
 fileAtt@author = systemfunc("echo $USER")
 fileAtt@input1 = ifile1
 fileattdef(fout1,fileAtt)
 ;
 ; Define Dimensions
 ;
 dsizes = getfiledimsizes( in1 );
 ndims  = dimsizes(dsizes);
 ; Get dimension names
 dimnames = getfiledims( in1 );
 ;
 is_unlim = new( (/ ndims /), "logical" );
 is_unlim = False;
 do i = 0, ndims-1
    if ( dimnames(i) .eq. "ncl_scalar" )then
       dimnames(i) = "scalar";
    end if
 end do
 filedimdef( fout1, dimnames, dsizes, is_unlim );

 vars = (/ "LANDMASK", "PCT_LAKE", "LAKEDEPTH", "LONGXY", "LATIXY" /)
 do i = 0, dimsizes(vars)-1
    filevardef ( fout1, (/ vars(i) /), (/ "float" /),  (/ "lat", "lon" /) );
 end do

 fout1->LANDMASK        = doubletofloat(lgm_mask1d)
 fout1->PCT_LAKE        = pct_lake
 fout1->LAKEDEPTH       = in1->LAKEDEPTH
 fout1->LONGXY          = lono
 fout1->LATIXY          = lato

;---ncks -A the other variables

 system("ncks -A -v EDGEE,EDGEN,EDGES,EDGEW,LAT,LON "+ifile1+" "+ofile1)

 print( " ");
 print( "Successfully created lake data file" );
 
end

