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
 print ("Starting:  change_lakes_pct.ncl")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")

 dd= systemfunc("date -u +%y%m%d")

; ==============================================
; read files
; ==============================================

 basedir    = "/gpfs/fs1/work/sbradley/CESM2.0_SETUP/CLM/LAKE_DATA"
 ifile1     = basedir+"/data/lanwat_dataset_presentday.nc"
 ifile2     = basedir+"/data/icemask_LGM_10min.nc"
 ifile3      = basedir+"/data/icedammedlakes_LGM_10min.nc"
 ofile1     = "./lanwat_dataset_LGM.nc"
 ofile2     = "./lanwat_dataset_LGM_flip.nc"

;---Read present-day lake and wetland dataset
 in1 = addfile(ifile1,"r")
 lato        = in1->LATIXY
 lono        = in1->LONGXY
 pct_lake    = in1->PCT_LAKE
 pct_wetland = in1->PCT_WETLAND
 nlato       = 360
 nlono       = 720

;---Remove lakes and wetland when overlapping with ice sheets
 in2 = addfile(ifile2,"r")
 lati        = in2->lat
 loni        = in2->lon
 lgm_ice     = in2->sftgif   
 lgm_mask    = in2->sftlf

;---The two ice-dammed lakes to be added -> 1, 0 elsewhere 
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



;---Set pct_lakes = 0 if there is ice

    do j=0,nlato-1
       do i=0,nlono-1

          if(lgm_ice1d(j,i) .gt. 0) then
             pct_lake(j,i) = 0.
             pct_wetland(j,i) = 0.
          end if

       end do
    end do


;---Add the two ice-dammed lakes

    do j=0,nlato-1
       do i=0,nlono-1
       
          if(lgm_lakes1d(j,i) .gt. 0) then
             pct_lake(j,i) = 100.
             pct_wetland(j,i) = 0.
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
 fileAtt@create_date = systemfunc("date")
 fileAtt@src_code = "lake_dataset_lgm21ka.ncl"
 fileAtt@author = systemfunc("echo $USER")
 fileAtt@input1 = ifile1
 fileattdef(fout1,fileAtt)

 fout1->LANDMASK        = lgm_mask1d
 fout1->PCT_LAKE        = (/pct_lake/)
 fout1->PCT_WETLAND     = (/pct_wetland/)
 fout1->LONGXY          = (/lono/)
 fout1->LATIXY          = (/lato/)

;---ncks -A the other variables

 system("ncks -A -v EDGEE,EDGEN,EDGES,EDGEW,LAT,LON "+ifile1+" "+ofile1)
 
end

