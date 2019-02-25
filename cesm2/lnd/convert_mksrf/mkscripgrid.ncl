;
;  mkscripgrid.ncl
;
;  Create SCRIP grid from LGM glacier dataset
;
;  Erik Kluzek
;
begin
  ; ===========================================================================================================
  ; Set a few constants needed later
  cdate  = systemfunc( "date -u +%y%m%d" );
  ldate  = systemfunc( "date -u" );
  printn = True;
  ;
  name = "glacier"
  nmask = "lgm21ka.wGIA"
  idate = "190225"
  res = "0.5x0.5"
  infilename = "mksrf_"+name+"_"+nmask+".c"+idate+".nc"
  outfilename = "SCRIPgrid_"+res+"_"+nmask+"_c"+cdate+".nc";

  system( "/bin/rm -f "+outfilename );
  print( "output file: "+outfilename );

  print( "Open: "+infilename )
  nc = addfile( infilename, "r" )
  lmask = nc->LANDMASK
  Mask2D = new( dimsizes(lmask), "integer", "No_FillValue" )
  Mask2D(:,:) = floattointeger(lmask)

  ;
  ; Create SCRIP grid file
  ;
  
  Opt = True
  Opt@Mask2D   = Mask2D
  ;Opt@GridCornerLat = latCorners
  ;Opt@GridCornerLon = lonCorners
  Opt@Title = "SCRIP grid file for "+name
  if (printn) then
    Opt@Debug = True
  end if
  curvilinear_to_SCRIP(outfilename, nc->LATIXY, nc->LONGXY, Opt)

  ;
  ; Add global attributes to file
  ;

  nc = addfile( outfilename, "w" );
  nc@history = ldate+": create using mkscripgrid.ncl";
  nc@Version = systemfunc( "git describe" );
  if ( printn )then
    print( "================================================================================================" );
    print( "Successfully created SCRIP grid file: "+outfilename);
  end if

  ; ===========================================================================================================

end
