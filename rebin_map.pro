;+
; Project     : SOHO-CDS
;
; Name        : REBIN_MAP
;
; Purpose     : Rebin an image map to new dimensions
;
; Category    : imaging
;
; Explanation : Rebin a map to user-specified dimensions and
;               compute new output pixel spacings
;
; Syntax      : gmap=rebin_map(map,gx,gy)
;
; Inputs      : MAP = image map structure
;               GX,GY = new dimensions
;
; Outputs     : GMAP = rebinned map
;
; History     : Written 22 August 1997, D. Zarro, ARC/GSFC
;               29 September 2008, Zarro (ADNET) 
;                - improved memory management
;               31 March 2012, Zarro (ADNET)
;                - made /interpolate default
;
; Contact     : dzarro@solar.stanford.edu
;-

function rebin_map,map,gx,gy,err=err,_extra=extra

;-- check inputs (valid map & dimensions)

if ~valid_map(map) or ~is_number(gx) then begin
 pr_syntax,'gmap=rebin_map(map,gx,gy)'
 if exist(map) then return,map else return,-1
endif
if not exist(gy) then gy=gx
ngx=nint(gx) & ngy=nint(gy)

if (ngx lt 2) or (ngy lt 2) then begin
 message,'both binning dimensions must be greater than 1',/cont
 return,map
endif

for i=0,n_elements(map)-1 do begin
 sz=size(map[i].data)
 nx=sz[1] & ny=sz[2]
 if (ngx eq nx) and (ngy eq ny) then begin
  message,'no rebinning necessary',/cont
  return,map
 endif
 tmap=rep_tag_value(map[i],congrid(map[i].data,ngx,ngy,/interp,_extra=extra),'data',/no_copy)
 tmap.dx=map[i].dx*nx/float(ngx)
 tmap.dy=map[i].dy*ny/float(ngy) 
 if i eq 0 then gmap=temporary(tmap) else $
  gmap=[temporary(gmap),temporary(tmap)]
endfor

return,gmap & end

