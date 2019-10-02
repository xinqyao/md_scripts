# this is a script to fit every frame to the first
# given a molecule id (default is top)

proc fitframes {{mol top} {sel "name CA"}} {
   set ref [atomselect $mol $sel frame 0]
   set fr  [atomselect $mol $sel]
   set mov [atomselect $mol "all"]
   set nsteps [molinfo $mol get numframes]
   for {set step 1} {$step < $nsteps} {incr step} {
       $fr frame $step
       $mov frame $step
       set transmat [measure fit $fr $ref]
       $mov move $transmat
   }
}
