# Alternative way of getting/setting view points for all molecules
proc get_vp {} {
  global M 
  set molid 0
  set M [list [molinfo $molid get center_matrix] [molinfo $molid get rotate_matrix] [molinfo $molid get scale_matrix] [molinfo $molid get global_matrix ] ]
#  return $M
}

proc set_vp {} {
  global M 
  foreach mol [molinfo list] {
    molinfo $mol set rotate_matrix   [lindex $M 1]
    molinfo $mol set center_matrix   [lindex $M 0]
    molinfo $mol set scale_matrix    [lindex $M 2]
    molinfo $mol set global_matrix   [lindex $M 3]
  }
}
