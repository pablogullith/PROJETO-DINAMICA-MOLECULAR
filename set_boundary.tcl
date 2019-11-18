set n [molinfo top get numframes]
for {set i 0} {$i<$n} {incr i} {
    animate goto $i
    pbc set "10.0 10.0 10.0"
}
pbc box
