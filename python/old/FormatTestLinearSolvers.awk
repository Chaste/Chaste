BEGIN {USER=0}
IF $1=="RUN:" && USER==1 {print $2, $3; COUNTER=0; USER=1}
IF $1=="RUN:" && USER==0 {print $2, $3; COUNTER=0; USER=1}
COUNTER++ { }
IF COUNTER==13 && $2=="ERROR:" {print "Error\n"}
IF COUNTER==13 && $2!="ERROR:" {print $19"\n";}

#### if parallel comment out above two lines (if counter==..) and uncomment this
##IF $1=="avg:" {print $20"\n";}