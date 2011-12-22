BEGIN {USER=0; STATUS=0}
IF $1=="RUN:" && USER==1 {print; print$0; STATUS=1; USER=1}
IF $1=="RUN:" && USER==0 {print $0; STATUS=1; USER=1}
IF STATUS==1 && $1=="Passed" {print $0; STATUS=0}
IF STATUS==1 && match($0, "Error: Expected") {print $0; STATUS=0}
IF STATUS==1 && match($0, "KSP_") {print $0; STATUS=0}
IF $1=="user" && STATUS==0 {print $2; USER=0}
IF $1=="user" && STATUS==1 {print "Error (no status found)"; print $2; USER=0}