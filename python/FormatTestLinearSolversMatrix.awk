IF $1=="RUN:" {solver=$2; pc=$3; solvers[$2]=1; precons[$3]=1; tol=""; time=""; counter=0}
counter++{}
IF $(NF-1) == "!=" {tol=" "$(NF-2)" "$(NF-1)" "$(NF)}
IF $(NF) == "failed:" {time="Error"}
IF $1 == "Error" {time="Error"}
IF counter == 3 && time == "" {time = $1}
/KSP_DIV/ {time="Diverged"}
IF counter == 3 {results[solver"-"pc]=time" "tol}
END {printf("||  ||")
     for (pc in precons) printf("%s ||",pc)
     print ""
     for (solver in solvers) {
       printf("|| %s ||", solver)
       for (pc in precons) printf("%s ||", results[solver"-"pc]);
       print ""
     }
}
