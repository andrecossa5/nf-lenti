#!/bin/bash
files=($@)
for f in ${files[@]}; do  
    if [[ $f =~ .T.txt ]]; then
        cat "$f" >> T_cells.txt

    elif [[ $f =~ .G.txt ]]; then  
        cat "$f" >> G_cells.txt

    elif [[ $f =~ .A.txt ]]; then 
        cat "$f" >> A_cells.txt

    elif [[ $f =~ .C.txt ]]; then 
        cat "$f" >> C_cells.txt

    elif [[ $f =~ .coverage.txt ]]; then 
        cat "$f" >> coverage_cells.txt
    fi
done


