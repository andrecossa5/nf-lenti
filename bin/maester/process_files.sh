#!/bin/bash
files=($@)
for f in "${files[@]}"; do  
    if [[ $f =~ .T.txt ]]; then  
        cat "$f" >> T_cells.txt
    elif [[ $f =~ .G.txt ]]; then  
        cat "$f" >> G_cells.txt
    elif [[ $f =~ .A.txt ]]; then 
        echo $f 
        cat "$f" >> A_cells.txt
    fi
done


