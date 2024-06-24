#!/bin/bash
files=("$@")
for f in "${files[@]}"; do  
    echo $f
    if [[ "$f" =~ .T\\.txt$ ]]; then  
        cat "$f" >> T_cells.txt
    elif [[ "$f" =~ .G\\.txt$ ]]; then  
        cat "$f" >> G_cells.txt
    elif [[ "$f" =~ .A\\.txt$ ]]; then  
        cat "$f" >> A_cells.txt
    fi
done


