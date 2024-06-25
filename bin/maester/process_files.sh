#!/bin/bash
files=($@)
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
for f in ${files[@]}; do  
    if [[ $f =~ .T.txt ]]; then 
        echo $f  
        cat "$f" >> T_cells.txt

    elif [[ $f =~ .G.txt ]]; then 
        echo $f  
        cat "$f" >> G_cells.txt

    elif [[ $f =~ .A.txt ]]; then 
        echo $f 
        cat "$f" >> A_cells.txt

    elif [[ $f =~ .C.txt ]]; then 
        echo $f 
        cat "$f" >> C_cells.txt

    elif [[ $f =~ .coverage.txt ]]; then 
        echo $f 
        cat "$f" >> coverage_cells.txt
    fi
done


