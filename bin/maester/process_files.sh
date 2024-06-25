#!/bin/bash
files=($1)
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo "it's started"
echo ${files}
for f in ${files[@]}; do  
    echo $f
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


