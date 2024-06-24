#!/bin/bash
files=("$@")
for f in "${files[@]}"; do  
    if [[ "$f" =~ .T\\.txt$ ]]; then  
        cat "$f" >> T_cells.txt
    elif [[ "$f" =~ .G\\.txt$ ]]; then  
        cat "$f" >> G_cells.txt
    # Add other conditions similarly
    fi
done

# Rest of your processing here...
