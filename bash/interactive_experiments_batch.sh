#!/bin/bash

while read type tau alt null grid times seed; do
    filename='job_'"$type"'_'"$tau"'_'"$alt"'_'"$null"'_'"$grid"'_'"$times"'.job'
    touch $filename
    echo "#!/bin/bash" >> $filename
    echo "#SBATCH -p RM" >> $filename
    echo "#SBATCH -N 1" >> $filename
    echo "#SBATCH -n 50" >> $filename
    echo "#SBATCH --mem 128GB" >> $filename
    echo "#SBATCH -t 12:00:00" >> $filename
    echo "#SBATCH --mail-type=ALL" >> $filename
    echo "Rscript interactive_testing_simulations.R "$type $tau $alt $null $grid $times $seed >> $filename
    sbatch $filename
done < bash/params.txt
