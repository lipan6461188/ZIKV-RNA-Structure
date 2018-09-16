
### Download PsBL (http://olddriver.website/PsBL/) and compile to produce paris_background and paris_heatmap program

paris_background -in 59.sam -chr KU501215.1 -out /tmp/59.norm.matrix -method estimate
paris_heatmap -in /tmp/59.norm.matrix -file_type matrix -out_pdf /tmp/59.pdf -bins 400 -min 0 -max 3

paris_background -in 766.sam -chr AY632535.2 -out /tmp/766.norm.matrix -method estimate
paris_heatmap -in /tmp/766.norm.matrix -file_type matrix -out_pdf /tmp/766.pdf -bins 400 -min 0 -max 3

