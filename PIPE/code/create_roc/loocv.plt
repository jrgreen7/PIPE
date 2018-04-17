set terminal post color
set key right bottom Left
set output "loocv.eps"
set ylabel "Sensitivity (True Positive Rate)"
set xlabel "1-Specificity (False Positive Rate)"
set mytics 2 
set mxtics 2

set style line 1 lt 1 lw 1 pt 1 lc rgb "red" 
set style line 2 lt 1 lw 1 pt 1 lc rgb "blue" 

plot    "PIPE_score/roc_curves/n1.roc" using 1:2 smooth bezier with lines ls 1 title "Traditional PIPE score 1",\
        "PIPE_score/roc_curves/n2.roc" using 1:2 smooth bezier with lines ls 1 title "Traditional PIPE score 2",\
        "PIPE_score/roc_curves/n3.roc" using 1:2 smooth bezier with lines ls 1 title "Traditional PIPE score 3",\
        "SW_score/roc_curves/n1.roc" using 1:2 smooth bezier with lines ls 2 title "Similarity-Weighted score 1",\
        "SW_score/roc_curves/n2.roc" using 1:2 smooth bezier with lines ls 2 title "Similarity-Weighted score 2",\
        "SW_score/roc_curves/n3.roc" using 1:2 smooth bezier with lines ls 2 title "Similarity-Weighted score 3"
