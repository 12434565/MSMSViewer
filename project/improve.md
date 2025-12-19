1. why gaussian_sim_match(sigma=0.4) and ppm(<100) have same result
    1. it is important to know b or y which one is more similar. now stratigy is if b exist, skip y.
    2. i haven't code for mutiple b or y.
2. use class
    1. aaseq need class function like get aaseq.b and aaseq.y
    2. msms also 
        1. need class function like get the peak.text
        2. some other msms index(unnecessary)
3. init could be more clear
4. <S>html website(?)</S>
5. <S><5% no add b and y</S>
6. df.loc[df["mzs"] == mzs_i, "ppm"] change to df.loc[i, "ppm"]