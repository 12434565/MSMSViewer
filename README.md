# V2
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
   
# V3
## need to correct 
1. <s>df.loc[df["mzs"] == mzs_i, "ppm"] -> 不要用这种配对</s>
2. <S>y坐标变长,或者换成新的</S>
3. <S>先画黑的再画红的</S>
4. <S>用y盖b</S>
5. <s>ppm改成50</s>
6. <s>一个画布上两张图</S>

1. 加differences作为第三种比较情况()
---
# V4
1. change ppm and gs b and y compare part
2. add --help to chack the input
