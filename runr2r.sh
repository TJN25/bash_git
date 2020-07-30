#!/bin/sh

esl-reformat pfam $1 | egrep -v ^"#=GS|#=GF" > ./.stk0
r2r --GSC-weighted-consensus ./.stk0 ./.stk 3 0.95 0.75 0.60 4 0.95 0.75 0.60 0.4 0.05
r2r ./.stk $1.svg
r2r ./.stk $1.pdf
