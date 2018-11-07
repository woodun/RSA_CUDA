

echo "" >> runtime.txt


./main nodiv_nopre1.txt nondivpairs_nopre1.txt 10000 >> runtime.txt
./main div_nopre1.txt divpairs_nopre1.txt 10000 >> runtime.txt

./main nodiv_nopre0.txt nondivpairs_nopre0.txt 10000 >> runtime.txt
./main div_nopre0.txt divpairs_nopre0.txt 10000 >> runtime.txt


echo "" >> runtime.txt

./main nodiv_pre1.txt nondivpairs_pre1.txt 10000 >> runtime.txt
./main div_pre1.txt divpairs_pre1.txt 10000 >> runtime.txt

./main nodiv_pre0.txt nondivpairs_pre0.txt 10000 >> runtime.txt
./main div_pre0.txt divpairs_pre0.txt 10000 >> runtime.txt

