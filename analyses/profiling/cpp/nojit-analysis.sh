b1=nodefine
b2=1define
b3=12-define
b4=12-define-Histo2D
b5=12-define-6-Histo2D

common=distributions-allhistos-nojit

for b in $b1 $b2 $b3 $b4 $b5
do
    igprof -d -pp -z -o igprof.pp.gz "./$common-$b" distill_DS1_d45b_56t.root  
    igprof-analyse -d -v -g igprof.pp.gz >& igreport_perf-${b}.res
    igprof-analyse --sqlite -d -v -g igprof.pp.gz | sqlite3 ig-${b}.sql3    
done
