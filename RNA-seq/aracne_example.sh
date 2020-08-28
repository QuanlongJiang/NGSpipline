#########################################################################

java -Xmx5G -jar aracne.jar -e TPM_log_qn.PC.symbol.txt -o ./ -t worm_TF.txt --pvalue 1E-8 --seed 2 --calculateThreshold

for i in {1..1000}; 
do java -Xmx5G -jar aracne.jar -e TPM_log_qn.PC.symbol.txt -o ./ -t worm_TF.txt --pvalue 1E-8 --seed $i --threads 30
done

java -Xmx5G -jar aracne.jar -o ./ --consolidate --nobonferroni

