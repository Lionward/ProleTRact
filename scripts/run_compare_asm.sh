tandemtwist_ccs_vcf="/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/Ashkenazi/TandemTwist/CCS/HG002.vcf.gz"

asm1="/confidential//tGenVar/vntr/output_maryam/tools/run_all_tools/output_not_complete/hgsvc/TandemTwist/asm/cut_reads_NA24385_h1.tsv.fasta"
asm2="/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output_not_complete/hgsvc/TandemTwist/asm/cut_reads_NA24385_h2.tsv.fasta"

out="output/comapre_tandemtwist_ccs_with_asm.csv"

tool="TandemTwist"
log="logs/comapre_tandemtwist_ccs_with_asm.log"

python run_seq_comparison.py --vcf $tandemtwist_ccs_vcf --asm1 $asm1 --asm2 $asm2 --out $out --tool $tool --log $log
