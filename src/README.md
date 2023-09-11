# Statistics for the manuscript
"*In 177 (96%) of the 184 Streptomyces adpA genes that have a UUA in the linker region of the main ORF, an upstream ORF has a UUA 5 or 6 nt upstream the main ORF start codon*"
```bash
# 177
./print_lines_with_unique_id_from_2_files.pl  --invert  ../data/adpA_with_ATG.internal_TTA.txt  ../data/adpA_with_ATG.uORF_TTA.txt  |  wc -l 

# 184
wc -l  ../data/adpA_with_ATG.internal_TTA.txt
```

"*Interestingly, in the other 34  \adpA sequences analyzed, \MARK{23\% of the total}, where there is no UUA in the main ORF ‘linker’, the proportion with a corresponding UUA-containing upstream ORF is much lower, \red{56\% (i.e. 19 of the 34)}.*"

```bash
# 34
wc -l  ../data/adpA_with_ATG.wo_internal_TTA.txt

# Total number of adpA with ATG = 218  =>  34 / 218 = 15%
wc -l  ../data/adpA_with_ATG.txt

# 19  =>  19 / 34 = 56%
print_lines_with_unique_id_from_2_files.pl  --invert  ../data/adpA_with_ATG.wo_internal_TTA.txt  ../data/adpA_with_ATG.uORF_TTA.txt  |  wc -l 
```

"*For the 177 sequences with a main ORF UUA the upstream ORF UUA is located 5nt (65 cases) or 6nt (11 cases)  5’ of the main ORF start*"
```bash
# 177
grep '>' ../data/adpA_with_ATG.internal_TTA.uORF_TTA.fasta  | wc -l

# 153
grep_fasta.pl  --seq_re  'TTA.{5}ATG.{6}$'  --pattern_str '.*'  ../data/adpA_with_ATG.internal_TTA.uORF_TTA.fasta  |  grep '>' | wc -l

# 24
grep_fasta.pl  --seq_re  'TTA.{6}ATG.{6}$'  --pattern_str '.*'  ../data/adpA_with_ATG.internal_TTA.uORF_TTA.fasta  |  grep '>' | wc -l
```