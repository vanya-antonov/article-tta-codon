
# Data related to clustering
all_orgs.tsv
```
select c.id, c.name,
  (select t.num from cof_params t where t.name = 'num_orgs' and t.parent_id=c.id) AS num_orgs,
  (select t.num from cof_params t where t.name = 'tta__num_tta_genes' and t.parent_id=c.id) AS num_tta_genes,
  (select t.num from cof_params t where t.name = 'num_fs' and t.parent_id=c.id) AS num_fs,
  (select t.num from cof_params t where t.name = 'tta__num_tta_fs_genes' and t.parent_id=c.id) AS num_tta_fs_genes,
   c.descr
from cofs c
order by num_tta_genes desc;
```

all_cofs.tsv
```
select o.id, o.name, o.genus, o.kingdom, o.phylum,
  (select t.num from org_params t where t.name = 'num_annotated_genes' and t.parent_id=o.id) AS num_annotated_genes,
  (select t.num from org_params t where t.name = 'num_fshifts' and t.parent_id=o.id) AS num_fshifts,
  (select t.num from org_params t where t.name = 'tta__num_tta_genes' and t.parent_id=o.id) AS num_tta_genes,
  (select t.num from org_params t where t.name = 'tta__num_tta_fs_genes' and t.parent_id=o.id) AS num_tta_fs_genes,
  (select t.num from org_params t where t.name = 'virus_host_name' and t.parent_id=o.id) AS virus_host_name,
  (select t.num from org_params t where t.name = 'virus_host' and t.parent_id=o.id) AS virus_host
from orgs o
where name like 'Streptomyces%'
order by o.kingdom, o.name;
```

# Data related to adpA genes
```bash
# Co-occurrence of the TTA codons in the main ORF and in the uORF:
# adpA_with_ATG = 218
# adpA_with_ATG.uORF_TTA = 196
# adpA_with_ATG.internal_TTA = 184
# adpA_with_ATG.wo_internal_TTA  = 34
../src/grep_fasta.pl  --seq_re  'ATG.{6}$'  --pattern_str '.*'  adpA_with_ATG.around_ATG_start_51_9.fna  \
|  grep '>' | perl -pe 's/^>(\S+?)\:.*/$1/' \
> adpA_with_ATG.txt

../src/grep_fasta.pl  --seq_re  'TTA.{5,6}ATG.{6}$'  --pattern_str '.*'  adpA_with_ATG.around_ATG_start_51_9.fna  \
|  grep '>' | perl -pe 's/^>(\S+?)\:.*/$1/' \
> adpA_with_ATG.uORF_TTA.txt

../src/print_lines_with_unique_id_from_2_files.pl  --invert adpA_with_ATG.txt  adpA_with_internal_TTA.txt  > adpA_with_ATG.internal_TTA.txt
../src/print_lines_with_unique_id_from_2_files.pl           adpA_with_ATG.txt  adpA_with_internal_TTA.txt  > adpA_with_ATG.wo_internal_TTA.txt

# TTA in both the main ORF and uORF = 177
../src/print_lines_with_unique_id_from_2_files.pl  --invert  ../data/adpA_with_ATG.internal_TTA.txt  ../data/adpA_with_ATG.uORF_TTA.txt  >  adpA_with_ATG.internal_TTA.uORF_TTA.txt
```

