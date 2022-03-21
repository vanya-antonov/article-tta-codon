
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
