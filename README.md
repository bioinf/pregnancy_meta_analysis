# Meta-analysis of genome-wide association data from UK Biobank and FinnGen highlights risk loci for pregnancy complications

## Code

### Launch meta-analysis
To tun meta-analysis using METAL:
```./1.1_parallel.sh data/file_mapping.csv analysis_n```
And then using MTAG:
```./1.1_parallel_mtag.sh data/file_mapping.csv analysis_n```
(should be launched in the specific MATG's environment)

(Also, to postprocess all files, you can use: `./1.3_process_mtag_output.py`)

For this script you also need tsv-mapping file with 4 columns:
* `name_of_trait`
* `path_to_ukb_file`
* `path_to_finngen_file`
* `n_samples_in_finngen_file` (from manifest)

To draw Manhattan and Q-Q plots for all FinnGen traits:
```1.2_parallel_draw_fg.sh```

Auxilary scripts:
* `1.0.0_pipe_go_R6.py` - script for meta-analysis using metal tool for specific trait
* `1.0.0_pipe_go_R6_mtag.py` - the same analysis, but with MTAG tool (this script uses output of previous one).
* `1.0.1_DRAW.R` - script for drawing sketches of Manhattan and QQ plots


### Genetic correlations

1) Run the same meta-analysis on lots of pairs of traits.
2) The considered traits should be located in one directory (`PREG_DATA` in our case) and other in the directories, starting with `analysis_n_*`.


In case you re-launch calculations of correlations, first, clean-up old correlations:
```
find . -type f -name '*.cors' -delete
find . -type f -name '*.cors.full' -delete
find . -type f -name '*.labels' -delete
find . -type f -name '*.overlap' -delete
find . -type f -name '*.progress' -delete

```


Launch in parallel for lots of pairs of traits: ```2.1.0_parallel_prepare_traits_for_ldak.sh``` (here `2.0.0_prepare_traits_for_ldak.py` is aucilary script).

Launch for considered traits in `PREG_DATA` directory: ```./2.1.1_prepare_pregnancy_for_ldak.py```

3) Launch ldak itself:

`./2.2_parallel_launch_ldak.sh`

4) Assemble all correlations:

```
TEMP_T=("GEST_DIABETES1" "I9_HYPTENSPREG1" "O15_PRETERM1")

for t in "${TEMP_T[@]}" ; do for i in analysis_n* ; do ls ${i}/data/${t}*.cors 2> /dev/null ; cat ${i}/data/${t}*.cors 2> /dev/null | grep Cor_All | awk '{  if ($2 > 2.57*$3)  print }' | grep -v nan ; done | grep -B 1 Cor_All > data/cor_${t}.txt ; done

for t in "${TEMP_T[@]}" ; do for i in analysis_n* ; do ls ${i}/data/${t}*.cors 2> /dev/null ; cat ${i}/data/${t}*.cors 2> /dev/null | grep Cor_All | awk '{ print }' | grep -v nan ; done | grep -B 1 Cor_All > data/cor_full_${t}.txt ; done

```

As a result we will have `cor_full_*.txt` with all not-na genetic correlations for specific trait and `cor_*.txt` files with filtered by significance genetic correlations.

5) Draw genetic correlation plot:

* Launch `2.3.1_make_table_for_r.ipynb` to prepare table
* Launch `2.3.2_draw_gen_cor.R` - to draw the plot

### Select and annotate top snps

`3.1_making_top_snp_table_FG.ipynb` and `3.2_making_top_snp_table_META.ipynb` -- selecting and annotation of top SNPs for FinnGen-only and meta-analysis respectively.


### Final Manhattan and Q-Q plots

`4_final_mh_qq.R` - for drawing final versions of Q-Q and Manhattan plots.


## Images
All images are located in `img` directory. 
* `img/final_gen_cor.pdf` - forest plot with genetic correlations
* `img/QQplot.pval__*.pdf` - Q-Q plots of significant traits
    * `img/QQplot.pval__FG_*.pdf` - for FinnGen data only.
    * `img/QQplot.pval__MET_*.pdf` - for meta-analysis data only.
* `img/Rectangular-Manhattan..pval__*.pdf` - Manhattan plots of significant traits
    * `img/Rectangular-Manhattan..pval__FG_*.pdf` - for FinnGen data only.
    * `img/Rectangular-Manhattan..pval__MET_*.pdf` - for meta-analysis data only.
    
    
## Data

All data is located in `data` directory:
* Selected summary statistics:
    * `data/f_special/` - directory with filtered FinnGen GWAS summary statistics (only selected as significant).
    * `data/f_special/` - directory with meta-analysis summary statistics (only selected as significant).
* Genetic correlations:
    * `data/cor_full_<trait>.txt` - file with all genetic correlations for selected traits
    * `data/cor_<trait>.txt` - file with significant genetic correlations for selected traits
    * `data/feature.csv` - annotated table with significant genetic correlations 
* Annotated SNPs:
    * `data/finn_top.csv` - significant annotated summstats from FinnGen GWAS.
    * `data/finn_top_short.csv` - significant and filtered (selected 1 per loci) annotated summstats from FinnGen GWAS.
    * `data/meta_top.csv`- significant annotated summstats from meta-analysis.
    * `data/meta_top_short.csv` - significant and filtered (selected 1 per loci) annotated summstats from meta-analysis.
* Other:
    * `data/file_mapping` - mapping of selected 24 traits files and N_samples for finngen
    *  All summary statistics can be found [here](https://drive.google.com/drive/folders/1J0KvO_G50yOtg8wxNKAPZryUCs3OLH9w?usp=sharing):
        * `maf_fg_*.tsv` - FinnGen summary statistics filtered by MAF.
        * `extended_*.TBL` - summary statistics from meta-analysis by METAL.




