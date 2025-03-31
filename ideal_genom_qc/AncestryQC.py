import os
import psutil
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do, delete_temp_files
from ideal_genom_qc.get_references import Fetcher1000Genome, FetcherLDRegions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class ReferenceGenomicMerger():

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name:str, high_ld_regions: Path, reference_files: dict, built: str = '38', rename_snps: bool=False) -> None:

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_regions, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str):
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'")
        
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_regions.exists():
            raise FileNotFoundError("high_ld_regions does not exist")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
        self.high_ld_regions = high_ld_regions
        self.reference_files = reference_files
        self.renamed_snps    = rename_snps

        self.reference_AC_GT_filtered= None
        self.study_AC_GT_filtered    = None
        self.pruned_reference        = None
        self.pruned_study            = None
        self.reference_fixed_chr     = None
        self.reference_fixed_pos     = None
        self.reference_flipped       = None
        self.reference_cleaned       = None

        pass

    def execute_rename_snpid(self) -> None:
        
        if not self.renamed_snps:
            logger.info(f"STEP: Rename SNPs. `rename_snps` set to {self.renamed_snps}. Skipping renaming of SNPs in the study data")
            return
        else:
            logger.info(f"STEP: Rename SNPs. `rename` set to {self.renamed_snps}. Renaming SNPs in the study data to the format chr_pos_a1_a2")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        df_bim = pd.read_csv(self.input_path / (self.input_name+ '.bim'), sep="\t", header=None, dtype={0: str})
        df_bim.columns = ['chrom', 'snp', 'cm', 'pos', 'a1', 'a2']

        df_link = pd.DataFrame(columns=['old_id', 'new_id'])
        
        df_link['old_id'] = df_bim['snp']
        df_link['new_id'] = df_bim['chrom'].astype(str) + "_" + df_bim['pos'].astype(str) + "_" + df_bim['a1'].str[0] + "_" + df_bim['a2'].str[0]

        df_link.to_csv(self.input_path / "update_names.txt", sep="\t", header=False, index=False)
        logger.info(f"STEP: Renaming SNPs in the study data: created file {self.input_path / 'update_names.txt'}")

        # PLINK2 command
        plink2_cmd = f"plink2 --bfile {str(self.input_path / self.input_name)} --update-name {str(self.input_path / "update_names.txt")} --make-bed --out {str(self.input_path / (self.input_name+ '-renamed'))} --memory {memory} --threads {max_threads}"

        # Execute PLINK2 command
        shell_do(plink2_cmd, log=True)

        return

    def execute_filter_prob_snps(self)->None:

        logger.info("STEP: Filtering A->T and C->G SNPs from study and reference data.")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10
        
        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        if self.renamed_snps:

            # find A->T and C->G SNPs in study data
            filtered_study = self._filter_non_AT_or_GC_snps(target_bim=self.output_path / f"{self.input_name+'-renamed'}.bim", output_filename=self.input_name)
            logger.info("STEP: Filtering problematic SNPs from the study data: filtered study data")

        else:
            # find A->T and C->G SNPs in study data
            filtered_study = self._filter_non_AT_or_GC_snps(target_bim=self.input_path / f"{self.input_name}.bim", output_filename=self.input_name)
            logger.info("STEP: Filtering problematic SNPs from the study data: filtered study data")


        # find A->T and C->G SNPs in reference data
        filtered_reference = self._filter_non_AT_or_GC_snps(target_bim=self.reference_files['bim'], output_filename=self.reference_files['bim'].stem)
        logger.info("STEP: Filtering problematic SNPs from the study data: filtered reference data")

        self.reference_AC_GT_filtered= self.output_path / f"{self.reference_files['bim'].stem}-no_ac_gt_snps"
        self.study_AC_GT_filtered    = self.output_path / f"{self.input_name}-no_ac_gt_snps"

        with open(filtered_study, 'r') as f:
            logger.info(f"STEP: Filtering problematic SNPs from the study data: {len(f.readlines())} SNPs filtered")
        with open(filtered_reference, 'r') as f:
            logger.info(f"STEP: Filtering problematic SNPs from the reference data: {len(f.readlines())} SNPs filtered")

        if self.renamed_snps:
        
            # PLINK command: generate cleaned study data files
            plink_cmd1 = f"plink --bfile  {str(self.output_path / (self.input_name+'-renamed'))} --chr 1-22 --exclude {str(filtered_study)} --keep-allele-order --threads {max_threads} --make-bed --out {str(self.study_AC_GT_filtered)}"
        else:
            # PLINK command: generate cleaned study data files
            plink_cmd1 = f"plink --bfile  {str(self.input_path / self.input_name)} --chr 1-22 --exclude {str(filtered_study)} --keep-allele-order --threads {max_threads} --make-bed --out {str(self.study_AC_GT_filtered)}"

        # PLINK command: generate cleaned reference data files
        plink_cmd2 = f"plink --bfile  {self.reference_files['bim'].with_suffix('')} --biallelic-only strict --chr 1-22 --exclude {str(filtered_reference)} --keep-allele-order --allow-extra-chr --memory {memory} --threads {max_threads} --make-bed --out {str(self.reference_AC_GT_filtered)}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return
    
    def execute_ld_pruning(self, ind_pair:list) -> None:

        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        
        logger.info("STEP: LD-based pruning of study and reference data")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # PLINK command: generates prune.in and prune.out files from study data
        plink_cmd1 = f"plink --bfile {str(self.study_AC_GT_filtered)} --exclude range {self.high_ld_regions} --keep-allele-order --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --out {str(self.output_path / self.input_name)}"

        # PLINK command: prune study data and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {str(self.study_AC_GT_filtered)} --extract {str((self.output_path / self.input_name).with_suffix('.prune.in'))} --keep-allele-order --threads {max_threads} --make-bed --out {str((self.output_path / (self.input_name+'-pruned')))}"

        # PLINK command: generates a pruned reference data files
        plink_cmd3 = f"plink --bfile {str(self.reference_AC_GT_filtered)} --extract {str((self.output_path / self.input_name).with_suffix('.prune.in'))} --keep-allele-order --make-bed --threads {max_threads} --out {str((self.output_path / (self.reference_files['bim'].stem+'-pruned')))}"

        self.pruned_reference = self.output_path / (self.reference_files['bim'].stem+'-pruned')
        self.pruned_study = self.output_path / self.output_path / (self.input_name+'-pruned')

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return
    
    def execute_fix_chromosome_mismatch(self) -> dict:

        logger.info("STEP: Fixing chromosome mismatch between study data and reference panel")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # File paths
        study_bim = self.pruned_study.with_name(self.pruned_study.name + ".bim")
        reference_bim = self.pruned_reference.with_name(self.pruned_reference.name + ".bim")

        to_update_chr_file = self._find_chromosome_mismatch(study_bim, reference_bim)

        self.reference_fixed_chr = self.output_path / f"{self.reference_files['bim'].stem}-updateChr"

        with open(to_update_chr_file, 'r') as f:
            logger.info(f"STEP: Fixing chromosome mismatch between study data and reference panel: {len(f.readlines())} SNPs to update")

        # PLINK command
        plink_cmd = f"plink --bfile {self.pruned_reference} --update-chr {to_update_chr_file} 1 2 --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_fixed_chr}"

        # Execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_fix_possition_mismatch(self) -> None:

        logger.info("STEP: Fixing position mismatch between study data and reference panel")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # File paths
        study_bim = self.pruned_study.with_name(self.pruned_study.name + ".bim")
        reference_bim = self.pruned_reference.with_name(self.pruned_reference.name + ".bim")

        to_update_pos_file = self._find_position_mismatch(study_bim, reference_bim)

        self.reference_fixed_pos = self.output_path / f"{self.reference_files['bim'].stem}-updatePos"

        with open(to_update_pos_file, 'r') as f:
            logger.info(f"STEP: Fixing position mismatch between study data and reference panel: {len(f.readlines())} SNPs to update")

        # PLINK command
        plink_cmd = f"plink --bfile {self.reference_fixed_chr} --update-map {to_update_pos_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_fixed_pos}"

        # Execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_fix_allele_flip(self) -> None:

        logger.info("STEP: Allele flipping between study data and reference panel")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # File paths
        study_bim = self.pruned_study.with_name(self.pruned_study.name + ".bim")
        reference_bim = self.pruned_reference.with_name(self.pruned_reference.name + ".bim")

        to_flip_file = self.output_path / f"{self.reference_files['bim'].stem}.toFlip"
        self._find_allele_flip(study_bim, reference_bim, to_flip_file)

        self.reference_flipped = self.output_path / f"{self.reference_files['bim'].stem}-flipped"

        with open(to_flip_file, 'r') as f:
            logger.info(f"STEP: Allele flipping between study data and reference panel: {len(f.readlines())} SNPs to flip")

        # plink command
        plink_cmd = f"plink --bfile {self.reference_fixed_pos} --flip {to_flip_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_flipped}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return

    def execute_remove_mismatches(self) -> None:

        logger.info("STEP: Removing mismatched SNPs from reference data")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # File paths
        study_bim = self.pruned_study.with_name(self.pruned_study.name + ".bim")
        reference_bim = self.pruned_reference.with_name(self.pruned_reference.name + ".bim")

        mismatches_file = self.output_path / f"{self.reference_files['bim'].stem}.toRemove"
        self._find_allele_flip(study_bim, reference_bim, mismatches_file)

        self.reference_cleaned = self.output_path / f"{self.reference_files['bim'].stem}-cleaned"

        with open(mismatches_file, 'r') as f:
            logger.info(f"STEP: Removing mismatched SNPs from reference data: {len(f.readlines())} SNPs to remove")

        # plink command
        plink_cmd = f"plink --bfile {self.reference_flipped} --exclude {mismatches_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_cleaned}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_merge_data(self) -> None:

        logger.info("STEP: Merging study and reference data")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # plink command
        plink_cmd = f"plink --bfile {self.pruned_study} --bmerge {str(self.reference_cleaned.with_suffix('.bed'))} {str(self.reference_cleaned.with_suffix('.bim'))} {str(self.reference_cleaned.with_suffix('.fam'))} --keep-allele-order --threads {max_threads} --make-bed --out {self.output_path / (self.output_name+'-merged')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return

    def _filter_non_AT_or_GC_snps(self, target_bim: Path, output_filename: str) -> Path:

        df = pd.read_csv(
            target_bim, sep="\t", header=None, usecols=[1, 4, 5], names=["SNP", "A1", "A2"]
        )

        output_file = self.output_path / f"{output_filename}.ac_get_snps"

        filtered_snps = df[df[['A1', 'A2']].apply(lambda x: ''.join(sorted(x)) in {"AT", "TA", "GC", "CG"}, axis=1)]
        
        filtered_snps[["SNP"]].to_csv(output_file, index=False, header=False)

        return output_file
    
    def _find_chromosome_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Find mismatches where rsID is the same but chromosome differs
        mismatch_df = reference_df.merge(study_df[["chr", "rsid"]], on="rsid", suffixes=("_ref", "_study"))
        chromosome_mismatch_df = mismatch_df[mismatch_df["chr_ref"] != mismatch_df["chr_study"]]

        # Exclude chromosomes X and Y from updates
        mismatch_df = mismatch_df[~mismatch_df["chr_study"].astype(str).isin(["X", "Y"])]

        to_update_chr_file = self.output_path / "all_phase3.toUpdateChr"

        # Save the mismatch data to a file
        chromosome_mismatch_df[["chr_study", "rsid"]].to_csv(to_update_chr_file, sep="\t", header=False, index=False)

        return to_update_chr_file
    
    def _find_position_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary from file1 with column 2 as key and column 4 as value
        a = dict(zip(study_df['rsid'], study_df['pos_bp']))

        # Filter rows in reference_df where column 2 exists in 'a' and the values isn column 4 differ
        filtered = reference_df[reference_df['rsid'].map(a).notna() & (reference_df['pos_bp'] != reference_df['rsid'].map(a))]

        # Print the result to a file
        to_update_pos_file = self.output_path / f"{self.reference_files['bim'].stem}.toUpdatePos"
        filtered[['rsid', 'pos_bp']].to_csv(to_update_pos_file, sep="\t", header=False, index=False)

        return to_update_pos_file
    
    def _find_allele_flip(self, study_bim: Path, reference_bim: Path, output_filename: Path) -> None:

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary with the composite key from file1
        a = {f"{row['chr']}{row['rsid']}{row['pos_bp']}": f"{row['allele1']}{row['allele2']}" for _, row in study_df.iterrows()}

        # Filtering the rows in file2 based on the conditions
        filtered = reference_df[
            reference_df.apply(
                lambda row: (
                    f"{row['chr']}{row['rsid']}{row['pos_bp']}" in a and 
                    a[f"{row['chr']}{row['rsid']}{row['pos_bp']}"] not in {f"{row['allele1']}{row['allele2']}", f"{row['allele2']}{row['allele1']}"}
                ), axis=1
            )
        ]

        # Save the second column of filtered rows to a file
        filtered['rsid'].to_csv(output_filename, sep="\t", header=False, index=False)

        return
    
class GenomicOutlierAnalyzer:

    def __init__(self, input_path: Path, input_name: str, merged_file: Path, reference_tags: Path, output_path: Path, output_name: str) -> None:

        self.merged_file = merged_file
        self.reference_tags = reference_tags
        self.output_path= output_path
        self.output_name= output_name
        self.input_path = input_path
        self.input_name = input_name

        self.einvectors = None
        self.eigenvalues = None
        self.ancestry_fails = None
        self.population_tags = None

        pass

    def execute_pca(self, pca: int = 10, maf: float = 0.01) -> None:

        if not isinstance(pca, int):
            raise TypeError("pca should be an integer")
        if pca <= 0:
            raise ValueError("pca should be a positive integer")
        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be a float between 0 and 0.5")

        logger.info("STEP: Performing principal component decomposition")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # PLINK command: generate PCA for reference data
        plink_cmd = f"plink --bfile {str(self.merged_file)} --keep-allele-order --maf {maf} --out {str(self.output_path / (self.output_name+'-pca'))} --pca {pca} --memory {memory} --threads {max_threads}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        self.einvectors = self.output_path / (self.output_name+'-pca.eigenvec')
        self.eigenvalues = self.output_path / (self.output_name+'-pca.eigenval')

        return
    
    def find_ancestry_outliers(self, ref_threshold: float, stu_threshold: float, reference_pop: str, num_pcs: int = 2, fails_dir: Path = Path()) -> None:

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be a float")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be a float")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer")
        if num_pcs <= 0:
            raise ValueError("num_pcs should be a positive integer")
        if not isinstance(fails_dir, Path):
            raise TypeError("fails_dir should be a Path object")
        
        if not fails_dir.exists():
            logger.info("STEP: Identifying ancestry outliers: `fails_dir` does not exist.")
            logger.info(f"STEP: Identifying ancestry outliers: ancestry outliers will be saved in {self.output_path}")
            fails_dir = self.output_path
        
        logger.info("STEP: Identifying ancestry outliers")

        df_tags = pd.read_csv(self.reference_tags, sep="\t", usecols=['#IID', 'SuperPop'])
        df_tags['ID'] = '0'
        df_tags = df_tags[['ID', '#IID', 'SuperPop']]
        df_tags = df_tags.rename(columns={'ID': 'ID1', '#IID': 'ID2', 'SuperPop': 'SuperPop'})

        df = pd.read_csv(self.einvectors, sep=r"\s+",engine='python', header=None)
        logger.info("STEP: Identifying ancestry outliers: read eigenvec file")

        df = df[[0, 1]]
        df = df.rename(columns = {0: 'ID1', 1:'ID2'})

        df = pd.merge(df, df_tags, on=['ID1', 'ID2'], how='left')
        df['SuperPop'] = df['SuperPop'].fillna('StPop', inplace=False)

        df.to_csv((self.output_path / (self.output_name + 'pop_tags.csv')), sep='\t', index=False)

        self.population_tags = self.output_path / (self.output_name + 'pop_tags.csv')

        # filter samples who are ethnicity outliers
        ancestry_fails = self._find_pca_fails(
            output_path  = fails_dir,
            df_tags      = df,
            ref_threshold= ref_threshold,
            stu_threshold= stu_threshold,
            reference_pop= reference_pop,
            num_pcs      = num_pcs
        )

        self.ancestry_fails = ancestry_fails

        return
    
    def execute_drop_ancestry_outliers(self, output_dir: Path = Path()) -> None:

        logger.info("STEP: Dropping ancestry outliers from the study data")

        if not isinstance(output_dir, Path):
            raise TypeError("output_dir should be a Path object")
        
        if not output_dir.exists():
            logger.info("STEP: Dropping ancestry outliers from the study data: `output_dir` does not exist.")
            logger.info(f"STEP: Dropping ancestry outliers from the study data: ancestry outliers will be saved in {self.output_path}")
            output_dir = self.output_path

        with open(self.ancestry_fails, 'r') as f:
            logger.info(f"STEP: Dropping ancestry outliers from the study data: {len(f.readlines())} samples identified as ancestry outliers")

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {str(self.input_path / self.input_name)} --allow-no-sex --remove {str(self.ancestry_fails)} --make-bed --out {str(output_dir / (self.output_name+'-ancestry-cleaned'))}"

        # execute PLINK command
        shell_do(plink_cmd2, log=True)

        return
    
    def draw_pca_plot(self, plot_dir: Path = Path(), plot_name: str = 'pca_plot.jpeg') -> None:

        logger.info("STEP: Generating PCA plots")

        if not isinstance(plot_dir, Path):
            raise TypeError("plot_dir should be a Path object")
        if not isinstance(plot_name, str):
            raise TypeError("plot_name should be a string")
        
        if not plot_dir.exists():
            logger.info('STEP: Generating PCA plots: `plot_dir` does not exist.')
            logger.info(f'STEP: Generating PCA plots: pca plots will be saved in {self.output_path}')
            plot_dir = self.output_path

        # add population tags to pca output
        df_tags = pd.read_csv(self.population_tags, sep='\t')
        df_tags['ID1'] = df_tags['ID1'].astype(str)

        # load .eigenvec file and keep the first three principal components
        df_eigenvec = pd.read_csv(
            self.einvectors,
            header=None,
            sep   =r"\s+",
            engine='python'
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        # merge to get data with tagged populations
        df = pd.merge(df_eigenvec, df_tags, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        scatter_plot= sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop', ax=ax, marker='.', s=70)
        scatter_fig = scatter_plot.get_figure()
        scatter_fig.savefig(plot_dir / f'2D-{plot_name}', dpi=400)

        # generates a 3D scatter plot
        fig2= plt.figure()
        ax  = fig2.add_subplot(111, projection='3d')

        for s in df['SuperPop'].unique():
            ax.scatter(
                xs=df.pc_1[df.SuperPop==s],
                ys=df.pc_2[df.SuperPop==s],
                zs=df.pc_3[df.SuperPop==s], 
                label=s
            )
        ax.legend()
        plt.savefig(plot_dir / f'3D-{plot_name}', dpi=400)
        plt.close()

        return
    
    def _set_population_tags(self, psam_path: Path, study_fam_path: Path) -> pd.DataFrame:

        # Read population information from the .psam file
        df_psam = pd.read_csv(
            psam_path,
            sep='\t',
            usecols=['#IID', 'SuperPop']
        )

        # Set an ID column and rename columns for consistency
        df_psam['ID'] = 0
        df_psam = df_psam[['ID', '#IID', 'SuperPop']]
        df_psam.columns = ['ID1', 'ID2', 'SuperPop']

        # read individual IDs from the study .fam file
        df_fam = pd.read_csv(
            study_fam_path,
            sep=' ',
            header=None,
            index_col=False
        )

        # select relevant columns, assign a placeholder population tag, and rename columns
        df_fam = df_fam[df_fam.columns[:2]].copy()
        df_fam['SuperPop'] = 'StPop'
        df_fam.columns = ['ID1', 'ID2', 'SuperPop']

        # concatenate the two DataFrames to merge the information
        return pd.concat([df_fam, df_psam], axis=0)
    
    def _find_pca_fails(self, output_path: Path, df_tags: pd.DataFrame, ref_threshold: int, stu_threshold: int, reference_pop: str, num_pcs: str = 2) -> str:

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be an integer or float value")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be an integer or float value")
        if stu_threshold<=0:
            raise ValueError("stu_threshold should be a positive value")
        if ref_threshold<=0:
            raise ValueError("ref_threshold should be a positive value")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer value")
        if num_pcs<1:
            raise ValueError("num_pcs should be a positive integer")

        # filters reference subjects
        mask1 = (df_tags['SuperPop']==reference_pop)
        # filters subjects from study data
        mask2 = (df_tags['SuperPop']=='StPop')

        # generates two data frames with filtered subjects
        df_ref = df_tags[mask1].reset_index(drop=True)
        df_stu = df_tags[mask2].reset_index(drop=True)

        # read .eigenvec file
        df_eigenvec = pd.read_csv(
            self.einvectors,
            header=None,
            sep   =r"\s+",
            engine='python'
        )

        if num_pcs>df_eigenvec.shape[1]-2:
            raise ValueError("num_pcs should be less than or equal to the number of principal components in the .eigenvec file")
        
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:2+num_pcs]].copy()

        # renames columns for consistency
        new_col_names = []
        for k in range(2+num_pcs):
            if k<2:
                new_col_names.append(f"ID{k+1}")
            else:
                new_col_names.append(f"pc_{k-1}")
        df_eigenvec.columns = new_col_names

        # merge filtered subjects with its principal components
        df_ref = df_ref.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)
        df_stu = df_stu.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)

        # computes mean and standard deviation by columns in reference data
        mean_ref= df_ref[df_ref.columns[2:]].mean()
        std_ref = df_ref[df_ref.columns[2:]].std()

        # creates empty data frame
        outliers_1 = pd.DataFrame(columns=df_ref.columns)
        outliers_1[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `ref_threshold` std deviations from the reference mean
        for col in outliers_1.columns[2:]:
            outliers_1[col] = (np.abs(df_stu[col] - mean_ref[col]) > ref_threshold*std_ref[col])

        outliers_1['is_out'] = (np.sum(outliers_1.iloc[:,2:], axis=1) >0)

        df_1 = outliers_1[outliers_1['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        # computes mean and standard deviation by columns in study data
        mean_stu= df_stu[df_stu.columns[2:]].mean()
        std_stu = df_stu[df_stu.columns[2:]].std()

        # creates empty data frame
        outliers_2 = pd.DataFrame(columns=df_ref.columns)
        outliers_2[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `stu_threshold` std deviation from the study mean
        for col in outliers_2.columns[2:]:
            outliers_2[col] = (np.abs(df_stu[col] - mean_stu[col]) > stu_threshold*std_stu[col])

        outliers_2['is_out'] = (np.sum(outliers_2.iloc[:,2:], axis=1) >0)

        df_2 = outliers_2[outliers_2['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        df = pd.merge(df_1, df_2, on=['ID1', 'ID2'])

        ancestry_fails = output_path / (self.output_name + '_fail-ancestry-qc.txt')

        logger.info(f"STEP: Identifying ancestry outliers: {df.shape[0]} samples identified as ancestry outliers")

        # save samples considered as ethnicity outliers
        df.to_csv(
            ancestry_fails,
            header=None,
            index =False,
            sep   ='\t'
        )

        return ancestry_fails

class AncestryQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str, high_ld_file: Path, reference_files: dict = dict(), recompute_merge: bool = True, built: str = '38', rename_snps: bool = False) -> None:

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(recompute_merge, bool):
            raise TypeError("recompute_merge should be a boolean")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'")
        if not isinstance(rename_snps, bool):
            raise TypeError("rename_snps should be a boolean")
        
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions()
            ld_fetcher.get_ld_regions()

            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")
        
        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
        self.reference_files = reference_files
        self.high_ld_regions = high_ld_file
        self.recompute_merge = recompute_merge
        self.built = built
        self.rename_snps = rename_snps

        if not reference_files:

            logger.info(f"No reference files provided. Fetching 1000 Genomes reference data for built {self.built}")

            fetcher = Fetcher1000Genome(built=self.built)
            fetcher.get_1000genomes()
            fetcher.get_1000genomes_binaries()

            self.reference_files = {
                'bim': fetcher.bim_file,
                'bed': fetcher.bed_file,
                'fam': fetcher.fam_file,
                'psam': fetcher.psam_file
            }

        self.results_dir = self.output_path / 'ancestry_qc_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.merging_dir = self.results_dir / 'merging'
        self.merging_dir.mkdir(parents=True, exist_ok=True)

        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

        self.fail_samples_dir = self.results_dir / 'fail_samples'
        self.fail_samples_dir.mkdir(parents=True, exist_ok=True)

        self.clean_files = self.results_dir / 'clean_files'
        self.clean_files.mkdir(parents=True, exist_ok=True)

        pass

    def merge_reference_study(self, ind_pair: list = [50, 5, 0.2]) -> None:

        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not self.recompute_merge:
            logger.info("STEP: Merging study and reference data: recompute_merge is set to False. Skipping merging step")
            logger.info(f"STEP: Merging study and reference data: merged data is expected to be in {self.merging_dir}")
            return

        rgm = ReferenceGenomicMerger(
            input_path= self.input_path,
            input_name= self.input_name,
            output_path= self.merging_dir, 
            output_name= self.output_name,
            high_ld_regions =self.high_ld_regions, 
            reference_files = self.reference_files,
            rename_snps=self.rename_snps
        )

        #rgm.execute_rename_snpid()
        rgm.execute_filter_prob_snps()
        rgm.execute_ld_pruning(ind_pair=ind_pair)
        rgm.execute_fix_chromosome_mismatch()
        rgm.execute_fix_possition_mismatch()
        rgm.execute_fix_allele_flip()
        rgm.execute_remove_mismatches()
        rgm.execute_merge_data()

        return
    
    def _clean_merging_dir(self) -> None:

        for file in self.merging_dir.iterdir():
            if file.is_file() and '-merged' not in file.name:
                file.unlink()

        return
    
    def run_pca(self, ref_population: str, pca: int = 10, maf: float = 0.01, num_pca: int = 10, ref_threshold: float = 4, stu_threshold: float = 4) -> None:

        goa = GenomicOutlierAnalyzer(
            input_path= self.input_path, 
            input_name= self.input_name,
            merged_file= self.merging_dir / (self.output_name + '-merged'),
            reference_tags= self.reference_files['psam'],
            output_path= self.results_dir, 
            output_name= self.output_name
        )

        logger.info(f"STEP: Running PCA analysis: `ref_population` = {ref_population}")
        logger.info(f"STEP: Running PCA analysis: `pca` = {pca}")
        logger.info(f"STEP: Running PCA analysis: `maf` = {maf}")
        logger.info(f"STEP: Running PCA analysis: `num_pca` = {num_pca}")
        logger.info(f"STEP: Running PCA analysis: `ref_threshold` = {ref_threshold}")
        logger.info(f"STEP: Running PCA analysis: `stu_threshold` = {stu_threshold}")

        goa.execute_pca(pca=pca, maf=maf)
        goa.find_ancestry_outliers(
            ref_threshold=ref_threshold, 
            stu_threshold=stu_threshold, 
            reference_pop=ref_population, 
            num_pcs      =num_pca, 
            fails_dir    =self.fail_samples_dir
        )
        goa.execute_drop_ancestry_outliers(output_dir=self.clean_files)
        goa.draw_pca_plot(plot_dir=self.plots_dir)

        return
