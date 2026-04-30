import os
from dataclasses import dataclass
import pandas as pd, numpy as np
from typing import Literal, Optional
from . import utils

@dataclass
class MainResults:
    data: Optional[pd.DataFrame] = None
    sig: Optional[pd.DataFrame] = None
    sig_symbol: Optional[pd.DataFrame] = None

class ResultDataLoader:
    def __init__(self, folder: str, dataset:str, apoe_identifier: Literal['APOE4', 'APOE2'], 
    pathology_identifier: Literal['Aβ', 'AD'], summarize_oncol: str = 'Protein_id', 
    fig_foldrer: Optional[str] = None, sensitivity: Optional[str] = ''):
        self.dataset = dataset
        self.folder = folder
        self.summarize_oncol = summarize_oncol
        self.fig_folder = fig_foldrer

        self.apoe_identifier = apoe_identifier
        self.pathology_identifier = pathology_identifier

        self.neg_mark = 'Aβ-' if pathology_identifier == 'Aβ' else 'CU'
        self.pos_mark = 'Aβ+' if pathology_identifier == 'Aβ' else 'AD'
        self.sensitivity = sensitivity

        self._load_core_results()
        self._load_mediation()
        self._load_auxiliary()
        self._calculate_sets()
        self._build_id_maps()

    def _load(self, filename_base: str | list[str]) -> MainResults:
        """Try loading from one or more base filenames. Returns empty if none exist."""
        if isinstance(filename_base, str):
            filename_base = [filename_base]

        for base in filename_base:
            filepath = os.path.join(self.folder, f"{base}.csv")
            if os.path.exists(filepath):
                try:
                    data, sig, sig_symbol = utils.load_result_df(filepath, oncol=self.summarize_oncol, sig_col='p_adjusted')
                    return MainResults(data, sig, sig_symbol)
                except Exception as e:
                    print(f"[ERROR] Failed to load {filepath}: {e}")
                    continue
        print(f"[WARN] None of the files found: {[f'{self.folder}/{b}.csv' for b in filename_base]}")
        return MainResults(None, None, None)

    def _try_load_with_suffixes(self, base: str, suffixes: list[str]) -> MainResults:
        """Try base+suffix combinations like base_suffix.csv"""
        base_candidates = [f"{base}_{suffix}" for suffix in suffixes]
        return self._load(base_candidates)

    def _load_core_results(self):
        # Example: apoe2protein, apoe2protein_adjad or adjab
        self.apoe_result = self._load("apoe2protein")
        self.apoe_adj_result = self._try_load_with_suffixes("apoe2protein", ["adjad", "adjab"])

        self.pathology_result = self._load(["ad2protein", 'ab2protein'])
        self.pathology_adj_result = self._load(["ad2protein_adjapoe", 'ab2protein_adjapoe'])

    def _load_mediation(self):
        self.protein_mediation = self._load("protein_mediation")
        self.pathology_mediation = self._load(["ad_mediation", "ab_mediation"])

        if self.apoe_result.sig == None or len(self.apoe_result.sig) == 0 or self.protein_mediation.data is None or self.pathology_mediation.data is None:
            print("[WARN] Skipping _calculate_sets due to missing results.")
            self.path1_proteins = None
            self.path2_proteins = None
            self.path1_full = None
            self.path1_part = None
            self.path2_full = None
            self.path2_part = None
            return
        
        self.mediation_dict = utils.devide_bymediation(self.protein_mediation.data, self.pathology_mediation.data, self.summarize_oncol)
        self.path1_proteins = self.mediation_dict['path1_proteins']
        self.path2_proteins = self.mediation_dict['path2_proteins']
        self.path1_full = self.mediation_dict['path1_full']
        self.path1_part = self.mediation_dict['path1_part']
        self.path2_full = self.mediation_dict['path2_full']
        self.path2_part = self.mediation_dict['path2_part']

    def _load_auxiliary(self):
        self.ad_ine3e3 = self._load("ad_ine3e3")
        self.neg_allage = self._load(["neg_allage", "cu_allage"])
        self.neg_younger = self._load(["neg_younger", "cu_younger"])
        self.neg_older = self._load(["neg_older", "cu_older"])
        self.pos_allage = self._load(["pos_allage", "ad_allage"])

        self.apoe_age_interaction = self._load(["neg_allage_interaction", "cu_allage_interaction"])
        self.apoe_pathology_interaction = self._load(["apoe_ab_interaction", "apoe_ad_interaction"])

        self.protein_regulation = self._load(["protein_regulation"])

    def _calculate_sets(self):
        """Only calculate sets if all required inputs are available."""

        if any(result.sig is None for result in [
            self.apoe_result,
            self.apoe_result,
            self.apoe_adj_result,
            self.pathology_result,
            self.pathology_adj_result,
        ]) or len(self.apoe_result.sig) == 0 or self.path1_proteins is None or self.path2_proteins is None:
            print("[WARN] Skipping _calculate_sets due to missing results.")
            return

        # Before taking mediation into account
        self.apoe_solo = self.apoe_result.sig & self.apoe_adj_result.sig - self.pathology_result.sig - self.pathology_adj_result.sig
        self.pathology_solo = self.pathology_result.sig & self.pathology_adj_result.sig - self.apoe_result.sig - self.apoe_adj_result.sig
        self.shared = self.apoe_result.sig & self.apoe_adj_result.sig & self.pathology_result.sig & self.pathology_adj_result.sig
        self.undivided = self.apoe_result.sig - self.apoe_solo - self.path1_proteins - self.path2_proteins
        # After taking mediation into account
        self.apoe_solo_final = self.apoe_solo - self.path1_proteins - self.path2_proteins
        self.nonspecific = self.apoe_result.sig - self.path1_proteins - self.path2_proteins - self.apoe_solo_final

        assign_category = pd.DataFrame({self.summarize_oncol: self.apoe_result.data[self.summarize_oncol].tolist()})
        assign_category['Category'] = assign_category[self.summarize_oncol].apply(lambda x: f'{self.apoe_identifier}=protein=>{self.pathology_identifier} total' if x in self.path1_full else
                                                f'{self.apoe_identifier}=>protein=>{self.pathology_identifier} partital' if x in self.path1_part else
                                                f'{self.apoe_identifier}=>{self.pathology_identifier}=>protein total' if x in self.path2_full else
                                                f'{self.apoe_identifier}=>{self.pathology_identifier}=>protein partial' if x in self.path2_part else
                                                f'{self.apoe_identifier}-specific' if x in self.apoe_solo else
                                                f'{self.pathology_identifier}-specific' if x in self.pathology_solo else
                                                'Non-specific' if x in self.apoe_result.sig else np.nan)
        self.assigned_category = assign_category.copy()


    def _build_id_maps(self):
        df = self.apoe_result.data
        id_col = self.summarize_oncol

        oncol = self.summarize_oncol
        if 'symbol' in df.columns:
            sym = df['symbol'].fillna(df[oncol])
            self.id2symbol = dict(zip(df[id_col], sym))

        if 'label' in df.columns:
            lab = df['label'].fillna(df[oncol])
            self.id2label = dict(zip(df[id_col], lab))

        if 'apt_name' in df.columns:
            self.id2apt = dict(zip(df[id_col], df['apt_name']))


            

    def _extract_sup_dfs(self):
        df_dict = {
            'apoe_result':
            utils.sup_of_onedf(self.apoe_result.data, self.dataset, f'{self.apoe_identifier} associated proteins (Without {self.pathology_identifier} adjustment)'),
            'pathology_result':
            utils.sup_of_onedf(self.pathology_result.data, self.dataset, f'{self.pathology_identifier} associated proteins (Without {self.apoe_identifier} adjustment)'),
            'apoe_adj_result':
            utils.sup_of_onedf(self.apoe_adj_result.data, self.dataset, f'{self.apoe_identifier} associated proteins (With {self.pathology_identifier} adjustment)'),
            'pathology_adj_result':
            utils.sup_of_onedf(self.pathology_adj_result.data, self.dataset, f'{self.pathology_identifier} associated proteins (With {self.apoe_identifier} adjustment)'),
            
            'protein_mediation':
            utils.sup_of_onedf(self.protein_mediation.data, self.dataset, f'{self.apoe_identifier} => protein => {self.pathology_identifier}'),
            'pathology_mediation':
            utils.sup_of_onedf(self.pathology_mediation.data, self.dataset, f'{self.apoe_identifier} => {self.pathology_identifier} => protein '),

            'ad_ine3e3':
            utils.sup_of_onedf(self.ad_ine3e3.data, self.dataset, f'{self.pathology_identifier} associated proteins in e3e3 carriers'),

            'neg_allage':
            utils.sup_of_onedf(self.neg_allage.data, self.dataset, f'{self.apoe_identifier} associated proteins in {self.neg_mark}'),
            'neg_younger':
            utils.sup_of_onedf(self.neg_younger.data, self.dataset, f'{self.apoe_identifier} associated proteins in younger {self.neg_mark}'),
            'neg_older':
            utils.sup_of_onedf(self.neg_older.data, self.dataset, f'{self.apoe_identifier} associated proteins in older {self.neg_mark}'),

            'apoe_age_interaction':
            utils.sup_of_onedf(self.apoe_age_interaction.data, self.dataset, f'{self.apoe_identifier} * age effect on proteins in {self.neg_mark}'),

            'pos_allage':
            utils.sup_of_onedf(self.pos_allage.data, self.dataset, f'{self.apoe_identifier} associated proteins in {self.pos_mark}'),

            'apoe_pathology_interaction':
            utils.sup_of_onedf(self.apoe_pathology_interaction.data, self.dataset, f'{self.apoe_identifier} * {self.pathology_identifier} effect on proteins'),
        }

        on_cols = list(df_dict.values())[0][1]
        
        valid_dfs = [df for df,_a, in df_dict.values() if df is not None and not df.empty]

        if not valid_dfs or len(valid_dfs) == 0:
            print("[WARN] No valid DataFrames to merge.")
            return None

        merged = valid_dfs[0]
        

        for i, df in enumerate(valid_dfs[1:], start=1):
            merged = pd.merge(merged, df, on=on_cols, how='outer', suffixes=('', f'_x{i}'))

        first_col = merged.columns[0]
        header_mask = merged[first_col].isin([
            "OLINK proteins, corresponding genes and labels used for plots and annotation",
            "Somamers, corresponding genes and labels used for plots and annotation",
            'Gene expression from bulk RNA-seq'
        ])
        protein_mask = merged[first_col].isin(['Protein_id', 'apt_name', 'Gene expression'])

        header_rows = merged[header_mask]
        protein_rows = merged[protein_mask]
        other_rows = merged[~(header_mask | protein_mask)]

        merged = pd.concat([header_rows, protein_rows, other_rows], ignore_index=True)

        # Assaign every protein a category
        if self.dataset.split('-')[-1] != 'RNAseq' and self.assigned_category is not None:
            assigned_category = self.assigned_category.copy()
            if self.summarize_oncol not in merged.columns.tolist():
                if 'apt_name' in merged.columns.tolist() and self.id2apt is not None:
                    assigned_category['apt_name'] = assigned_category[self.summarize_oncol].apply(lambda x: self.id2apt[x])
                    assigned_category = dict(zip(assigned_category['apt_name'], assigned_category['Category']))
                    merged['Category'] = merged['apt_name'].apply(lambda x: assigned_category[x] if x in assigned_category.keys() else 'Category')
                else:
                    print("[WARN] cannot transfer back to summarize column according to the merged df")
            else:
                assigned_category = dict(zip(assigned_category[self.summarize_oncol], assigned_category['Category']))
                merged['Category'] = merged[self.summarize_oncol].apply(lambda x: assigned_category[x] if x in assigned_category.keys() else 'Category')
        else:
            print("[WARN] pass exporting category column, cause assigned_category is None")
        
        if 'Category' in merged.columns.tolist():
            top = merged.iloc[:2]
            rest = merged.iloc[2:].sort_values(by='Category', ascending=False)
            merged = pd.concat([top, rest], ignore_index=True)
            
            
        return merged

