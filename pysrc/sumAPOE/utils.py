import pandas as pd

def load_result_df(file_path: str, oncol='Protein_id', sig_col='p_adjusted') -> tuple[pd.DataFrame, set]:
    """
    Load a results DataFrame from a CSV file.
    Modify on some cols

    Parameters:
    file_path (str): The path to the CSV file containing the results.
    
    Returns:
    pd.DataFrame: The loaded DataFrame.
    """
    df = pd.read_csv(file_path)

    if 'Protein_id' in df.columns.tolist() and 'CE__ANML' in df['Protein_id'].tolist():
        df = df[df['Protein_id'] != 'CE__ANML']

    if 'p' in df.columns:
        df['significance'] = df.apply(lambda row: '*+' if (row['p']<0.05) & (row['p_adjusted']<0.05) else('*' if row['p']<0.05 else ''), axis=1)
    elif 'acme_p' in df.columns:
        df['significance'] = df.apply(lambda row: '*+' if (row['acme_p']<0.05) & (row['p_adjusted']<0.05) else('*' if row['acme_p']<0.05 else ''), axis=1)

    sig = set(df[df[sig_col] < 0.05][oncol].tolist())

    if 'symbol' in df.columns:
        sig_symbol = set(df[df[sig_col] < 0.05]['symbol'].tolist())
    else:
        sig_symbol = None
    return df,sig, sig_symbol

def devide_bymediation(protein_mediation: pd.DataFrame, ab_mediation: pd.DataFrame, oncol: str) -> dict:
    """
    Divide proteins by mediation results.
    """
    protein_mediation_sig = set(protein_mediation[protein_mediation['p_adjusted'] < 0.05][oncol])
    ab_mediation_sig = set(ab_mediation[ab_mediation['p_adjusted'] < 0.05][oncol])

    tmp = pd.merge(protein_mediation, ab_mediation, on=oncol, suffixes=('_protein', '_ab'))
    tmp['path1'] = abs(tmp['prop_m_protein']) > abs(tmp['prop_m_ab'])
    tmp['path2'] = abs(tmp['prop_m_protein']) < abs(tmp['prop_m_ab'])
    path1 = set(tmp[tmp['path1'] == True][oncol].tolist())
    path2 = set(tmp[tmp['path2'] == True][oncol].tolist())

    path1_proteins = set()
    path2_proteins = set()
    for x in protein_mediation_sig | ab_mediation_sig:
        if x in protein_mediation_sig and x in ab_mediation_sig:
            if x in path1:
                path1_proteins.add(x)
            elif x in path2:
                path2_proteins.add(x)
            else:
                print(x)
        elif x in protein_mediation_sig:
            path1_proteins.add(x)
        elif x in ab_mediation_sig:
            path2_proteins.add(x)
            
    ab_full_mediation = ab_mediation[(ab_mediation['p_adjusted']<0.05) & (ab_mediation['ade_p_adjusted']>=0.05)]
    ab_partial_mediation = ab_mediation[(ab_mediation['p_adjusted']<0.05) & (ab_mediation['ade_p_adjusted']<0.05)]

    protein_full_mediation = protein_mediation[(protein_mediation['p_adjusted']<0.05) & (protein_mediation['ade_p_adjusted']>=0.05)]
    protein_partial_mediation = protein_mediation[(protein_mediation['p_adjusted']<0.05) & (protein_mediation['ade_p_adjusted']<0.05)]

    path1_full = set(protein_full_mediation[protein_full_mediation[oncol].isin(path1_proteins)][oncol].tolist())
    path1_part = set(protein_partial_mediation[protein_partial_mediation[oncol].isin(path1_proteins)][oncol].tolist())
    path2_full = set(ab_full_mediation[ab_full_mediation[oncol].isin(path2_proteins)][oncol].tolist())
    path2_part = set(ab_partial_mediation[ab_partial_mediation[oncol].isin(path2_proteins)][oncol].tolist())

    return {
        'path1_full': path1_full,
        'path1_part': path1_part,
        'path2_full': path2_full,
        'path2_part': path2_part,
        'path1_proteins': path1_proteins,
        'path2_proteins': path2_proteins
    }

def extract_mainmodel_results(pc1, Results) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ Extracts main model results for a given list of proteins (pc1) from the Results object.
    Returns two DataFrames for heatmap plot: models and annotations.
    Args:
        pc1 (list[str]): List of protein symbols to filter results.
        Results (ResultDataLoader): An instance of ResultDataLoader containing the analysis results.
        Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames for heatmap plot:
            - models: Contains standardized beta values for the specified proteins in each model.
            - annotations: Contains p-values for the specified proteins in each model.
    """
    models_df = pd.DataFrame({Results.summarize_oncol: list(pc1)})
    
    def extract_and_merge(source_result, beta_col: str, p_col: str):
        if source_result.data is None:
            return pd.DataFrame(columns=[Results.summarize_oncol, beta_col, p_col])
        
        tmp = source_result.data[source_result.data[Results.summarize_oncol].isin(pc1)][
            [Results.summarize_oncol, 'standardized_beta', 'significance']
        ].copy()
        tmp.rename(columns={'standardized_beta': beta_col, 'significance': p_col}, inplace=True)
        return tmp

    tmp = extract_and_merge(Results.apoe_result, 'β2', 'p2')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.pathology_result, 'β4', 'p4')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.apoe_adj_result, 'β6', 'p6')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.pathology_adj_result, 'β7', 'p7')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.neg_allage, f'β({Results.neg_mark})', f'p({Results.neg_mark})')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.apoe_age_interaction, f'β({Results.apoe_identifier}*age in {Results.neg_mark})', 
                            f'p({Results.apoe_identifier}*age in {Results.neg_mark})')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.pos_allage, f'β({Results.pos_mark})', f'p({Results.pos_mark})')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    tmp = extract_and_merge(Results.apoe_pathology_interaction, f'β({Results.apoe_identifier}*AD)', 
                            f'p({Results.apoe_identifier}*AD)')
    models_df = pd.merge(models_df, tmp, on=[Results.summarize_oncol], how='left')

    models_df['abs_beta1'] = models_df['β2'].abs()
    models_df = models_df.sort_values(by='abs_beta1', ascending=False)

    models_df.set_index(Results.summarize_oncol, inplace=True)

    models = models_df[['β2', 'β4', 'β6', 'β7', f'β({Results.neg_mark})', f'β({Results.apoe_identifier}*age in {Results.neg_mark})', f'β({Results.pos_mark})', f'β({Results.apoe_identifier}*AD)']]
    annotations = models_df[['p2', 'p4', 'p6', 'p7', f'p({Results.neg_mark})', f'p({Results.apoe_identifier}*age in {Results.neg_mark})', f'p({Results.pos_mark})', f'p({Results.apoe_identifier}*AD)']]

    models.columns = [f'{Results.apoe_identifier}_{x}' for x in models.columns]
    annotations.columns = [f'{Results.apoe_identifier}_{x}' for x in annotations.columns]

    return models, annotations


def prepare_heatmap_inputs(
    pc1: list[str],
    Results1,
    Results2,
    target_sufix1,
    target_sufix2,
    beta_threshold: float = 0.1
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function extracts results from two ResultDataLoader instances, merges them,
    and prepare the data for heatmap plot.
    """
    models1, annotations1 = extract_mainmodel_results(pc1, Results1)
    models2, annotations2 = extract_mainmodel_results(pc1, Results2)

    models = pd.merge(models1, models2, left_index=True, right_index=True, how='outer')
    annotations = pd.merge(annotations1, annotations2, left_index=True, right_index=True, how='outer')
    
    β2_col1 = f"{Results1.apoe_identifier}_β{target_sufix1}"
    β2_col2 = f"{Results2.apoe_identifier}_β{target_sufix2}"
    models['bold'] = models[β2_col1] * models[β2_col2]
    models['bold'] = models['bold'].apply(lambda x: 'opposite' if x < 0 else 'same' if x > 0 else None)

    models['type'] = models.index.to_series().apply(
        lambda x: f"{Results1.apoe_identifier.lower()} mediator" if x in Results1.path1_proteins
        else (f"{Results2.apoe_identifier.lower()} mediator" if x in Results2.path1_proteins else '')
    )

    ref2bold = dict(zip(models.index.tolist(), models['bold'].tolist()))
    ref2type = dict(zip(models.index.tolist(), models['type'].tolist()))

    models['mean'] = (models[β2_col1].abs() + models[β2_col2].abs()) / 2
    models = models[models['mean'] > beta_threshold]
    models['mean'] = 0 - models['mean']
    models.sort_values(['bold', 'mean'], ascending=True, inplace=True)

    beta_cols = [col for col in models.columns if col.startswith(f"{Results1.apoe_identifier}_β") or col.startswith(f"{Results2.apoe_identifier}_β")]
    pval_cols = [col for col in annotations.columns if col.startswith(f"{Results1.apoe_identifier}_p") or col.startswith(f"{Results2.apoe_identifier}_p")]

    annotations = annotations.loc[models.index, pval_cols]
    models = models[beta_cols]

    assert all(models.index == annotations.index), "Row index mismatch between models and annotations"
    assert models.shape[1] == annotations.shape[1], "Number of β and p-value columns must match"
    def suffix(col): return col.split('_')[-1].replace('β', '').replace('p', '')
    assert all(
        suffix(m) == suffix(p) for m, p in zip(models.columns, annotations.columns)
    ), "Column suffixes (e.g., β2 vs p2) do not align properly"

    return models, annotations, ref2bold, ref2type

def return_df4onesymbol(symbol, ResultClass) -> pd.DataFrame:
    """ Extracts results for a specific protein symbol from all models in a analysis results.
    Args:
        symbol (str): The protein symbol to filter results by.
        ResultClass (ResultDataLoader): An instance of ResultDataLoader containing the analysis results.
    Returns:
        pd.DataFrame: A DataFrame containing results for the specified protein symbol across various models.
    """
    apoe = ResultClass.apoe_identifier
    dataset = ResultClass.dataset

    apoe2protein         = {'df': pd.DataFrame() if ResultClass.apoe_result.data is None else ResultClass.apoe_result.data,
                            'model': f'APOE[{apoe}]'} 
    ad2protein           = {'df': pd.DataFrame() if ResultClass.pathology_result.data is None else ResultClass.pathology_result.data,
                            'model': f'ADorAβ[{apoe}]'}
    apoe2protein_adjad   = {'df': pd.DataFrame() if ResultClass.apoe_adj_result.data is None else ResultClass.apoe_adj_result.data,
                            'model': f'APOE_adjADorAβ[{apoe}]'}
    ad2protein_adjapoe   = {'df': pd.DataFrame() if ResultClass.pathology_adj_result.data is None else ResultClass.pathology_adj_result.data,
                            'model': f'ADorAβ_adjAPOE[{apoe}]'}
    # cu_allage            = {'df': pd.DataFrame() if ResultClass.neg_allage.data is None else ResultClass.neg_allage.data,
    #                         'model': f'APOE in CUorAβ-[{apoe}]'}
    # apoe_age_interaction = {'df': pd.DataFrame() if ResultClass.apoe_age_interaction.data is None else ResultClass.apoe_age_interaction.data,
    #                         'model': f'APOE*age in CUorAβ-[{apoe}]'}
    # protein_mediation    = {'df': pd.DataFrame() if ResultClass.protein_mediation.data is None else ResultClass.protein_mediation.data,
    #                         'model': f'APOE=>protein=>ADorAβ[{apoe}]'}
    # pathology_mediation  = {'df': pd.DataFrame() if ResultClass.pathology_mediation.data is None else ResultClass.pathology_mediation.data,
    #                         'model': f'ADorAβ=>protein=>APOE[{apoe}]'}
    
    models_df = pd.DataFrame()
    for data in [apoe2protein, ad2protein, apoe2protein_adjad, ad2protein_adjapoe, 
               #cu_allage, apoe_age_interaction, 
               #protein_mediation, pathology_mediation
               ]:
        df = data['df']
        if df.empty or len(df) == 0 or 'symbol' not in df.columns:
            continue
        if symbol in df['symbol'].tolist():
            df = df[df['symbol'] == symbol]
            cols = [ResultClass.summarize_oncol, 'significance','p', 'p_adjusted','symbol']
            if 'standardized_beta' in data['df'].columns.tolist():
                cols.append('standardized_beta')
            if 'prop_m' in data['df'].columns.tolist():
                cols.append('prop_m')

            df = df[cols]
            if 'prop_m' in data['df'].columns.tolist():
                df.rename(columns={'prop_m':'standardized_beta'}, inplace=True)

            df['model'] = data['model']

            models_df = pd.concat([models_df, df], join='outer')
            models_df['id'] = models_df[ResultClass.summarize_oncol]
            models_df['dataset'] = dataset

    return models_df

def mediation_heatmap_data(Results, proteins:list, top_n=20, abs_cutoff=0.1):
    oncol = Results.summarize_oncol

    tmp1 = Results.protein_mediation.data
    tmp2 = Results.pathology_mediation.data

    def modify_df(df, oncol, sub1, sub2):
        df = df[df[oncol].isin(proteins)][[oncol, 'prop_m', 'significance']]
        df.rename(columns={'prop_m': sub1, 'significance': sub2}, inplace=True)
        if 'label' not in df.columns:
            df['label'] = df[oncol].apply(lambda x: Results.id2label[x])

        return df

    mediation_df = pd.DataFrame({oncol: list(proteins)})
    mediation_df['label'] = mediation_df[oncol].apply(lambda x: Results.id2label[x])

    if tmp1 is not None:
        tmp1 = modify_df(tmp1, oncol, 'Path1', 'p1')
        mediation_df = pd.merge(mediation_df, tmp1, on=[oncol, 'label'], how='outer')

    if tmp2 is not None:
        tmp2 = modify_df(tmp2, oncol, 'Path2', 'p2')
        mediation_df = pd.merge(mediation_df, tmp2, on=[oncol, 'label'], how='outer')

    def get_type(x):
        if x in Results.path1_proteins:
            return 'protein total' if x in Results.path1_full else 'protein partial'
        elif x in Results.path2_proteins:
            return 'ad total' if x in Results.path2_full else 'ad partial'
        else:
            return 'none'

    mediation_df['type'] = mediation_df[oncol].apply(get_type)

    path1 = mediation_df[mediation_df[oncol].isin(Results.path1_proteins)]
    path1['Path1_abs'] = path1['Path1'].abs()
    path1 = path1[path1['Path1_abs'] > abs_cutoff]
    path1 = path1.nlargest(top_n, 'Path1_abs')

    path2 = mediation_df[mediation_df[oncol].isin(Results.path2_proteins)]
    path2['Path2_abs'] = path2['Path2'].abs()
    path2 = path2[path2['Path2_abs'] > abs_cutoff]
    path2 = path2.nlargest(top_n, 'Path2_abs')

    concated = pd.concat([path1, path2], join='inner')
    concated.set_index(oncol, inplace=True)

    mediations = concated[['Path1', 'Path2']]
    mediations.columns = [
        f'{Results.apoe_identifier} → protein → {Results.pathology_identifier}',
        f'{Results.apoe_identifier} → {Results.pathology_identifier} → protein'
    ]
    mediations = mediations * 100

    annotations = concated[['p1', 'p2']]

    return mediations, annotations, dict(zip(mediation_df[oncol].tolist(), mediation_df['type'].tolist()))


def sup_of_onedf(df, dataset, df_title):
    if df is None or df.empty:
        return None, None
    
    dataset = dataset.split('-')[-1]

    if 'uniprot' in df.columns.tolist():
        df = df.rename(columns = {'uniprot':'UniProt'})
    
    if 'uniprot.x' in df.columns.tolist():
        df = df.rename(columns = {'uniprot.x':'UniProt'})

    if dataset == 'SomaLogic':
        head1_part1 = ['Somamers, corresponding genes and labels used for plots and annotation', '', '', '']
        head2_part1 = ['apt_name', 'UniProt', 'EntrezGeneSymbol', 'label']
        columns_part1 = ['apt_name', 'UniProt', 'symbol', 'label']
    elif dataset == 'OLINK':
        head1_part1 = ['OLINK proteins, corresponding genes and labels used for plots and annotation', '', '', '']
        head2_part1 = ['Protein_id', 'UniProt', 'EntrezGeneSymbol', 'label']
        columns_part1 = ['Protein_id', 'UniProt', 'symbol', 'label']
    elif dataset == 'RNAseq':
        head1_part1 = ['Gene expression from bulk RNA-seq', '', '']
        head2_part1 = ['Gene expression', 'EntrezGeneSymbol', 'label']
        columns_part1 = ['Protein_id', 'symbol', 'label']
    elif dataset == 'TMT':
        head1_part1 = ['Proteins measured by TMT-MS', '', '', '']
        head2_part1 = ['Protein_id', 'UniProt', 'EntrezGeneSymbol', 'label']
        columns_part1 = ['Protein_id', 'UniProt', 'symbol', 'label']
    else:
        print(f"[WARN] Cannot output dataset other than SomaLogic, OLINK, TMT, RNAseq")
        return None, None

    if 'prop_m' in df.columns.tolist():
        columns_part2 = [
            "acme", "acme_p", "p_adjusted", 
            "ade", "ade_p", "ade_p_adjusted", 
            "prop_m", "prop_m_p", "prop_m_p_adjusted"
        ]
        head1_part2 = [df_title] + [''] * (len(columns_part2) - 1)
        head2_part2 = columns_part2
    else:
        columns_part2 = ['coef','standardized_beta', 'se', 'resid', 'p', 'p_adjusted']
        head1_part2 = [df_title] + [''] * (len(columns_part2) - 1)
        head2_part2 = ['coef', 'std.beta','se', 'resid', 'p', 'p_adjusted']

    try:
        tmp = df[columns_part1 + columns_part2].copy()
    except KeyError as e:
        print(df.columns.tolist())
        print(f"[ERROR] Missing expected columns: {e}")
        return None, None

    head1 = head1_part1 + head1_part2
    head2 = head2_part1 + head2_part2

    header_df = pd.DataFrame([head1, head2], columns=columns_part1 + columns_part2)
    tmp = pd.concat([header_df, tmp], ignore_index=True)

    return tmp, columns_part1

def pick_and_rename_multi(df, key_col, cols_map, prefix):
    """
    cols_map: dict like {"beta": "std_beta", "p": "pvalue", "fdr": "p_adjusted"}
    prefix: dataset name, will create columns like f"{prefix}__beta"
    """
    if df is None or df.empty:
        return None

    available = {new: old for new, old in cols_map.items() if old in df.columns}
    if key_col not in df.columns or len(available) == 0:
        return None

    out = df[[key_col] + list(available.values())].copy()
    out.rename(columns={old: f"{prefix}__{new}" for new, old in available.items()}, inplace=True)
    return out

def best_per_symbol(df, symbol_col="symbol", p_col="p_adjusted"):
    """For duplicated proteins (gene symbol), select the row with smallest p_col."""
    if df is None or df.empty:
        return df
    if symbol_col not in df.columns:
        raise ValueError(f"Missing symbol_col={symbol_col}")
    if p_col not in df.columns:
        raise ValueError(f"Missing p_col={p_col} in df.columns: {list(df.columns)[:20]} ...")

    tmp = df.copy()
    tmp[p_col] = pd.to_numeric(tmp[p_col], errors="coerce")

    g = (tmp.dropna(subset=[p_col])
            .sort_values(p_col)
            .groupby(symbol_col, as_index=False)
            .head(1))
    if g.empty:
        g = (tmp.groupby(symbol_col, as_index=False).head(1))
    return g.reset_index(drop=True)