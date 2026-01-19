import pandas as pd
import pyranges as pr

# === 1. Ler os arquivos ===
te_df = pd.read_csv("Psta-TEs-TEanno.bed", sep="\t", header=None,
                    names=["chr", "start", "end", "TE_type", "score", "strand"])

# Ler mantendo a primeira coluna como índice
meth_df = pd.read_csv("allo_CpG_DMR_significativas.csv", index_col=0)

# === 2. Adicionar ID único baseado no índice original ===
meth_df['original_id'] = meth_df.index

# === 3. Converter para PyRanges ===
te_gr = pr.PyRanges(te_df.rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"}))
meth_gr = pr.PyRanges(meth_df.rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"}))

# === 4. Fazer interseção ===
overlap = meth_gr.join(te_gr, how="left")
overlap_df = overlap.as_df()

# === 5. Consolidar TEs por DMR ===
te_annotations = overlap_df.groupby('original_id')['TE_type'].apply(
    lambda x: ';'.join(sorted(x.dropna().unique())) if x.notna().any() else 'no_TE'
).reset_index(name='TE_annot')

te_counts = overlap_df.groupby('original_id')['TE_type'].apply(
    lambda x: len(x.dropna().unique())
).reset_index(name='TE_count')

# Combinar anotações
te_info = pd.merge(te_annotations, te_counts, on='original_id')

# === 6. Juntar com dados originais ===
final_df = meth_df.merge(te_info, on='original_id', how='left')

# Preencher valores faltantes
final_df['TE_annot'] = final_df['TE_annot'].fillna('no_TE')
final_df['TE_count'] = final_df['TE_count'].fillna(0)

# Remover coluna temporária
final_df = final_df.drop('original_id', axis=1)

# === 7. Salvar resultado ===
final_df.to_csv("allo_CpG_DMR_significativas_annot_consolidated.csv")

print("✅ Arquivo salvo com sucesso!")
print(f"Total de DMRs: {len(final_df)}")
print(f"DMRs com TEs: {len(final_df[final_df['TE_count'] > 0])}")
