import pandas as pd

# Caminho para o arquivo de anotação e o arquivo de metilação
annotation_file = "RANCOL_sorted_file.bed"
methylation_file = "methylation.txt"

# Carregar o arquivo de anotação com nomes de colunas únicos
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None,
                             names=["chr", "start", "end", "dot1", "dot2", "strand", "source", "feature", "dot3", "dot4", "ID"])

# Carregar o arquivo de metilação
methylation_df = pd.read_csv(methylation_file, sep='\t', header=None, names=["chr", "start"])

# Converter a coluna "chr" para string para garantir que o método str.replace funcione
methylation_df["chr"] = methylation_df["chr"].astype(str)
annotation_df["chr"] = annotation_df["chr"].astype(str)

# Remover o prefixo "chr" dos cromossomos, se houver
methylation_df["chr"] = methylation_df["chr"].str.replace("chr", "", regex=False)
annotation_df["chr"] = annotation_df["chr"].str.replace("chr", "", regex=False)

# Exibir as primeiras linhas para verificar o conteúdo
print("Verificando as primeiras linhas do arquivo de metilação:")
print(methylation_df.head())

# Limpar espaços e garantir que a coluna "start" do arquivo de metilação seja numérica
methylation_df["start"] = pd.to_numeric(methylation_df["start"], errors="coerce")

# Remover linhas com valores NaN na coluna "start" (se existirem)
methylation_df = methylation_df.dropna(subset=["start"])

# Verificar as linhas restantes após limpeza
print("Após a limpeza, as primeiras linhas do arquivo de metilação:")
print(methylation_df.head())

# Garantir que as colunas "start" de ambos os DataFrames estejam no mesmo tipo (inteiro)
annotation_df["start"] = annotation_df["start"].astype(int)

# Realizar o merge dos dados de metilação com a anotação com base no cromossomo e posição
merged_df = pd.merge(methylation_df, annotation_df, how="left", left_on=["chr", "start"], right_on=["chr", "start"])

# Verificar o merge para garantir que as correspondências estejam sendo feitas
print("DataFrame após merge:")
print(merged_df.head())

# Marcar as linhas sem anotação
merged_df["feature"].fillna("not found", inplace=True)

# Selecionar e salvar as colunas relevantes: cromossomo, posição e anotação da feature
output_df = merged_df[["chr", "start", "feature"]]

# Salvar o arquivo de saída com a anotação
output_file = "methylation_anotated.txt"
output_df.to_csv(output_file, sep='\t', index=False, header=False)

print(f'Anotação concluída. O arquivo foi salvo em: {output_file}')
