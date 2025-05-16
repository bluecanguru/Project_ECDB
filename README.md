# Projeto_ECDB
## 🧬Introdução

### Tipo de Cancro  
Esse estudo faz parte do The Cancer Genome Atlas (TCGA), um grande projeto colaborativo que visa caracterizar geneticamente diversos tipos de cancro usando tecnologias de alto rendimento, como RNA-Seq para expressão genética.  
Este dataset está focado no *Glioma de baixo grau* (Low Grade Glioma - LGG), um tipo de tumor cerebral com evolução mais lenta do que o glioblastoma, mas que pode ser fatal em vários casos. 
Apesar da sua progressão mais lenta, o LGG apresenta uma elevada heterogeneidade molecular, tornando-se relevante para estudos de estratificação de pacientes e identificação de subtipos tumorais.  

### Relevância científica  
Este tipo de dados permite identificar *genes diferencialmente expressos* entre:  
- Subtipos tumorais
- Grupos com ou sem mutações específicas
- Amostras normais e tumorais

No caso do estudo, estão disponíveis apenas as 514 amostras tumorais de glioma de baixo grau (LGG).  

Essas análises podem contribuir para:  
- Descoberta de biomarcadores moleculares  
- Estratificação de pacientes com base em perfis de expressão  
- Identificação de possíveis alvos terapêuticos  

Além disso, os dados podem ser utilizados em análises de enriquecimento funcional (como GSEA ou over-representation analysis) para responder a perguntas como:  

> Quais as vias biológicas que estão mais ativas ou reprimidas em determinados grupos de pacientes?  

Também é possível associar perfis de expressão a desfechos clínicos, como tempo de sobrevida ou resposta ao tratamento, contribuindo para avanços na medicina personalizada.  

O dataset é compatível com técnicas de machine learning, permitindo:  
- Classificação de pacientes  
- Seleção de features relevantes  
- Predição de prognóstico 

O acesso aberto e a documentação clara disponíveis pelo cBioPortal reforçam a reprodutibilidade científica e a integração com ferramentas computacionais, como APIs e bibliotecas em R ou Python.


## ⚡Dados   
No *cBioPortal*, este estudo inclui:  
- Dados clínicos (idade, sexo, sobrevida, etc.)  
- Dados de Expressão Génica (RNA-Seq)  
- Dados de Mutação Somática  
- Dados de CNV (alterações no número de cópias)  
- Dados de Metilação 
- Dados de Fusões 

Todos esses dados são integráveis, permitindo análises multi-ómicas para uma visão mais ampla dos mecanismos moleculares envolvidos no desenvolvimento e progressão do glioma.  

### Dados de Expressão Génica (*RNA-Seq*)  
Os dados de expressão génica foram processados usando o método RSEM, com normalização por lote (batch-normalized) a partir da plataforma *Illumina HiSeq RNASeqV2*.  

A matriz de expressão contém *milhares de genes* (geralmente cerca de 20.000), permitindo análises em larga escala de:  
- Expressão diferencial  
- Coexpressão  
- Redes regulatórias 

## 🔍Fonte dos Dados
Os dados provêm do estudo *"LGG TCGA PanCancer Atlas 2018"*, disponível no portal cBioPortal: [🔗 cBioPortal](https://www.cbioportal.org/study/summary?id=lgg_tcga_pan_can_atlas_2018).  

## 📖Estrutura Repositório
    📂 data/                # Raw and processed data files
    📂 scripts/             # R scripts for data analysis
    📂 docs/                # Documentation and reports (HTML, R Markdown)
    📄 README.md            # Project overview and instructions

## 🛠 Setup
1. Clone the repository:
```bash
git clone https://github.com/bluecanguru/Project_ECDB
```
2. Install required dependencies:
```bash
pip install -r requirements.txt
```

## ⚙️Dependências
```{r}
install.packages(c(
  "ggplot2",      # Visualization
  "gplots",       # Enhanced plots
  "gridExtra",    # Combine multiple plots
  "DT",           # Interactive datatables
  "dplyr",        # Data manipulation
  "plyr",         # Data manipulation (older)
  "stringr",      # String operations
  "knitr",        # RMarkdown rendering
  "data.table",   # Efficient data handling
  "viridis",      # Color palettes
  "openxlsx",     # Excel export
  "corrplot",     # Correlation plots
  "mgcv",         # GAM models
  "ggpubr",       # Publication-ready plots
  "GGally",       # Pair plots (extension of ggplot2)
  "htmltools",    # To print text in specific format
  "gprofiler2"    # g:Profiler interface (CRAN but bio-focused)
  "caret",        # Modelation
  "randomForest", # Random Forest
  "smotefamily",  # SMOTE
  "xgboost",      # xgboost
  "e1071",        # SVM
  "pROC",         # ROC CURVE
  "MLmetrics",    # F1_score
  "pander",       # Tables
  "tibble",       # Interactive datatables
  "Rtsne",        # UMAP
  "uwot",         # T-SNE
  "factoextra"    # T-SNE
))
```

```{r}
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR",             # RNA-seq analysis
  "limma",             # Linear models for microarray/RNA-seq
  "clusterProfiler",   # Functional enrichment
  "org.Hs.eg.db",      # Human gene annotation
  "AnnotationDbi",     # Annotation infrastructure
  "biomaRt",           # Interface to BioMart databases
  "enrichplot"         # Visualization of enrichment results
))
```

## 📝Contribuição
- [Cátia Rosário](https://github.com/bluecanguru)
- [Vanessa Rodriguez](https://github.com/VaneBR)
- [André Dias](https://github.com/Diasf333)

## 📜Licença
This project is open source and available under the MIT License.