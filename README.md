# MirrorMaker

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14603319.svg)](https://doi.org/10.5281/zenodo.14603319)

🧬 A Bacterial rRNA Operon Extraction and Curation Pipeline  
For microbial genomics researchers extracting 16S–23S rRNA operons from genome annotation files.  
GFF/FNA-based operon detection → GTDB taxonomy annotation → Redundancy & ambiguity filtering

---

### 📁 Directory Structure

```
MirrorMaker/
├── input/
│   ├── fna_sample100/             # FNA genome files (*.fna)
│   ├── gff_sample100/             # GFF annotation files (*.gff)
│   └── gtdb_taxonomy_r220.tsv     # GTDB taxonomy mapping file
├── output/
│   ├── 1.Operon_Filtered.txt          # Operons filtered by length
│   ├── 2.Operon_Extraction_Summary.txt # GFF processing summary
│   ├── 3.Operon_Extraction_Results.txt # Extracted operon sequences
│   ├── 4.Operon_Curation_Detail.txt   # With taxonomy & duplication flags
│   └── 5.Operon_Curation_Results.txt  # Final curated output
└── main.py           # Main pipeline script
```

---

### ✅ Requirements

- Python 3.7+
- Required modules:
  - `pandas`
  - `glob`
  - `datetime`
  - `re`

```bash
pip install pandas
```

---

### 🚀 How to Run

```bash
python main.py
```

---

### 🔧 Pipeline Overview

#### 1. GFF Parsing & Operon Detection
- Detects rRNA operons composed of 16S and 23S features (strand-aware)
- Accepts only operons within 3,500–7,000 bp range

#### 2. Sequence Extraction
- Retrieves sequence segments from corresponding `.fna` files
- Performs reverse-complement if strand is negative (`-`)

#### 3. GTDB Taxonomy Annotation
- Uses `gtdb_taxonomy_r220.tsv` to annotate taxonomic information
- Fields extracted: Phylum, Class, Order, Family, Genus, Species

#### 4. Redundancy Filtering
- Identifies and flags/removes:
  - **Intra-genome duplicated** sequences
  - **Intra-species duplicated** sequences
  - **Inter-species duplicated** sequences

#### 5. Ambiguous Nucleotide Filtering
- Filters out sequences with more than 5 ambiguous bases (`N`)

#### 6. Final Sequence Curation
- Drops entries missing taxonomy info
- Reassigns operon CopyNo after removals
- Saves final curated list

---

### 📄 Input File Formats

#### GFF File

- Must follow standard 9-column GFF format
- Should include `rRNA` features with `16S` and `23S` in attributes

#### FNA File

- FASTA format genome file
- Header lines must include accession (e.g. `>accession...`)

#### GTDB Taxonomy File

- Tab-separated format:

```
Genbank_ID <TAB> taxonomy_string
```

Example:

```
GCA_000599625.1  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
```

---

### 📤 Output Files

| File                                | Description                                   |
|-------------------------------------|-----------------------------------------------|
| `1.Operon_Filtered.txt`             | Length-filtered (discarded) operons           |
| `2.Operon_Extraction_Summary.txt`   | Summary of genome and GFF file status         |
| `3.Operon_Extraction_Results.txt`   | Extracted operon data with sequences          |
| `4.Operon_Curation_Detail.txt`      | Operon data + taxonomy + duplication flags    |
| `5.Operon_Curation_Results.txt`     | Final clean set for downstream analysis       |


##### Example: `1.Operon_Filtered.txt`

```tsv
GenbankNo	Accession	CopyNo	Start	End	Strand	Length
GCA_000599625.1	CP007390.1	operon_1	4715033	4723387	+	8355
GCA_027286205.1	CP114738.1	operon_1	462430	472704	+	10275
GCA_000185245.1	CP002336.1	operon_1	1014367	1145719	-	131353
```
- Note: These operons were excluded because their lengths were either too short or too long (outside the 3,500–7,000 bp threshold).*

##### Example: `2.Operon_Extraction_Summary.txt`
```txt
### Operon Extraction Summary ###
Total genomes processed: 1
Genomes with no operons: 0
Genomes with empty GFF files: 0
```

---
### 📚 Citation

If you use this pipeline in your research or publication, please cite:

> *Zenodo*. [https://doi.org/10.5281/zenodo.14603319](https://doi.org/10.5281/zenodo.14603319)

---

### 📬 Contact

Contact: **Jisol Lee** – [leejisol@snu.ac.kr](mailto:leejisol@snu.ac.kr)
