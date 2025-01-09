import os
import glob
import datetime
import pandas as pd
import traceback
import re

####################################################################################
# Utility Functions
####################################################################################

def WriteLog(functionname, msg):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}][{functionname}] {msg}")

def create_folder(folder_path):
    try:
        os.makedirs(folder_path, exist_ok=True) 
        print(f"\nDirectory '{folder_path}' is ready.\n")
    except Exception as e:
        print(f"Error creating directory '{folder_path}': {e}")

####################################################################################
# Main Class
####################################################################################

class MIrRORv2_Bacteria:
    def __init__(self):
        self.setup_paths()

    def setup_paths(self):
        base_path = os.path.dirname(os.path.abspath(__file__))
        self.path_input_fna = f"{base_path}/input/fna_sample100/"
        self.path_input_gff = f"{base_path}/input/gff_sample100/GCA*.gff"
        self.path_input_gtdb_taxon = f"{base_path}/input/gtdb_taxonomy_r220.tsv"
        self.path_output_1_filtered = f"{base_path}/output/1.Operon_Filtered.txt"
        self.path_output_2_summary = f"{base_path}/output/2.Operon_Extraction_Summary.txt"
        self.path_output_3_extraction = f"{base_path}/output/3.Operon_Extraction_Results.txt"
        self.path_output_4_curation = f"{base_path}/output/4.Operon_Curation_Detail.txt"
        self.path_output_5_final = f"{base_path}/output/5.Operon_Curation_Results.txt"

    def read_gff(self):
        try:
            gff_files = glob.glob(self.path_input_gff)
            data, data_filter = [], []
            empty_files, no_operons = [], []

            for gff_file in gff_files:
                GCA_name = "GCA_" + os.path.basename(gff_file).split("_")[1]
                if os.path.getsize(gff_file) == 0:
                    empty_files.append(GCA_name)
                    continue

                temp_data = self.process_gff_file(gff_file, GCA_name, data_filter)
                if temp_data:
                    data.extend(temp_data)
                else:
                    no_operons.append(GCA_name)

            self.save_operon_data(data, data_filter, empty_files, no_operons)
            
        except Exception as e:
            traceback.print_exc()
            WriteLog("Read gff files", f"Error: {e}")

    def process_gff_file(self, gff_file, GCA_name, data_filter):
        temp_data = []
        operon_number, operon_num_fil = 1, 1
    
        with open(gff_file, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines) - 1):
                line = lines[i].strip()
                next_line = lines[i + 1].strip()
    
                if not self.is_valid_line_pair(line, next_line):
                    continue
    
                columns, next_columns = line.split('\t'), next_line.split('\t')
                accession, feature_type, start, end, strand, attributes = (
                    columns[0], columns[2], int(columns[3]), int(columns[4]), columns[6], columns[8]
                )
                next_accession, next_feature_type, next_start, next_end, next_strand, next_attributes = (
                    next_columns[0], next_columns[2], int(next_columns[3]), int(next_columns[4]), next_columns[6], next_columns[8]
                )
    
                if strand == '+' and feature_type == 'rRNA' and next_feature_type == 'rRNA' and '16S' in attributes and '23S' in next_attributes:
                    if accession == next_accession and strand == next_strand:
                        operon_start, operon_end = start, next_end
                        length = operon_end - operon_start + 1
                        if 3500 <= length <= 7000:
                            temp_data.append([GCA_name, accession, f"operon_{operon_number}", operon_start, operon_end, strand, length])
                            operon_number += 1
                        else:
                            data_filter.append([GCA_name, accession, f"operon_{operon_num_fil}", operon_start, operon_end, strand, length])
                            operon_num_fil += 1
    
                elif strand == '-' and feature_type == 'rRNA' and next_feature_type == 'rRNA' and '23S' in attributes and '16S' in next_attributes:
                    if accession == next_accession and strand == next_strand:
                        operon_start, operon_end = start, next_end
                        length = operon_end - operon_start + 1
                        if 3500 <= length <= 7000:
                            temp_data.append([GCA_name, accession, f"operon_{operon_number}", operon_start, operon_end, strand, length])
                            operon_number += 1
                        else:
                            data_filter.append([GCA_name, accession, f"operon_{operon_num_fil}", operon_start, operon_end, strand, length])
                            operon_num_fil += 1
    
        return temp_data

    def is_valid_line_pair(self, line, next_line):
        return all([
            not line.startswith("#"),
            len(line.split('\t')) >= 9,
            not next_line.startswith("#"),
            len(next_line.split('\t')) >= 9
        ])

    def save_operon_data(self, data, data_filter, empty_files, no_operons):
        if data:
            df_operon = pd.DataFrame(data, columns=["GenbankNo", "Accession", "CopyNo", "Start", "End", "Strand", "Length"])
            df_operon.to_csv(self.path_output_3_extraction, sep='\t', index=False)
    
        if data_filter:
            df_filter = pd.DataFrame(data_filter, columns=["GenbankNo", "Accession", "CopyNo", "Start", "End", "Strand", "Length"])
            df_filter.to_csv(self.path_output_1_filtered, sep='\t', index=False)
    
        with open(self.path_output_2_summary, 'w') as f:
            f.write(f"### Operon Extraction Summary ###\n")
            f.write(f"Total genomes processed: {len(data) + len(data_filter)}\n")
            f.write(f"Genomes with no operons: {len(no_operons)}\n")
            f.write(f"Genomes with empty GFF files: {len(empty_files)}\n")

    def extract_sequences(self):
        try:
            df_operon = pd.read_csv(self.path_output_3_extraction, sep='\t')
            df_operon['Sequence'] = df_operon.apply(lambda row: self.get_sequence(row), axis=1)
            df_operon.to_csv(self.path_output_3_extraction, sep='\t', index=False)
        except Exception as e:
            traceback.print_exc()
            WriteLog("extract_sequences", f"Error: {e}")
    
    def get_sequence(self, row):
        fna_files = glob.glob(f"{self.path_input_fna}/{row['GenbankNo']}*.fna")
        if not fna_files:
            return 'N/A'
        with open(fna_files[0], 'r') as fna_file:
            sequence = self.extract_sequence_from_fna(fna_file, row['Accession'])
            return sequence[row['Start'] - 1: row['End']] if sequence else 'N/A'
    
    def extract_sequence_from_fna(self, fna_file, accession):
        capture = False
        sequence_parts = []
        for line in fna_file:
            if line.startswith(">"):
                if capture:
                    break
                if accession in line:
                    capture = True
            elif capture:
                sequence_parts.append(line.strip())
        return "".join(sequence_parts)

    def reverse_and_complement(self):
        try:
            df = pd.read_csv(self.path_output_3_extraction, sep='\t')
            df['Sequence'] = df.apply(lambda row: self.reverse_complement(row['Sequence']) if row['Strand'] == '-' else row['Sequence'], axis=1)
            df.to_csv(self.path_output_3_extraction, sep='\t', index=False)
        except Exception as e:
            traceback.print_exc()
            WriteLog("reverse_and_complement", f"Error: {e}")

    def reverse_complement(self, sequence):
        return sequence[::-1].translate(str.maketrans('ATCGatcg', 'TAGCtagc'))

    def append_gtdb_taxonomy(self):
        taxa_info = {}
    
        with open(self.path_input_gtdb_taxon, 'r', encoding='utf-8-sig') as f:  # encoding 수정
            next(f)
            for line in f:
                col = line.strip().split('\t')
                GCA_name = col[0]
                gtdb_taxa = col[1]
                taxa_list = gtdb_taxa.split(";")
                taxa_info[GCA_name] = {
                    'Taxonomy': gtdb_taxa,
                    'Phylum': taxa_list[1].split('__')[1] if len(taxa_list) > 1 else "",
                    'Class': taxa_list[2].split('__')[1] if len(taxa_list) > 2 else "",
                    'Order': taxa_list[3].split('__')[1] if len(taxa_list) > 3 else "",
                    'Family': taxa_list[4].split('__')[1] if len(taxa_list) > 4 else "",
                    'Genus': taxa_list[5].split('__')[1] if len(taxa_list) > 5 else "",
                    'Species_with_suffix': taxa_list[6].split('__')[1] if len(taxa_list) > 6 else ""
                }
    
        df = pd.read_csv(self.path_output_3_extraction, sep='\t', dtype=str)
        for key in ['Taxonomy', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species_with_suffix']:
            df[key] = None
    
        for idx, row in df.iterrows():
            GCA_name = row['GenbankNo']
            if GCA_name in taxa_info:
                for key, value in taxa_info[GCA_name].items():
                    df.at[idx, key] = value
    
        df.to_csv(self.path_output_4_curation, sep='\t', index=False)
        WriteLog("Step_01", "GTDB Taxonomy appended successfully.")

    def curation_redundant_sequence(self):
        df = pd.read_csv(self.path_output_4_curation, sep='\t')
        df['IntraGenomeDup'] = None
        df['IntraSpeciesDup'] = None
        df['InterSpeciesDup'] = None

        def rm_species_segmentation(Species_with_suffix):
            if not isinstance(Species_with_suffix, str):
                return ''
            return re.sub(r'_[A-Z]$', '', re.sub(r'_[A-Z](?=\s|$)', '', Species_with_suffix))

        df.insert(df.columns.get_loc('Species_with_suffix') + 1, 'Species', df['Species_with_suffix'].apply(rm_species_segmentation))

        for genbank_no, group in df.groupby('GenbankNo'):
            duplicates = group[group.duplicated('Sequence', keep=False)]
            for seq, sub_group in duplicates.groupby('Sequence'):
                sorted_idx = sub_group['CopyNo'].str.extract('(\d+)', expand=False).astype(int).sort_values().index
                df.loc[sorted_idx[1:], 'IntraGenomeDup'] = 'Rm'

        for taxonomy, group in df.groupby('Species'):
            duplicates = group[group.duplicated('Sequence', keep=False)]
            for idx in duplicates.index:
                df.at[idx, 'IntraSpeciesDup'] = 'Rm'

        for sequence, group in df.groupby('Sequence'):
            if len(group['Species'].unique()) > 1:
                for idx in group.index:
                    df.at[idx, 'InterSpeciesDup'] = 'Rm'

        df.to_csv(self.path_output_4_curation, sep='\t', index=False)
        WriteLog("Step_02", "Redundant sequence curation completed successfully.")

    def curation_ambiguous_nucleotides(self):
        df = pd.read_csv(self.path_output_4_curation, sep='\t')
        df['AmbigNT'] = None
        df['AmbigNT>5'] = None

        df['Sequence'] = df['Sequence'].fillna('')
        df['AmbigNT'] = df['Sequence'].apply(lambda x: x.upper().count('N'))
        df['AmbigNT>5'] = df['AmbigNT'].apply(lambda x: 'Rm' if x > 5 else '')
        df.to_csv(self.path_output_4_curation, sep='\t', index=False)
        WriteLog("Step_03", "Ambiguous nucleotide curation completed successfully.")
        
    def curation_sequence_removal(self):
        df = pd.read_csv(self.path_output_4_curation, sep='\t')
        removal_columns = ['IntraGenomeDup', 'InterSpeciesDup', 'AmbigNT>5'] #'IntraSpeciesDup'
        df[removal_columns] = df[removal_columns].fillna('').astype(str)
        df = df[~df[removal_columns].apply(lambda x: x.str.contains('Rm', na=False)).any(axis=1)]
        
        taxonomy_columns = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        df = df.dropna(subset=taxonomy_columns, how='any')
        
        columns_to_drop = ['IntraGenomeDup', 'IntraSpeciesDup', 'InterSpeciesDup', 'AmbigNT', 'AmbigNT>5', 'Species_with_suffix']
        df = df.drop(columns=columns_to_drop, errors='ignore')

        df['CopyNo'] = df.groupby('GenbankNo').cumcount().apply(lambda x: f"operon_{x + 1}")
        
        df.to_csv(self.path_output_5_final, sep='\t', index=False)
        WriteLog("Step_04", "Sequence removal and renumbering completed successfully.")
                 
####################################################################################
# Run Script
####################################################################################

if __name__ == '__main__':
    create_folder(os.path.dirname(os.path.abspath(__file__)))

    mirror = MIrRORv2_Bacteria()
    mirror.read_gff()
    mirror.extract_sequences()
    mirror.reverse_and_complement()
    mirror.append_gtdb_taxonomy()
    mirror.curation_redundant_sequence()
    mirror.curation_ambiguous_nucleotides()
    mirror.curation_sequence_removal()

    print("\nOperon Extraction and Curation Complete \n")
    
