This readme file and repository was created by Daur Meretukov, MRes in Cancer Informatics student, Imperial College of London, 2024

MRes Cancer informatics, Project 1

HCT116 immunopeptidomics project at Institute of Cancer Research, Functional Proteomics team


Supervisors: James Wright| Sreejan Bandyopadhyay | Jyoti Choudhary

Files description:

Folder 'bash_scripts' contain 1 .txt file that covers all data manipulation and analysis files required to reproduce project results, including plots and tables. 

  1) hct116_project_bash_scripts.txt - all bash scripts used to run jobs for Trinity, BLAST, Kallisto, Immunopeptidomics analysis algorithm on SLURM-based high-performance computing cluster

Folder 'python_project' contain 3 .py files that covers all data manipulation and analysis files required to reproduce project results, including plots and tables. 

  1) HCT116_project_part_1_assembly_and_blast_core_database.py - Trinity FASTA files and BLAST output data manipulation, resulting in core database for further analysis/manipulation
  2) HCT116_project_part_2_RiboSeq_Translation_to_AA.py - translation of transcripts nucleic acid sequences into amino acid sequences with three framed ORF split
  3) HCT116_project_part_3_Immunopeptidome_and_plots.py - RNAseq, RiboSeq, immunopeptidome analysis + plots
