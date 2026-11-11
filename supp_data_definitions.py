FAM='the ten families RDase, HyDase, GSTDase, OxDase, DehyDase, MetDase, FlaHase, VanHase, NHFeHase, SAMHase'
C={}  # (file, sheet) -> {column: definition}
common_meta={
'sample_id':'BioSample accession that identifies the sample; the runs of one BioSample are pooled.',
'host_system':'Host-associated system of the sample: Chicken (chicken gut), Coral, Rhizosphere or Sponge.',
'study':'Study grouping used for the analyses (one or more BioProjects from the same publication or project).',
'BioProject':'NCBI/ENA BioProject accession of the study.',
'GE':'Genome equivalents of the sample, estimated from the reads with MicrobeCensus; the unit that abundances are divided by (RPKG).',
'platform':'Sequencing platform of the runs (instrument_platform in the ENA record).',
}
C[(1,'samples')]={**common_meta,
'runs':'Accessions of the sequencing runs pooled in the sample.',
'n_runs':'Number of sequencing runs pooled in the sample.',
'reads_processed':'Reads (fragments) processed by Salmon when mapping to the target genes, summed over the runs of the sample.',
'Gb':'Sequencing yield in gigabases (billions of bases), from the base count in the ENA record.',
'lat':'Latitude of the sampling location in decimal degrees (blank where the sample record gives none).',
'lon':'Longitude of the sampling location in decimal degrees (blank where the sample record gives none).',
'country':'Country or sea region named in the sample record (blank where none is named).',
'reads_on_target_genes':'Reads assigned to the target genes of any of the ten families (the sum of the ten family read counts in Supplementary Data 4, sheet family_reads).'}
C[(1,'runs')]={**common_meta,
'run':'Sequencing run accession (SRA/ENA/DDBJ).',
'sample_accession':'BioSample accession of the run, as given in the ENA record.',
'scientific_name':'Organism/taxon label of the sample in the ENA record (for example "sponge metagenome").',
'host':'Host organism named in the sample record (blank where none is given).',
'isolation_source':'Isolation source named in the sample record.',
'library_strategy':'Library strategy label in the ENA record (WGS = whole-genome shotgun).',
'library_source':'Library source label in the ENA record (blank where the record gives none).',
'instrument_platform':'Sequencing platform in the ENA record.',
'reads_processed':'Reads (fragments) processed by Salmon for this run.',
'GE':'Genome equivalents of this run (MicrobeCensus).',
'fastp_reads_in':'Reads before quality filtering with fastp (blank where fastp statistics were not kept).',
'fastp_reads_out':'Reads remaining after fastp quality filtering.',
'q30_after':'Fraction of bases with Phred quality of at least 30 after filtering.',
'collection_date':'Collection date in the sample record.',
'study_title':'Title of the BioProject, as published in the record.'}
C[(3,'genes')]={
'hit_id':'Identifier of the curated gene (genome identifier plus contig and gene number).',
'genome_id':'Identifier of the genome (MAG or isolate genome) that carries the gene.',
'host_system':'Host-associated system of the genome catalogue the genome comes from: Chicken, Coral, Rhizosphere or Sponge.',
'enzyme_family':'Enzyme family of the gene (one of '+FAM+').',
'GTDB_lineage':'GTDB taxonomic lineage of the genome (domain to species).',
'lineage_source':'Where the lineage comes from: genome classification = the genome\'s own GTDB classification; species representative = the lineage of the species-representative genome of its GTDB species.',
'identity_to_nearest_HaloCycDB_reference_pct':'Amino-acid identity (%) of the gene to its closest characterized reference enzyme in HaloCycDB, from DIAMOND (blank where no reference aligned).',
'query_coverage_pct':'Percentage of the gene\'s protein covered by the alignment to that reference.',
'DIAMOND_evalue':'DIAMOND E-value of the alignment to the closest reference.',
'aligned_to_reference':'TRUE if the gene aligned to a reference enzyme of its family at E <= 10; FALSE if it did not.',
'subfamily':'Subfamily label (unassigned = no label; no_motif_signal = family without subfamily resolution).',
'Pfam_architecture':'Pfam-A domains found in the protein (hmmscan with the gathering cut-off, i-Evalue <= 1e-5), in order along the sequence.',
'DeepTMHMM_topology':'DeepTMHMM prediction: SP = signal peptide, TM = transmembrane, GLOB = globular, SP+TM = signal peptide and transmembrane.',
'length_aa':'Length of the protein in amino acids.'}
fam_cols={f:'RPKG of the family: reads per kilobase of target gene per genome equivalent in the sample (all target genes of the family pooled).' for f in ['RDase','HyDase','GSTDase','OxDase','DehyDase','MetDase','FlaHase','VanHase','NHFeHase','SAMHase']}
C[(4,'family_RPKG')]={**common_meta,'reads_processed':'Reads (fragments) processed by Salmon, summed over the runs of the sample.',**fam_cols}
C[(4,'family_reads')]={**common_meta,'reads_processed':'Reads (fragments) processed by Salmon, summed over the runs of the sample.',**{f:'Number of reads assigned to the target genes of the family.' for f in fam_cols}}
C[(4,'subfamily_RPKG')]={'sample_id':common_meta['sample_id'],'family':'Enzyme family.','subfamily':'Subfamily label.','RPKG':'Reads per kilobase of target gene per genome equivalent for the subfamily in the sample.','reads':'Number of reads assigned to the subfamily in the sample.'}
C[(5,'metatranscriptomes')]={
'run_accession':'Run accession of the metatranscriptome.','host_system':'Host-associated system of the sample.','study':common_meta['study'],
'origin':'earlier = run from the first selection (147 runs); added = run added later (43 runs); both were processed in the same way.',
**{f:'RPKM of the family: reads per kilobase of target gene per million reads processed (first 30 million reads of the run).' for f in fam_cols},
'million':'Millions of reads processed (at most 30).'}
nullc={'null_rpk_sub':'Reads per kilobase of the null genes (genes without a halogenase or dehalogenase function, drawn from genomes of the same system) in the read subsample (see nproc_sub) that was mapped to the null genes.',
'nproc_sub':'Number of reads processed in the read subsample (about 10 million read pairs) that was mapped to the null genes; nproc_full divided by nproc_sub is the factor that scales the null-gene counts to full depth.',
'nproc_full':'Number of reads processed when the run was mapped to the target genes at full depth.',
'null_RPKG':'Null-gene abundance of the run on the same scale as the family RPKG: null_rpk_sub x (nproc_full / nproc_sub) / GE. Used to normalize family abundance (null-normalized abundance).'}
C[(5,'surroundings')]={'run':'Run accession of the surrounding-environment metagenome.',
'system':'Host system the surroundings are compared with (Chicken, Coral, Rhizosphere or Sponge); a metagenome used for two systems has two rows.',
'group':'ambient_external = surrounding metagenome from a project outside the study; ambient_matched = surrounding metagenome from the same project as host samples.',
'BioProject':common_meta['BioProject'],'GE':common_meta['GE'],**fam_cols,**nullc,
'null_source':'Internal label of the batch in which the null-gene mapping was run (earlier for all rows).'}
C[(5,'host_null')]={'run':'Run accession of the host metagenome.','system':'Host system of the run.','study':common_meta['study'],'GE':common_meta['GE'],**nullc}
st={'safe_id':'Identifier of the predicted structure (family, subfamily, host system and number; used as the file name of the model).',
'gene':'Enzyme family of the structure.','subfamily':'Subfamily label.','environment':'Host system of the gene: Chicken, Coral, Rhizosphere or Sponge.',
'original_hit_id':'Identifier of the curated gene in Supplementary Data 3 (hit_id).',
'best_tm_score':'TM-score (0 to 1; above 0.5 = same fold) of the best Foldseek match in the searched databases.',
'percent_identity':'Sequence identity (%) of the structural alignment to that best match.',
'best_database':'Database that held the best match: pdb100 (Protein Data Bank), afdb50, afdb-swissprot, afdb-proteome (AlphaFold database sets) or mgnify_esm30 (ESM Metagenomic Atlas).',
'best_target':'Identifier of the best Foldseek match.',
'candidate_novel_fold':'TRUE if the best TM-score indicates no confident fold match (none of the 529 structures).',
'high_divergence_conserved_fold':'TRUE if the fold is conserved (best TM-score >= 0.9) but sequence identity to the match is below 40%.',
'plddt':'Mean pLDDT (0 to 100) of the predicted structure: confidence of the model.','n_residues':'Number of residues (chain length).',
'n_pockets':'Number of pockets detected by fpocket.','volume':'fpocket pocket volume (cubic angstrom).','druggability':'fpocket druggability score (0 to 1; higher = more druggable).',
'hydrophobicity':'fpocket hydrophobicity score of the pocket.','polarity':'fpocket polarity score of the pocket.',
'pfam_architecture':'Pfam-A domains of the protein, in order.','n_distinct_pfam':'Number of different Pfam domains in the protein.'}
C[(6,'structures_529')]={k:st[k] for k in ['safe_id','gene','subfamily','environment','original_hit_id','best_tm_score','percent_identity','best_database','candidate_novel_fold','high_divergence_conserved_fold','plddt','n_residues']}
C[(6,'four_groups_plddt_pockets')]={'safe_id':st['safe_id'],'gene':st['gene'],'group':'Subfamily group of the structure set (family__subfamily).','environment':st['environment'],'n_residues':st['n_residues'],
'mean_plddt':'Mean pLDDT over the residues.','median_plddt':'Median pLDDT over the residues.','min_plddt':'Lowest per-residue pLDDT.','max_plddt':'Highest per-residue pLDDT.',
'n_pockets_detected':st['n_pockets'],'druggability_score':st['druggability'],'pocket_volume':st['volume'],'hydrophobicity_score':st['hydrophobicity'],'polarity_score':st['polarity']}
for sh in ['GSTDase_plddt_pockets','MetDase_NHFeHase_plddt_pockets']:
    C[(6,sh)]={k:st[k] for k in ['safe_id','gene','subfamily','environment','plddt','n_pockets','volume','druggability','hydrophobicity','polarity','n_residues']}
for sh in ['fold_search_flagship','fold_search_indicator','fold_search_GSTDase','fold_search_MetDase_NHFeHase']:
    C[(6,sh)]={k:st[k] for k in ['safe_id','gene','subfamily','environment','original_hit_id','best_database','best_target','best_tm_score','percent_identity','candidate_novel_fold','high_divergence_conserved_fold']}
C[(6,'pfam_domains')]={k:st[k] for k in ['safe_id','gene','subfamily','environment','n_residues','pfam_architecture','n_distinct_pfam']}
C[(6,'p2rank_pockets')]={'safe_id':st['safe_id'],'Structure set':'Structure set of the model, named by its enzyme families (RDase, HyDase, FlaHase, VanHase / GSTDase / MetDase, NHFeHase).'}
bg={'genome_id':'Identifier of the genome analyzed with antiSMASH.',
'selected_via':'Why the genome was analyzed: flagship = carries a flagship subfamily; hub = belongs to a hub genus; hub+flagship = both.',
'contig':'Contig on which the biosynthetic gene cluster (BGC) lies.','bgc_start':'Start position of the BGC region on the contig (base pairs).','bgc_end':'End position of the BGC region on the contig (base pairs).',
'product_class':'antiSMASH product class of the BGC (for example NRPS, terpene, RiPP).',
'nearest_target_hit':'Curated gene (hit_id) closest to the BGC on the same contig.',
'relationship':'within = the gene lies inside the BGC region; flank_<=10kb = the gene lies within 10 kb outside the region; distant = farther than 10 kb; "-" = no curated gene on the contig.',
'distance_bp':'Distance in base pairs between the BGC region and the nearest curated gene (0 where within).'}
C[(7,'primary_ge3')]=bg; C[(7,'sensitivity_ge2')]=bg
