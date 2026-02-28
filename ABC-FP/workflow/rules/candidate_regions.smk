rule generate_chrom_sizes_bed_file:
	input:
		chrom_sizes = config['ref']['chrom_sizes']
	output:
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}
		"""

## sort narrowPeaks
rule sort_narrowpeaks:
	input:
		narrowPeak = get_peak_files,
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		chrom_sizes = config['ref']['chrom_sizes']
	conda:
		"../envs/abcenv.yml"
	output:
		narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		# intersect to remove alternate chromosomes
		bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
		"""


rule make_candidate_regions:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
		accessibility = get_accessibility_files,
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed'),
	params:
		TSS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'TSS'],
		chrom_sizes = config['ref']['chrom_sizes'],
		regions_blocklist = config['ref']['regions_blocklist'],
		peakExtendFromSummit = config['params_candidate']['peakExtendFromSummit'],
		nStrongestPeak = config['params_candidate']['nStrongestPeaks'],
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Peaks"),
		scripts_dir = SCRIPTS_DIR,
	conda:
		"../envs/abcenv.yml"
	output: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
	resources:
		mem_mb=determine_mem_mb
	shell: 
		"""
		python {params.scripts_dir}/makeCandidateRegions.py \
			--narrowPeak {input.narrowPeak}\
			--accessibility {input.accessibility} \
			--outDir {params.output_dir} \
			--chrom_sizes {params.chrom_sizes} \
			--chrom_sizes_bed {input.chrom_sizes_bed} \
			--regions_blocklist {params.regions_blocklist} \
			--regions_includelist {params.TSS} \
			--peakExtendFromSummit {params.peakExtendFromSummit} \
			--nStrongestPeak {params.nStrongestPeak}
		"""
