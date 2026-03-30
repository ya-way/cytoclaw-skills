---
name: openclaw-fastqc
description: Run FastQC quality control for FASTQ files, interpret quality metrics, and propose next actions (trim, recheck, or resequence). Use when users mention FastQC, FASTQ QC, adapter contamination, per-base quality, or NGS pre-alignment quality checks.
---

# openclaw-fastqc

## Inputs

- FASTQ file path(s)
- output directory
- optional thread count

## Command

```bash
mkdir -p <out_dir>
fastqc -o <out_dir> -t <threads> <fastq_1> [<fastq_2> ...]
```

Optional aggregation:

```bash
multiqc <out_dir> -o <out_dir>/multiqc
```

## Output Contract

Return:

- pass/conditional pass/fail
- top 1-3 quality risks
- actionable next steps
- report artifact paths
